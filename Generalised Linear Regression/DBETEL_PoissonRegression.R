# Dependencies ------------------------------------------------------------
library(truncnorm)
library(MASS)
library(mvtnorm)

library(MCMCpack)
library(psych)
library(CVXR)
library(bdsmatrix)

library(expm)

library("loo")

library(doParallel)
library(foreach)

numCores <- detectCores()
numCores
registerDoParallel(numCores)

# Controls ----------------------------------------------------------------

# Optimization Controls (May change at each Instance) 
epsilon = 3 # Change to 9
lambda = 10^5

# Simulations Controls
nIter = 2000 # MAKE IT 2000
nrep = 50

# ELWS Functions ----------------------------------------------------------

LEL_WS <- function(X_Dat, mu1,sigma1, epsilon,lambda){ 
  X = X_Dat; N = nrow(X)
  # Calculation of M
  M = matrix(NA, nrow = 1, ncol = N)
  M = ((sweep(X, 2, mu1))^2) %*% c(1,1)
  # Optimization via DCP (Using Embedded Conic Solver)
  # Variable Declaration
  U <- Variable(rows = N,cols = 1)
  #U <- rep(1/N,N) # REMOVE
  # Objective Function
  ##Objective <- - sum(log(U))
  Objective <- - sum(CVXR::entr(U))
  
  # The RMW Calculation via CVXR functions
  #### Discrete approximation
  Univariate_term <- function(U, mu1, sigma1){
    # Datatype adjustment
    w = value(U)
    ## Numerical adjustments
    w[is.na(w)] = 10^(-5)
    
    M1 = cbind(X[ , 1], w); M2 = cbind(X[ , 2], w) 
    M1 = M1[order(M1[,1],decreasing=FALSE),];  M2 = M2[order(M2[,1],decreasing=FALSE),]
    M1 = cbind(M1, cumsum(M1[ ,2])); M2 = cbind(M2, cumsum(M2[ ,2]))
    
    ## Numerical adjustments
    M1[which(M1[ , 3] < 10^(-5) ), 3 ]  = 10^(-5); M2[which(M2[ , 3] < 10^(-5) ), 3 ]  = 10^(-5)
    M1[which(M1[ , 3] > 1 - 10^(-5) ), 3 ]  = 1- 10^(-5); M2[which(M2[ , 3] > 1 -  10^(-5) ), 3 ]  = 1 - 10^(-5)
    
    # Appx integral by sum
    m1 = sum((M1[ ,1] - qnorm(M1[ ,3], mu1[1], sigma1[1, 1])^2)*M1[ ,2])
    m2 = sum((M2[ ,1] - qnorm(M2[ ,3], mu1[2], sigma1[2, 2])^2)*M2[ ,2])
    ret = m1 + m2
    return(ret)
  }
  
  ####
  
  RMW_CVXR <-function(U, mu1, sigma1){
    Negative_Entropy = - sum( CVXR::entr(U) )
    OT = sum(multiply(M,U))
    Trace_Term = (psych::tr(sigma1)) 
    return((Negative_Entropy/lambda) + OT + Trace_Term + Univariate_term(U, mu1, sigma1))
  }
  
  
  # Constraints
  Constraints_WS <- list(sum(U) == 1 ,
                         U >=  10^(-20),
                         RMW_CVXR(U, mu1, sigma1) <= epsilon
  )
  
  # Setting up the Problem
  Prob_WS <- Problem(Minimize(Objective), Constraints_WS)
  #cat("Prob_WS is DCP:", is_dcp(Prob_WS), "\n")
  
  # Solving the Problem
  Result_WS <- CVXR::psolve(Prob_WS) 
  # Optimal Value
  ##Optimal = -Result_WS$value
  Optimal = sum(log(Result_WS$getValue(U)))
  O = Optimal
  # Finding P
  P = Sol = Result_WS$getValue(U)
  
  return(list("Optimal_Value" = O,"P" = P))
}

Test_Out = LEL_WS(X_Dat = ((t(sqrtm(solve(I))%*%t(Scores(beta_MLE))))), 
                  mu1 = c(0,0), sigma1 = (sqrtm(solve(I)))%*%I%*%(sqrtm(solve(I))) , 
                  epsilon = 10 ,lambda)

Test_Out$Optimal_Value

Test_Out$P

#Score_Now = Scores(beta_MLE)
#var(Score_Now)
#sqrtm(solve(I))%*%Score_Now[5,]
#(t(sqrtm(solve(I))%*%t(Scores(beta_MLE))))
#((t(sqrtm(solve(var(Score_Now)))%*%t(Scores(beta_MLE)))))
#Test_Out = LEL_WS(X_Dat = ((t(sqrtm(solve(I))%*%t(Scores(beta_MLE))))), mu1 = c(0,0), sigma1 = (sqrtm(solve(I)))%*%I%*%(sqrtm(solve(I))) , epsilon = 25 ,lambda)
#Test_Out$Optimal_Value


Result<- foreach(m = 1:nrep, .combine = rbind)%dopar%{ 
  res<- tryCatch({
    # Generate Data -----------------------------------------------------------
    set.seed(3076+(m-1))
    # Sample Size
    N = 100# MAKE IT 250
    # Purterbations
    perturbation = vector(length = N) 
    for(k in 1:N){
      if(runif(1) < 0.90){
        perturbation[k] = 0
      }else{
        perturbation[k] = rnorm(1,1,0.1) 
      }
    }
    # Generating from Simple Linear Regression
    beta_0 = 5; beta_1 = 1; z = rnorm(N,5,1)
    y_temp = exp(beta_0*rep(1,N) + beta_1*z + perturbation )
    y = vector(length = N)
    for(k in 1:N){
      y[k] = rpois(1, lambda = y_temp[k])
    }
    #y
    
    # MLE ---------------------------------------------------------------------
    fit <- glm(y ~ z, family = poisson(link=log))
    #summary(fit) 
    beta_MLE = fit$coefficients
    # Asymtotics of MLE
    X = cbind(rep(1,N),z)
    V = matrix(0, nrow = N, ncol = N)
    for(k in 1:N){
      V[k,k] = exp(t(X[k,])%*%beta_MLE)
    }
    I = t(X)%*%V%*%X # Covariance of Score Function
    Asymtotic_CovBeta = solve(I) # Asymtotic Covariance of beta_MLE
    
    # Multivariate Delta Method
    A1 = log(sum(diag(V)))
    A2 = log(sum(diag(V)*(X[ ,2]^2)))
    A3 = sum(diag(V)*X[,2])/exp((A1+A2)/2)
    
    g_prime = matrix(NA, nrow = 5, ncol = 2 )
    g_prime[1,1] = 1; g_prime[1,2] = sum(diag(V)*X[ ,2])/sum(diag(V))
    g_prime[2,1] = 1; g_prime[2,2] = sum(diag(V)*(X[ ,2]^3))/sum(diag(V)*(X[ ,2]^2))
    g_prime[3,1] = 0
    g_prime[3,2] = (sum(diag(V)*(X[,2])^2) - sum(diag(V)*X[,2])*0.5*(g_prime[1,2] + g_prime[2,2] ))/exp((A1+A2)/2) 
    g_prime[4, ] = c(1,0); g_prime[5, ] = c(0,1)
    
    AsymtoticCov_g = g_prime%*%Asymtotic_CovBeta%*%t(g_prime)
    
    Prop = mvrnorm(n = 1, mu = c(A1,A2,A3,beta_MLE) , Sigma = AsymtoticCov_g)
    Prop[3] = min(Prop[3],1)
    Beta_Prop = Prop[4:5]
    Sigma_Prop = rbind(c(exp(Prop[1]),Prop[3]*exp((Prop[1]+Prop[2])/2)),c(Prop[3]*exp((Prop[1]+Prop[2])/2),exp(Prop[2])))
    
    # Function Calculate Score Function at beta
    Scores <-function(beta){
      Score = matrix(NA, nrow = N, ncol = 2)
      Score[,1] = y - exp(X%*%beta)
      Score[,2] = X[,2]*(y - exp(X%*%beta))
      return(Score)
    }
    #Scores(beta_MLE)
    #colSums(Scores(beta_MLE))
    
    
    
    
    # Sampler -----------------------------------------------------------------
    loglik_mat = matrix(NA, nrow = nIter , ncol = N)
    # Storage
    D = 10 # Scale up Proposal Variance
    Accepted = 0
    burn_in = floor(nIter/2)
    Beta_Storage = matrix(NA, nrow = nIter, ncol = 2)
    # Starting point at MLE
    Beta_Init = beta_MLE; Sigma_Init = I
    # Preparation for Proposal
    V_Init_Vec = exp(X%*%Beta_Init)
    A1_Init = log(sum(V_Init_Vec))
    A2_Init = log(sum(V_Init_Vec*(X[ ,2]^2)))
    A3_Init = sum(V_Init_Vec*X[,2])/exp((A1_Init+A2_Init)/2)
    
    for(iter in 1:nIter){ 
      # Proposal
      Prop = mvrnorm(n = 1, mu = c(A1_Init,A2_Init,A3_Init,Beta_Init) , Sigma = D*AsymtoticCov_g)
      Prop[3] = min(Prop[3],1)
      Beta_Prop = Prop[4:5]
      Sigma_Prop = rbind(c(exp(Prop[1]),Prop[3]*exp((Prop[1]+Prop[2])/2)),c(Prop[3]*exp((Prop[1]+Prop[2])/2),exp(Prop[2])))
      # MH Ratio
      llTemp_top = LEL_WS(X_Dat = ((t(sqrtm(solve(I))%*%t(Scores(Beta_Prop))))), mu1 = c(0,0), sigma1 = (sqrtm(solve(I)))%*%Sigma_Prop%*%(sqrtm(solve(I))) , epsilon  ,lambda)
      llTemp_bot = LEL_WS(X_Dat = ((t(sqrtm(solve(I))%*%t(Scores(Beta_Init))))), mu1 = c(0,0), sigma1 = (sqrtm(solve(I)))%*%Sigma_Init%*%(sqrtm(solve(I))) , epsilon  ,lambda)
      
      loglik_mat[iter, ] = llTemp_bot$P
      
      ll_top = llTemp_top$Optimal_Value
      ll_bot = llTemp_bot$Optimal_Value
      
      lPrior_top = dnorm(Prop[1], mean = 0, sd = 100, log = TRUE) +
        dnorm(Prop[2], mean = 0, sd = 100, log = TRUE) + 
        dnorm(Beta_Prop[1], mean = 0, sd = 100, log = TRUE) +
        dnorm(Beta_Prop[2], mean = 0, sd = 100, log = TRUE)
      lPrior_bot = dnorm(A1_Init, mean = 0, sd = 100, log = TRUE) +
        dnorm(A2_Init, mean = 0, sd = 100, log = TRUE) + 
        dnorm(Beta_Init[1], mean = 0, sd = 100, log = TRUE) +
        dnorm(Beta_Init[2], mean = 0, sd = 100, log = TRUE)
      
      lProp_top = dmvnorm(c(A1_Init, A2_Init,A3_Init,Beta_Init), mean = Prop , sigma = D*AsymtoticCov_g, log = TRUE)            
      if(lProp_top == -Inf){lProp_top = -10^10}
      lProp_bot = dmvnorm(Prop , mean = c(A1_Init, A2_Init,A3_Init,Beta_Init) , sigma = D*AsymtoticCov_g, log = TRUE)            
      if(lProp_bot == -Inf){lProp_bot = -10^10}
      
      top = ll_top + lPrior_top + lProp_top
      bot = ll_bot + lPrior_bot + lProp_bot
      p_accept = min(1,exp(top-bot))
      if(runif(1) < p_accept){
        A1_Init = Prop[1]; A2_Init = Prop[2]; A3_Init = Prop[3]
        Beta_Init = Beta_Prop; Sigma_Init = Sigma_Prop
        loglik_mat[iter, ] = llTemp_top$P
        Accepted = Accepted + 1
      }
      Beta_Storage[iter,] = Beta_Init
      #print(iter)
    }
    
    MS1 = loo(log(loglik_mat[(nIter/2):(nIter-1), ]))
    #Result[m, ] = c(beta_MLE,colMeans(Beta_Storage[((burn_in +1):nIter),]),MS1$elpd_loo,MS1$se_elpd_loo,(Accepted/nIter)*100)
    
    #c(beta_MLE, (summary(fit)$coefficients[ ,2])*sqrt(N),
    #  colMeans(Beta_Storage[((burn_in +1):nIter),]), apply(Beta_Storage[((burn_in +1):nIter),], 2, sd),
    #  MS1$elpd_loo,MS1$se_elpd_loo,
    #  (Accepted/nIter)*100)
    
    c(beta_MLE, 
      abs(beta_MLE - c(beta_0, beta_1)),
      2*1.96*(summary(fit)$coefficients[ ,2])*sqrt(N),
      as.numeric(abs(beta_MLE - c(beta_0, beta_1)) < 1.96*(summary(fit)$coefficients[ ,2])*sqrt(N)),
      colMeans(Beta_Storage[((burn_in +1):nIter),]), 
      abs(colMeans(Beta_Storage[((burn_in +1):nIter),])- c(beta_0, beta_1)),
      2*1.96*apply(Beta_Storage[((burn_in +1):nIter),], 2, sd),
      as.numeric(abs(colMeans(Beta_Storage[((burn_in +1):nIter),])- c(beta_0, beta_1)) < 1.96*apply(Beta_Storage[((burn_in +1):nIter),], 2, sd)),
      MS1$elpd_loo,MS1$se_elpd_loo,
      (Accepted/nIter)*100)
  }, error = function(e){rep(NA,17)}) # was 11
}

Result
colMeans(Result, na.rm = T)
#saveRDS(colMeans(Result, na.rm = T), paste0("S",N,"e",epsilon))

