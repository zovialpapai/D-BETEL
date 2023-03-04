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

#  Functions ----------------------------------------------------------

LEL_MCM <- function(y, X, Beta, v){ 

N = nrow(X); p = ncol(X)

# Variable Declaration
U <- Variable(rows = N,cols = 1)

# Objective Function
Objective <- - sum(CVXR::entr(U))
  
# Moment conditions
mcm <- function(U, y, X, Beta){
M = diag(U)%*%(y - exp(X%*%Beta))
ret = t(t(M)%*%X)
return(ret)
}
  
  
# Constraints
Constraints_mcm <- list(sum(U) == 1 ,
                        U >=  10^(-20),
                        mcm(U, y, X, Beta)== rep(v, p))
  
# Setting up the Problem
Prob_mcm <- Problem(Minimize(Objective), Constraints_mcm)
#cat("Prob_WS is DCP:", is_dcp(Prob_WS), "\n")
  
# Solving the Problem
Result_mcm <- CVXR::psolve(Prob_mcm) 
if(Result_mcm$status != "solver_error"){
# Optimal Value
##Optimal = -Result_WS$value
Optimal = sum(log(Result_mcm$getValue(U)))
O = Optimal
# Finding P
P = Sol = Result_mcm$getValue(U)
}else{
O = -10^10
# Finding P
P = Sol = rep(NA, N)
}
  
return(list("Optimal_Value" = O,"P" = P))
}

nrep = 50; # m = 1
Result<- foreach(m = 1:nrep, .combine = rbind)%dopar%{ 
res<- tryCatch({
# Generate Data -----------------------------------------------------------
set.seed(3076+(m-1))
# Sample Size
N = 250 # MAKE IT 250
# Purterbations
perturbation = vector(length = N) 
for(k in 1:N){
if(runif(1) < 0.9){
perturbation[k] = 0
}else{
perturbation[k] = rnorm(1,1,0.1) 
}
}
# Generating from Simple Linear Regression
Beta_0 = 5; Beta_1 = 1; 
z = rnorm(N,5,1); X = cbind(rep(1,N), z)
y_temp = exp(Beta_0*rep(1,N) + Beta_1*z + perturbation )
y = vector(length = N)
for(k in 1:N){
y[k] = rpois(1, lambda = y_temp[k])
}

# Testing -----------------------------------------------------------------
#Test_Out = LEL_MCM(y, X, Beta = c(5, 1), v = 0); Test_Out$Optimal_Value; Test_Out$P

#fit <- glm(y ~ z, family = poisson(link=log))
#Test_Out = LEL_MCM(y, X, Beta = fit$coefficients, v = 0); Test_Out$Optimal_Value; Test_Out$P
# MLE ---------------------------------------------------------------------
    fit <- glm(y ~ z, family = poisson(link=log)); out = summary(fit) 
    # Asymtotics of MLE
    Beta_MLE = Beta_Init =  fit$coefficients; v_Init = 0
    #Beta_Init[1] =  runif(1, Beta_0*0.99, Beta_0*1.01); Beta_Init[2] =  runif(1, Beta_1*0.99, Beta_1*1.01)
    Asymtotic_CovBeta = diag(out$coefficients[,2])
    
    # Storage and control
    D = 1/20# Scale up Proposal Variance
    nIter = 1000; burn_in = floor(nIter/2)
    
    Beta_Storage = matrix(NA, nrow = nIter, ncol = 2)
    Accepted = 0
# Sampler -----------------------------------------------------------------
    for(iter in 1:nIter){ 
      #  proposal
      Beta_Prop = mvrnorm(n = 1, mu = Beta_Init , Sigma = D*Asymtotic_CovBeta)
      v_Prop = rnorm(1, v_Init, 1)
      # Likelihood
      llTemp_top = LEL_MCM(y, X, Beta = Beta_Prop, v = 0)
      llTemp_bot = LEL_MCM(y, X, Beta = c(Beta_Init), v = 0)
      
      ll_top = llTemp_top$Optimal_Value
      ll_bot = llTemp_bot$Optimal_Value
      
      # Prior
      lPrior_top = dnorm(Beta_Prop[1], mean = 0, sd = 100, log = TRUE) +
                   dnorm(Beta_Prop[2], mean = 0, sd = 100, log = TRUE) +
                   dt(v_Prop, df = 2, ncp = 0, log = TRUE)
      lPrior_bot = dnorm(Beta_Init[1], mean = 0, sd = 100, log = TRUE) + 
                   dnorm(Beta_Init[2], mean = 0, sd = 100, log = TRUE) +
                   dt(v_Init, df = 2, ncp = 0, log = TRUE)
                   
      # Proposal
      lProp_top = dmvnorm(Beta_Init, mean = Beta_Prop , sigma = D*Asymtotic_CovBeta, log = TRUE)            
      if(lProp_top == -Inf){lProp_top = -10^10}
      lProp_bot = dmvnorm(Beta_Prop, mean = Beta_Init , sigma = D*Asymtotic_CovBeta, log = TRUE)            
      if(lProp_bot == -Inf){lProp_bot = -10^10}
      
      # M--H ratio
      top = ll_top + lPrior_top + lProp_top; print(top)
      bot = ll_bot + lPrior_bot + lProp_bot
      p_accept = min(1,exp(top-bot))
      #print(c(top, bot, p_accept))
      if(is.na(p_accept)){p_accept = 0}
      if(runif(1) < p_accept){
        Beta_Storage[iter, ] =  Beta_Init = Beta_Prop; 
        v_Init = v_Prop
        #loglik_mat[iter, ] = llTemp_top$P
        Accepted = Accepted + 1
      }
      Beta_Storage[iter,] = Beta_Init
      print(iter)
    }
    
    c(colMeans(Beta_Storage[((burn_in +1):nIter),]), 
      abs(colMeans(Beta_Storage[((burn_in +1):nIter),])- c(Beta_0, Beta_1)),
      2*1.96*apply(Beta_Storage[((burn_in +1):nIter),], 2, sd),
      as.numeric(abs(colMeans(Beta_Storage[((burn_in +1):nIter),])- c(Beta_0, Beta_1)) < 1.96*apply(Beta_Storage[((burn_in +1):nIter),], 2, sd)),
      (Accepted/nIter)*100)
}, error = function(e){rep(NA,9)}) 
}

Result
colMeans(Result, na.rm = T)
