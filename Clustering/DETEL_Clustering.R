# Dependency --------------------------------------------------------------
library("sn")
library("KScorrect")
library("loo")

library(truncnorm)
library(MASS)
library(mvtnorm)

library(MCMCpack)
library(psych)
library(CVXR)
library(bdsmatrix)

library(expm)
library("foreach")
library("doParallel")

# Controls ---------------------------------------------------------
# Data genration controls
N = 500; d = 2
Alpha = rep(5, d)
# Alpha = rep(2, d)
# N_grid = c(100, 200, 300, 500)
# BF = matrix(NA, nrow = length(N_grid), ncol = 3)
# colnames(BF) <- c("Sample size", "M1", "M2")


#### Simulation control
epsilon = 2.25
# epsilon_grid = c(2.25, 3, 2^2, 2^3, 2^4)
nrep = 1
nIter = 100

# M1: Functions ---------------------------------------------------------------
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
    if(is.null(Result_WS$value)){
      O = -10^10
      P = rep(NA, N)
    }else{
      Optimal = sum(log(Result_WS$getValue(U)))
      O = Optimal
      # Finding P
      P = Sol = Result_WS$getValue(U)
    }
    
    return(list("Optimal_Value" = O,"P" = P))
}
# Test_Out = LEL_WS(X_Dat = x, 
#                   mu1 = rep(0, d), sigma1 = diag(1, d) , 
#                   epsilon = 3,lambda = 100)
# Test_Out$Optimal_Value
#Test_Out$P


# M2 Functions ------------------------------------------------------------
LEL_WS_d2 <- function(X_Dat, mu1,sigma1, mu2, sigma2, w, epsilon,lambda){ 
  X = X_Dat; N = nrow(X)
  rowsum_vec = matrix(1, nrow = 2, ncol = 1)
  # Calculation of M
  M = matrix(NA, nrow = 2, ncol = N)
  M[1, ] = ((sweep(X, 2, mu1))^2) %*% c(1,1)
  M[2, ] = ((sweep(X, 2, mu2))^2) %*% c(1,1)
  # Optimization via DCP (Using Embedded Conic Solver)
  # Variable Declaration
  # U <- matrix(1/(2*N), nrow = N, ncol = 2)
  U <- Variable(rows = N, cols = 2)
  U_obj <- U%*%rowsum_vec
  # Objective Function
  ##Objective <- - sum(log(U))
  Objective <- - sum(CVXR::entr(U_obj))
  
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
  
  RMW_CVXR <-function(U, mu1, sigma1, mu2, sigma2, w){
    Negative_Entropy = - sum( CVXR::entr(U) )
    #OT = sum(multiply(M,U))
    OT = sum(M%*%U)
    # OT  = 0
    Trace_Term1 = (psych::tr(sigma1)); Trace_Term2 = (psych::tr(sigma2)) 
    ret = (Negative_Entropy/lambda) + OT + 
      w*Trace_Term1 + (1-w)*Trace_Term2
    w*Univariate_term(U, mu1, sigma1) + (1-w)*Univariate_term(U, mu2, sigma2)
    return(ret)
  }
  
  
  # Constraints
  Constraints_WS <- list(sum(U_obj) == 1 ,
                         U_obj >=  10^(-20),
                         RMW_CVXR(U, mu1, sigma1, mu2, sigma2, w) <= epsilon
  )
  
  # Setting up the Problem
  Prob_WS <- Problem(Minimize(Objective), Constraints_WS)
  #cat("Prob_WS is DCP:", is_dcp(Prob_WS), "\n")
  
  # Solving the Problem
  Result_WS <- CVXR::psolve(Prob_WS) 
  # Optimal Value
  ##Optimal = -Result_WS$value
  if(is.null(Result_WS$value)){
    O = -10^10
    P = matrix(NA, nrow = N, ncol = 2)
  }else if(Result_WS$value==Inf){
    O = -10^10
    P = matrix(NA, nrow = N, ncol = 2)
  }else{
    Optimal = sum(log(Result_WS$getValue(U)))
    O = Optimal
    # Finding P
    P = Sol = Result_WS$getValue(U)
  }
  
  
  return(list("Optimal_Value" = O,"P" = rowSums(P)))
}
# Test_Out = LEL_WS_d2(X_Dat = x, 
#                      mu1 = rep(0, d), sigma1 = diag(1, d) ,
#                      mu2 = rep(0, d), sigma2 = diag(1, d) ,
#                      w = 0.5,
#                      epsilon = 5,lambda = 100)
# Test_Out$Optimal_Value
# Test_Out$P


# Repeated simulations ----------------------------------------------------
M1_varying_epsilon = M2_varying_epsilon = matrix(NA, 
                                                 nrow = length(epsilon_grid),
                                                 ncol = 4)
colnames(M1_varying_epsilon)<-colnames(M2_varying_epsilon)<-c("epsilon","Marg-lik", "ELPD_LOO", "SE(ELPD_LOO)")

for(epsilon in epsilon_grid){
M1_Results = M2_Results =  matrix(NA, nrow = nrep, ncol = 3)
colnames(M1_Results) <- colnames(M2_Results) <- c("Marg-lik", "ELPD_LOO", "SE(ELPD_LOO)")

for(m in 1:nrep){
# Generate Data -----------------------------------------------------------
set.seed(123+(m-1))
x = rmsn(n = N, xi = rep(0, d), Omega = diag(1, d), alpha = Alpha,  tau = 0, dp = NULL)

# M1: Simulations -------------------------------------------------------------
# index = 1
# result = matrix(NA, nrow = nrep, ncol = 10)
# rept = 1
# numCores <- detectCores()
# numCores
# registerDoParallel(numCores)
# Result<- foreach(m = 1:nrep, .combine = rbind)%dopar%{
# res<- tryCatch({
# MLE 
mu_MLE = colMeans(x); var_MLE = var(x)

# Sampler
# Storage
loglik_mat = matrix(NA, nrow = nIter , ncol = N)
Accepted = 0
burn_in = floor(nIter/2)
MH_storage = rep(NA, nIter)
mu_Storage = matrix(NA, nrow = nIter, ncol = d)
Sigma_Storage = array(NA, c(d, d, nIter))

# Starting point at MLE
mu_Init = mu_MLE; Sigma_Init = var_MLE
D = 1 # Scale up Proposal Variance

# Iterate the sampler
for(iter in 1:nIter){ 
# Proposal
mu_prop = mvrnorm(n = 1, mu = mu_Init, Sigma = D*(Sigma_Init/N))
# Sigma_prop = MCMCpack::rwish(v = 100, S = Sigma_Init/100)
Sigma_prop = MCMCpack::rwish(v = 100, S = var_MLE/100)
# Pre-calculations
llTemp_top = LEL_WS(X_Dat = x, 
                    mu1 = mu_prop, sigma1 = Sigma_prop , 
                    epsilon,lambda)
llTemp_bot = LEL_WS(X_Dat = x, 
                    mu1 = mu_Init, sigma1 = Sigma_Init , 
                    epsilon,lambda)
loglik_mat[iter, ] = llTemp_bot$P
# M-H Ratio  
# DETEL Likelihood
ll_top = llTemp_top$Optimal_Value
ll_bot = llTemp_bot$Optimal_Value
# # Parametrix Likelihood
# ll_top = sum(mvtnorm::dmvnorm(x, mean = mu_prop, sigma = Sigma_prop, log = TRUE))
# ll_bot = sum(mvtnorm::dmvnorm(x, mean = mu_Init, sigma = Sigma_Init, log = TRUE))
# Prior   
# Good example
# lPrior_top = mvtnorm::dmvnorm(mu_prop, mean = rep(0, d), sigma = 1000*diag(d), log = TRUE) + 
#              log(MCMCpack::dwish(W = Sigma_prop, v = 2, S = diag(1, d)*1000))
# lPrior_bot = mvtnorm::dmvnorm(mu_Init, mean = rep(0, d), sigma = 1000*diag(d), log = TRUE) + 
#              log(MCMCpack::dwish(W = Sigma_Init, v = 2, S = diag(1, d)*1000))

# Other example
lPrior_top = mvtnorm::dmvnorm(mu_prop, mean = rep(0, d), sigma = 1000*diag(d), log = TRUE) + 
  log(MCMCpack::dwish(W = Sigma_prop, v = 2, S = diag(1, d)*1000))
lPrior_bot = mvtnorm::dmvnorm(mu_Init, mean = rep(0, d), sigma = 1000*diag(d), log = TRUE) + 
  log(MCMCpack::dwish(W = Sigma_Init, v = 2, S = diag(1, d)*1000))
# Proposal
# mu_prop = mvrnorm(n = 1, mu = mu_Init, Sigma = D*(Sigma_Init/N))
# Sigma_prop = rWishart(n = 1, df = 2, Sigma = Sigma_Init)[ , ,1]         
lProp_top = mvtnorm::dmvnorm(mu_Init , mean = mu_prop , sigma = D*(Sigma_prop/N), log = TRUE) +
            log(MCMCpack::dwish(W = Sigma_Init, v = 100, S = Sigma_prop/100))
if(lProp_top == -Inf){lProp_top = -10^10}
lProp_bot = mvtnorm::dmvnorm(mu_prop , mean = mu_Init , sigma = D*(Sigma_Init/N), log = TRUE) +
            log(MCMCpack::dwish(W = Sigma_prop, v = 100, S = Sigma_Init/100))
if(lProp_bot == -Inf){lProp_bot = -10^10}
        
top = ll_top + lPrior_top + lProp_top; if(is.na(top)){top = -10^10}
bot = ll_bot + lPrior_bot + lProp_bot
p_accept = min(1,exp(top-bot))
if(is.na(p_accept)){p_accept = 0}
if(runif(1) < p_accept){
mu_Init = mu_prop; Sigma_Init = Sigma_prop
loglik_mat[iter, ] = llTemp_top$P
Accepted = Accepted + 1
}
mu_Storage[iter,] = mu_Init
Sigma_Storage[ , , iter] = Sigma_Init
MH_storage[iter] = sum(log(loglik_mat[iter, ] ))
print(iter)
}
# Marginal likelihood calculation
post_prob = exp(MH_storage - max(MH_storage))/sum(exp(MH_storage - max(MH_storage)))
M1_Marglik = MH_storage[which(post_prob == max(post_prob))] - log(max(post_prob))

# ELPD Calculation      
MS1 = loo::loo(log(loglik_mat[(nIter/2):(nIter-1), ]))
#Result[m, ] = c(beta_MLE,colMeans(Beta_Storage[((burn_in +1):nIter),]),MS1$elpd_loo,MS1$se_elpd_loo,(Accepted/nIter)*100)
 
#  output     
c(#mu_MLE, 
  #(summary(fit)$coefficients[ ,2])*sqrt(N),
  #colMeans(mu_Storage[((burn_in +1):nIter),]), 
  #apply(mu_Storage[((burn_in +1):nIter),], 2, sd),
  unique(M1_Marglik), 
  MS1$elpd_loo, MS1$se_elpd_loo, 
  (Accepted/nIter)*100
  )
# }, error = function(e){rep(NA, 4)})
# }


# M2: Simulations ---------------------------------------------------------
# nIter = 1000
# result = matrix(NA, nrow = nrep, ncol = 10)
# for(rept in 1:nrep){
# Generate from Skew Normal(0,1,\alpha) 
# set.seed(123)
#Alpha = rep(3, d)
# x = rmsn(n = N, xi = rep(0, d), Omega = diag(1, d), alpha = Alpha,  tau = 0, dp = NULL)
  
# Initialization
mu1_Init = mu2_Init = colMeans(x); Sigma1_Init = Sigma2_Init = var(x); wt_Init = 0.5
# Record Keeping
MH_storage = matrix(0, nrow = nIter, ncol = 1) # Dunson Miller fills it up with priors
loglik_mat = matrix(NA, nrow = nIter , ncol = N)
Accepted = 0; D = 1
  
# Draw Samples
for(iter in 1:nIter){ 
  # Proposal
  mu1_prop = mvrnorm(n = 1, mu = mu1_Init, Sigma = D*(Sigma1_Init/N))
  mu2_prop = mvrnorm(n = 1, mu = mu2_Init, Sigma = D*(Sigma2_Init/N))
  Sigma1_prop = MCMCpack::rwish(v = 100, S = Sigma1_Init/100)
  Sigma2_prop = MCMCpack::rwish(v = 100, S = Sigma2_Init/100)
  wt_prop = runif(1)
  # Pre-calculations
  llTemp_top = LEL_WS_d2(X_Dat = x, 
                         mu1 = mu1_prop, sigma1 = Sigma1_prop,
                         mu2 = mu2_prop, sigma2 = Sigma2_prop  ,
                         w = wt_prop,
                         epsilon, lambda)
  llTemp_bot = LEL_WS_d2(X_Dat = x, 
                         mu1 = mu1_Init, sigma1 = Sigma1_Init,
                         mu2 = mu2_Init, sigma2 = Sigma2_Init  ,
                         w = wt_Init,
                         epsilon, lambda)
  loglik_mat[iter, ] = llTemp_bot$P
  # M-H Ratio  
  # DETEL Likelihood
  ll_top = llTemp_top$Optimal_Value
  ll_bot = llTemp_bot$Optimal_Value
  # # Parametric Likelihood
  
  # Prior    
  # Good examples
  # lPrior_top = mvtnorm::dmvnorm(mu1_prop, mean = rep(0, d), sigma = 1000*diag(d), log = TRUE) + 
  #              log(MCMCpack::dwish(W = Sigma1_prop, v = 2, S = diag(1, d)*1000)) +
  #              mvtnorm::dmvnorm(mu2_prop, mean = rep(0, d), sigma = 1000*diag(d), log = TRUE) + 
  #              log(MCMCpack::dwish(W = Sigma2_prop, v = 2, S = diag(1, d)*1000))
  # lPrior_bot = mvtnorm::dmvnorm(mu1_Init, mean = rep(0, d), sigma = 1000*diag(d), log = TRUE) + 
  #              log(MCMCpack::dwish(W = Sigma1_Init, v = 2, S = diag(1, d)*1000)) +
  #              mvtnorm::dmvnorm(mu2_Init, mean = rep(0, d), sigma = 1000*diag(d), log = TRUE) + 
  #              log(MCMCpack::dwish(W = Sigma2_Init, v = 2, S = diag(1, d)*1000))
  
  # Other examples
  lPrior_top = mvtnorm::dmvnorm(mu1_prop, mean = rep(0, d), sigma = 1000*diag(d), log = TRUE) +
               log(MCMCpack::dwish(W = Sigma1_prop, v = 2, S = diag(1, d)*1000)) +
               mvtnorm::dmvnorm(mu2_prop, mean = rep(0, d), sigma = 1000*diag(d), log = TRUE) +
               log(MCMCpack::dwish(W = Sigma2_prop, v = 2, S = diag(1, d)*1000))
  lPrior_bot = mvtnorm::dmvnorm(mu1_Init, mean = rep(0, d), sigma = 1000*diag(d), log = TRUE) +
               log(MCMCpack::dwish(W = Sigma1_Init, v = 2, S = diag(1, d)*1000)) +
               mvtnorm::dmvnorm(mu2_Init, mean = rep(0, d), sigma = 1000*diag(d), log = TRUE) +
               log(MCMCpack::dwish(W = Sigma2_Init, v = 2, S = diag(1, d)*1000))
  # Proposal
  lProp_top = mvtnorm::dmvnorm(mu1_Init , mean = mu1_prop , sigma = D*(Sigma1_prop/N), log = TRUE) +
              log(MCMCpack::dwish(W = Sigma1_Init, v = 100, S = Sigma1_prop/100)) +
              mvtnorm::dmvnorm(mu2_Init , mean = mu2_prop , sigma = D*(Sigma2_prop/N), log = TRUE) +
              log(MCMCpack::dwish(W = Sigma2_Init, v = 100, S = Sigma2_prop/100))
  if(lProp_top == -Inf){lProp_top = -10^10}
  lProp_bot = mvtnorm::dmvnorm(mu1_prop , mean = mu1_Init , sigma = D*(Sigma1_Init/N), log = TRUE) +
              log(MCMCpack::dwish(W = Sigma1_prop, v = 100, S = Sigma1_Init/100)) + 
              mvtnorm::dmvnorm(mu2_prop , mean = mu2_Init , sigma = D*(Sigma2_Init/N), log = TRUE) +
              log(MCMCpack::dwish(W = Sigma2_prop, v = 100, S = Sigma2_Init/100))
  if(lProp_bot == -Inf){lProp_bot = -10^10}
  
  top = ll_top + lPrior_top + lProp_top; if(is.na(top)){top = -10^10}
  bot = ll_bot + lPrior_bot + lProp_bot
  p_accept = min(1, exp(top-bot)); if(is.nan(p_accept)){p_accept = 0}
  if(runif(1) < p_accept){
    mu1_Init = mu1_prop; Sigma1_Init = Sigma1_prop
    mu2_Init = mu2_prop; Sigma2_Init = Sigma2_prop
    wt_Init = wt_prop
    
    loglik_mat[iter, ] = llTemp_top$P
    Accepted = Accepted + 1
  }
  MH_storage[iter] = sum(log(loglik_mat[iter, ] ))
  print(iter)
}
# Marginal likelihood calculation
post_prob = exp(MH_storage - max(MH_storage))/sum(exp(MH_storage - max(MH_storage)))
M2_Marglik = MH_storage[which(post_prob == max(post_prob))] - log(max(post_prob))

# ELPD Calculation      
MS2 = loo::loo(log(loglik_mat[(nIter/2):(nIter-1), ]))
#Result[m, ] = c(beta_MLE,colMeans(Beta_Storage[((burn_in +1):nIter),]),MS1$elpd_loo,MS1$se_elpd_loo,(Accepted/nIter)*100)

#  output     
c( unique(M2_Marglik), 
  MS2$elpd_loo,
  MS2$se_elpd_loo#, 
  #(Accepted/nIter)*100
  )
# }
# 
# result
# colMeans(result)

######## Output
m1_out = c(unique(M1_Marglik), MS1$elpd_loo, MS1$se_elpd_loo)
m2_out = c(unique(M2_Marglik), MS2$elpd_loo, MS2$se_elpd_loo)
if(length(m1_out)==3){
  M1_Results[m, ] = c(unique(M1_Marglik), MS1$elpd_loo, MS1$se_elpd_loo)
}
if(length(m2_out)==3){
  M2_Results[m, ] = c(unique(M2_Marglik), MS2$elpd_loo, MS2$se_elpd_loo)
}

}

M1_varying_epsilon[which(epsilon == epsilon_grid), ] = c(epsilon, colMeans(M1_Results, na.rm = TRUE))
M2_varying_epsilon[which(epsilon == epsilon_grid), ] = c(epsilon, colMeans(M2_Results, na.rm = TRUE))
}


# Results -----------------------------------------------------------------
M1_varying_epsilon
M2_varying_epsilon


