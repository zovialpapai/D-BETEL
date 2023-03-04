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

# controls ---------------------------------------------------------
# Skewsness of skew normal distribution
d = 2 # dimension
Alpha = rep(3, d) 
# temparature of likelihood
eta_grid = c(0.25, 0.5, 0.75, 1.0)
# Storage
final_output = array(NA, c(length(eta_grid)+1, 4, 3))

# Run ---------------------------------------------------------------------
for(eta in eta_grid){
# Data generation control
# alpha = 3, 3.5
d = 2 # dimension

#### Simulation control
nrep = 5
nIter = 1000

#### Sample size
# N = 500
N_grid = c(100, 200, 300, 500)
BF = matrix(NA, nrow = length(N_grid), ncol = 3)
colnames(BF) <- c("Sample size", "M1", "M2")
# Repeated simulation -----------------------------------------------------
for(N in N_grid){
# Storage
result = matrix(NA, ncol = 2, nrow = nrep)

for(m in 1:nrep){

# Data generation ---------------------------------------------------------

set.seed(123+(m-1))
x = rmsn(n = N, xi = rep(0, d), Omega = diag(1, d), alpha = Alpha,  tau = 0, dp = NULL)

# M1: Simulations -------------------------------------------------------------
# MLE
mu_MLE = colMeans(x); var_MLE = var(x)

# Sampler
# Storage
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

  # M-H Ratio
  # Parametric Likelihood
  ll_top = eta*sum(mvtnorm::dmvnorm(x, mean = mu_prop, sigma = Sigma_prop, log = TRUE))
  ll_bot = eta*sum(mvtnorm::dmvnorm(x, mean = mu_Init, sigma = Sigma_Init, log = TRUE))
  # Prior
  lPrior_top = mvtnorm::dmvnorm(mu_prop, mean = rep(0, d), sigma = 1000*diag(d), log = TRUE) +
    log(MCMCpack::dwish(W = Sigma_prop, v = 2, S = diag(1, d)*10))
  lPrior_bot = mvtnorm::dmvnorm(mu_Init, mean = rep(0, d), sigma = 1000*diag(d), log = TRUE) +
    log(MCMCpack::dwish(W = Sigma_Init, v = 2, S = diag(1, d)*10))
  # Proposal
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

    MH_storage[iter] = top
    Accepted = Accepted + 1
  }else{
    MH_storage[iter] = bot
  }
  mu_Storage[iter,] = mu_Init
  Sigma_Storage[ , , iter] = Sigma_Init

  print(iter)
}
# Marginal likelihood calculation
post_prob = exp(MH_storage - max(MH_storage))/sum(exp(MH_storage - max(MH_storage)))
M1_Marglik = MH_storage[which(post_prob == max(post_prob))] - log(max(post_prob))


# #  output
# c(unique(M1_Marglik), (Accepted/nIter)*100
# )

print(paste0("M1: repeatation ",m, " of ", nrep, " completed!"))




# M2: Simulations ---------------------------------------------------------

# Initialization
mu1_Init = mu2_Init = colMeans(x); Sigma1_Init = Sigma2_Init = var(x); wt_Init = 0.5
# Record Keeping
MH_storage = matrix(0, nrow = nIter, ncol = 1) # Dunson Miller fills it up with priors
Accepted = 0; D = 1

# Draw Samples
for(iter in 1:nIter){
  # Proposal
  mu1_prop = mvrnorm(n = 1, mu = mu1_Init, Sigma = D*(Sigma1_Init/N))
  mu2_prop = mvrnorm(n = 1, mu = mu2_Init, Sigma = D*(Sigma2_Init/N))
  Sigma1_prop = MCMCpack::rwish(v = 100, S = Sigma1_Init/100)
  Sigma2_prop = MCMCpack::rwish(v = 100, S = Sigma2_Init/100)
  wt_prop = runif(1)
  # M-H Ratio
  # # Parametric Likelihood
  ll_top = eta*sum(log(wt_prop*mvtnorm::dmvnorm(x, mean = mu1_prop, sigma = Sigma1_prop, log = FALSE) + (1- wt_prop)*mvtnorm::dmvnorm(x, mean = mu2_prop, sigma = Sigma2_prop, log = FALSE)))
  ll_bot = eta*sum(log(wt_Init*mvtnorm::dmvnorm(x, mean = mu1_Init, sigma = Sigma1_Init, log = FALSE) + (1- wt_Init)*mvtnorm::dmvnorm(x, mean = mu2_prop, sigma = Sigma2_Init, log = FALSE)))
  # Prior
  lPrior_top = mvtnorm::dmvnorm(mu1_prop, mean = rep(0, d), sigma = 1000*diag(d), log = TRUE) +
    log(MCMCpack::dwish(W = Sigma1_prop, v = 2, S = diag(1, d)*10)) +
    mvtnorm::dmvnorm(mu2_prop, mean = rep(0, d), sigma = 1000*diag(d), log = TRUE) +
    log(MCMCpack::dwish(W = Sigma2_prop, v = 2, S = diag(1, d)*10))
  lPrior_bot = mvtnorm::dmvnorm(mu1_Init, mean = rep(0, d), sigma = 1000*diag(d), log = TRUE) +
    log(MCMCpack::dwish(W = Sigma1_Init, v = 2, S = diag(1, d)*10)) +
    mvtnorm::dmvnorm(mu2_Init, mean = rep(0, d), sigma = 1000*diag(d), log = TRUE) +
    log(MCMCpack::dwish(W = Sigma2_Init, v = 2, S = diag(1, d)*10))
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

    MH_storage[iter] = top
    Accepted = Accepted + 1
  }else{
    MH_storage[iter] = bot
  }
}

# Marginal likelihood calculation
post_prob = exp(MH_storage - max(MH_storage))/sum(exp(MH_storage - max(MH_storage)))
M2_Marglik = MH_storage[which(post_prob == max(post_prob))] - log(max(post_prob))
print(paste0("M2: repeatation ",m, " of ", nrep, " completed!"))

# #  output
# c(unique(M2_Marglik), (Accepted/nIter)*100)


# COMPARISON --------------------------------------------------------------
result[m, ] = c(M1_Marglik, M2_Marglik)
}
out = colMeans(result)
prob = exp(out - max(out))/sum(exp(out - max(out)))

BF[which(N==N_grid), ] = c(N, prob)
}
BF
final_output[which(eta==eta_grid), , ] = BF
}
if(Alpha[1]==2.5){
  final_output[length(eta_grid) + 1, , ] = output_alpha_2_5 
  final_output_2_5 = final_output
}
if(Alpha[1]==3){
  final_output[length(eta_grid) + 1, , ] = output_alpha_3  
  final_output_3 = final_output
}
if(Alpha[1]==3.5){
  final_output[length(eta_grid) + 1, , ] = output_alpha_3_5 
  final_output_3_5 = final_output
}


# D-BETEL summaries -------------------------------------------------------
# Calculate BF 
calculate_BF <- function(x){
  x[-1] = exp(x[-1] - min(x[-1]))/sum(exp(x[-1] - min(x[-1])))
  return(x)
}
# D-BETEL Summaries: alpha = 2.5 
sample_size = c(100, 200, 300, 500)
M1 = c(-458.20, -1055.72, -1708.96, -3106.65)
M2 = c(-463.52, -1068.11, -1722.64, -3126.49)
output_alpha_2_5 = cbind(sample_size, M1, M2)
colnames(output_alpha_2_5) <- c("sample_size", "M1", "M2")

for(k in 1:nrow(output_alpha_2_5)){
  output_alpha_2_5[k, ] = calculate_BF(output_alpha_2_5[k, ])
}

# D-BETEL Summaries: alpha = 3 
sample_size = c(100, 200, 300, 500)
M1 = c(-457.89, -1055.52, -1708.05, -3106.65)
M2 = c(-463.92, -1066.61, -1720.44, -3126.49)
output_alpha_3 = cbind(sample_size, M1, M2)
colnames(output_alpha_3) <- c("sample_size", "M1", "M2")

for(k in 1:nrow(output_alpha_3)){
  output_alpha_3[k, ] = calculate_BF(output_alpha_3[k, ])
}

# D-BETEL Summaries: alpha = 3.5 
sample_size = c(100, 200, 300, 500)
M1 = c(-457.67, -1055.38, -1708.12, -3106.65)
M2 = c(-463.58, -1066.55, -1719.56, -3126.49)
output_alpha_3_5 = cbind(sample_size, M1, M2)
colnames(output_alpha_3_5) <- c("sample_size", "M1", "M2")

for(k in 1:nrow(output_alpha_3_5)){
  output_alpha_3_5[k, ] = calculate_BF(output_alpha_3_5[k, ])
}




# Plots: D-BETEL vs Fractional posterior-------------------------------------------------------------------

par(mfrow=c(1, 3))

plot(final_output_2_5[5,  ,1], final_output[5,  ,2], type = "l", lwd = 4, lty = 1,
     ylim = c(-0.2, 1.05), col = 1, 
     xlab = "Sample size", ylab = "P(Model 1 | data)")
for(k in 1:(length(eta_grid))){
lines(final_output_2_5[k,  ,1], final_output_2_5[k,  ,2], type = "l", lwd = 4, lty = 2,
      ylim = c(0, 1.05), col = k+1)  
}
legend("bottomleft", legend=c("D-BETEL", "Fractional posterior 0.25", "Fractional posterior 0.50", "Fractional posterior 0.75" ,"Standard posterior"),
       col=1:6, lty=c(1, rep(2, 4)), cex=0.8, lwd = 4, text.font = 2)

plot(final_output_3[5,  ,1], final_output_3[5,  ,2], type = "l", lwd = 4, lty = 1,
     ylim = c(0, 1.05), col = 1,
     xlab = "Sample size", ylab = "P(Model 1 | data)")
for(k in 1:(length(eta_grid))){
  lines(final_output_3[k,  ,1], final_output_3[k,  ,2], type = "l", lwd = 4, lty = 2,
        ylim = c(0, 1.05), col = k+1)  
}
# legend("bottomleft", legend=c("D-BETEL", "FracPost-0.25", "FracPost-0.50", "FracPost-0.75" ,"StandardPost"),
#        col=1:6, lty=c(1, rep(2, 4)), cex=0.8, lwd = 4)

plot(final_output_3_5[5,  ,1], final_output[5,  ,2], type = "l", lwd = 4, lty = 1,
     ylim = c(0, 1.05),  col = 1,
     xlab = "Sample size", ylab = "P(Model 1 | data)")
for(k in 1:(length(eta_grid))){
  lines(final_output_3_5[k,  ,1], final_output_3_5[k,  ,2], type = "l", lwd = 4, lty = 2,
        ylim = c(0, 1.05), col = k+1)  
}
# legend("bottomleft", legend=c("D-BETEL", "FracPost-0.25", "FracPost-0.50", "FracPost-0.75" ,"StandardPost"),
#        col=1:6, lty=c(1, rep(2, 4)), cex=0.8, lwd = 4)

  
  
  


