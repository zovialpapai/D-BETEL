# Packages ----------------------------------------------------------------
library(DEoptim)
library(ggplot2)
library(ggpubr)
library(corrplot)
library(XICOR)

library(rpart)
library(rpart.plot)
library(caret)

library(tibble)
library(cvms)


# Loading Data ------------------------------------------------------------
Data = read.csv("C:/Users/zovia/Downloads/Data_speech.csv")
dim(Data)
colnames(Data)

# The severity of depression ranged from 0 to 27 with a score
# from 0-4 considered none or minimal, 5-9 mild, 10-14 moderate,
# 15-19 moderately severe, and 20-27 severe. 


# Exploration -------------------------------------------------------------
Data$BiologicalSex = Data$BiologicalSex + 1 
#### Moderately imbalanced Data
sum(Data$BiologicalSex==1)/length(Data$BiologicalSex) 

#### Biological sex specific variations
quantile(Data[which(Data$BiologicalSex==1) ,]$PHQ8_Score)
quantile(Data[which(Data$BiologicalSex==2) ,]$PHQ8_Score)

# Histogram by group in ggplot2 (# No normalization)
Data$BiologicalSex = as.factor(Data$BiologicalSex)
ggplot(Data, aes(x = PHQ8_Score, fill = BiologicalSex)) + 
  geom_histogram(colour = "black",
                 lwd = 0.75, alpha = 0.7,
                 linetype = 1, bins = 10,
                 position = "identity") +
  labs(x = "PHQ8 Score") +
  ggtitle("Comparison across biological sex")

#### Correlation plot
Data$BiologicalSex = as.numeric(Data$BiologicalSex)
colnames(Data)[2:20] <- c(paste0("x", 1:17),"y", "A")
colnames(Data)

# Preliminaries -----------------------------------------------------------
#### Sub-populations
set.seed(3076)
Part_1 = Data[Data$A == 1, -1]
Part_2 = Data[Data$A == 2, -1 ]
colnames(Part_1)
# Group fairness by Post Processing ---------------------------------------------------------
#################################### Step 1: Fit LM (or the predictive model)
mod_Part_1 = lm(y ~ ., data = Part_1[, -19])
summary(mod_Part_1)

mod_Part_2 = lm(y ~ ., data =  Part_2[, -19])
summary(mod_Part_2)

################################## Step 2: Group fairness by reweighting
### Plotting outputs before reweighting
y_A =mod_Part_1$fitted.values
y_B = mod_Part_2$fitted.values
plot(sort(y_A), cumsum(rep(1/length(y_A), length(y_A))), 
     type = "l",
     xlab = "h(x)", ylab = "Probabilty",
     main = "cdf of h(x)", col = "red")
lines(sort(y_B), cumsum(rep(1/length(y_B), length(y_B))), 
      type = "l", col = "black")

### Objective function to carry out re-weighting
objective_fn <- function(w){
  #### Calculation of squared Wasserstein-2
  sorted_y_A = sort(y_A); sorted_y_B = sort(y_B)
  
  N1 = length(sorted_y_A); N2 = length(sorted_y_B)
  
  cumsum_y_A = cumsum(rep(1/N1, N1)); 
  cumsum_y_B = cumsum(w) # Change it to be weighted
  
  Quantile_y_A <- function(p){
    #print(sum(cumsum_y_A < p))
    return(sorted_y_A[sum(cumsum_y_A < p)])
  }
  #Quantile_y_A(0.2)
  Quantile_y_B <- function(p){
    #print(sum(cumsum_y_B < p))
    return(sorted_y_B[sum(cumsum_y_B < p)])
  }
  #Quantile_y_B(0.2)
  integrand <- function(p){
    return((Quantile_y_A(p) - Quantile_y_B(p))^2)
  }
  #integrand(0.5)
  
  # To avoid -Inf, Inf
  int_lb = max(min(cumsum_y_A), min(cumsum_y_B))
  int_ub = min(max(cumsum_y_A), max(cumsum_y_B))
  answer = integrate(Vectorize(integrand), int_lb, int_ub, subdivisions = 2000)$value
  answer = 2*answer # for univariate discrete case
  #### Objective function (to be minimised)
  return( (1-lambda)*((sum(w*log(w)))/log(N2)) + lambda*answer + 1e2*(round(sum(w), 3) - 1)^2)
  
}

### Testing the function 
## Case 1
lambda = 1
N2 = length(y_B)
w_example = rep(1/N2, N2)
objective_fn(w = w_example)


### Optimization for reweighting
maxIt <- 1000                     
Out= DEoptim(fn = objective_fn, 
             lower = rep(10^(-7), N2), 
             upper = rep(1 - 10^(-7), N2),
             control=list(NP = 2*N2, itermax = maxIt)
)
fair_w = Out$optim$bestmem 
sum(Out$optim$bestmem) # Check
#### Summary
objective_fn(rep(1/N2, N2))
objective_fn(fair_w)

# Group fairness in model fitting -----------------------------------------------------
### Starting value ("good" starting values based on previous step)
# Weights
starting_w = fair_w ####
# Regression parameters
colnames(Part_1)
mod_Part_1 = lm(y ~ ., data = Part_1[, -19])
starting_mod_Part_1  = coefficients(mod_Part_1)

mod_Part_2 = lm(y ~ ., data = Part_2[, -19])
starting_mod_Part_2 = coefficients(mod_Part_2)

### Objective funstion for fair model fitting
objective_fn <- function(param){
  # Parameters
  beta_1 = as.vector(param[1:18]); beta_2 = as.vector(param[19:36]); 
  w = as.vector(param[-c(1:36)])
  
  #### Calculation of Loss function
  y_1 = beta_1%*%t(as.matrix(cbind(rep(1, nrow(Part_1)),Part_1[,1:17])))
  y_2 = beta_2%*%t(as.matrix(cbind(rep(1, nrow(Part_2)),Part_2[,1:17])))
  
  l_1 = sum((Part_1[, 18] - y_1)^2)/length(y_1)
  l_2 = sum(w*((Part_2[, 18] - y_2)^2))
  
  #### Calculation of squared Wasserstein-2
  sorted_y_1 = sort(y_1); sorted_y_2 = sort(y_2)
  
  N1 = length(sorted_y_1); N2 = length(sorted_y_2)
  
  cumsum_y_1 = cumsum(rep(1/N1, N1)); 
  cumsum_y_2 = cumsum(w) # Change it to be weighted
  
  Quantile_y_1 <- function(p){
    #print(sum(cumsum_y_A < p))
    return(sorted_y_1[sum(cumsum_y_1 < p)])
  }
  #Quantile_y_A(0.2)
  Quantile_y_2 <- function(p){
    #print(sum(cumsum_y_B < p))
    return(sorted_y_2[sum(cumsum_y_2 < p)])
  }
  #Quantile_y_B(0.2)
  integrand <- function(p){
    return((Quantile_y_1(p) - Quantile_y_2(p))^2)
  }
  #integrand(0.5)
  
  # To avoid -Inf, Inf
  int_lb = max(min(cumsum_y_1), min(cumsum_y_2))
  int_ub = min(max(cumsum_y_1), max(cumsum_y_2))
  answer = integrate(Vectorize(integrand), int_lb, int_ub, subdivisions = 2000)$value
  answer = 2*answer # discrete univariate case
  #### Objective function (to be minimised)
  return( l_1 + l_2 + (1-lambda)*((sum(w*log(w)))/log(N2)) + lambda*answer + 1e2*(round(sum(w), 3) - 1)^2)
  
}

### Testing the functions
## Case 1
lambda = 1
N2 = length(y_B)
#w_example = rep(1/N2, N2)
#objective_fn(c(beta_A, beta_B, w_example))

## Case 2
#lambda = 1
#set.seed(1234)
#unnorm_w_example = runif(N2); w_example = unnorm_w_example/(sum(unnorm_w_example))
#objective_fn(c(beta_A, beta_B, w_example))

### Optimization to carry out fair model fitting
maxIt <- 1500    
param_Init = c(starting_mod_Part_1,
               starting_mod_Part_2,
               fair_w)
Init_box = cbind(0.9*param_Init, 1.1*param_Init)

Out= DEoptim(fn = objective_fn, 
             lower = apply(Init_box, 1, min), 
             upper = apply(Init_box, 1, max),
             control=list(NP = 2*(7 + 7 + N2), itermax = maxIt)
)
res = Out$optim$bestmem 

### Summary
Coeff_Part_1 = res[1:18]; Coeff_Part_2 = res[19:36]; fair_w_modfit = res[-c(1:36)]
sum(fair_w_modfit) # Check
objective_fn(res)

Out_fn <- function(param){
  # Parameters
  beta_1 = as.vector(param[1:18]); beta_2 = as.vector(param[19:36]); 
  w = as.vector(param[-c(1:36)])
  
  #### Calculation of Loss function
  y_1 = beta_1%*%t(as.matrix(cbind(rep(1, nrow(Part_1)),Part_1[,1:17])))
  y_2 = beta_2%*%t(as.matrix(cbind(rep(1, nrow(Part_2)),Part_2[,1:17])))
  
  l_1 = sum((Part_1[, 18] - y_1)^2)/length(y_1)
  l_2 = sum(w*((Part_2[, 18] - y_2)^2))
  
  #### Calculation of squared Wasserstein-2
  sorted_y_1 = sort(y_1); sorted_y_2 = sort(y_2)
  
  N1 = length(sorted_y_1); N2 = length(sorted_y_2)
  
  cumsum_y_1 = cumsum(rep(1/N1, N1)); 
  cumsum_y_2 = cumsum(w) # Change it to be weighted
  
  Quantile_y_1 <- function(p){
    #print(sum(cumsum_y_A < p))
    return(sorted_y_1[sum(cumsum_y_1 < p)])
  }
  #Quantile_y_A(0.2)
  Quantile_y_2 <- function(p){
    #print(sum(cumsum_y_B < p))
    return(sorted_y_2[sum(cumsum_y_2 < p)])
  }
  #Quantile_y_B(0.2)
  integrand <- function(p){
    return((Quantile_y_1(p) - Quantile_y_2(p))^2)
  }
  #integrand(0.5)
  
  # To avoid -Inf, Inf
  int_lb = max(min(cumsum_y_1), min(cumsum_y_2))
  int_ub = min(max(cumsum_y_1), max(cumsum_y_2))
  answer = integrate(Vectorize(integrand), int_lb, int_ub, subdivisions = 2000)$value
  answer = 2*answer # Discrete univariate case
  
  #### Parts of Objective function 
  return(list("y_1" = y_1, "y_2" = y_2, "w" = w,
              "l_1" = l_1, "l_2" = l_2, 
              "norm_entropy" = ((sum(w*log(w)))/log(N2)),
              "Wasserstein2" = answer))
}
Out_fn(res)

# Plotting ----------------------------------------------------------------
########################################## Vanilla (no fairness constraint)

Vanilla_data1 = data.frame(sort(mod_Part_1$fitted.values), cumsum(rep(1/length(y_A), length(y_A))))
colnames(Vanilla_data1) <- c("sort_y_1", "cumprob_1")

Vanilla_data2 = data.frame(sort(mod_Part_2$fitted.values), cumsum(rep(1/length(y_B), length(y_B))))
colnames(Vanilla_data2) <- c("sort_y_2", "cumprob_2")
sp1 <- ggplot()+
  geom_line(data = Vanilla_data1, aes(x = sort_y_1, y = cumprob_1), size = 1.5, color = "aquamarine3")+
  geom_line(data = Vanilla_data2, aes(x = sort_y_2, y = cumprob_2), size = 1.5, color = "coral1") +
  geom_point() + labs(x = "h") + labs(y = "cdf of h")+
  ggtitle("Unconstrained") +
  theme(
    # LABLES APPEARANCE
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #panel.background = element_rect(fill = "transparent",colour = NA),
    #plot.background = element_rect(fill = "transparent",colour = NA),
    plot.title = element_text(hjust = 0.5, size=12, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=12, face="bold", colour = "black"),
    axis.title.y = element_text(size=12, face="bold", colour = "black"),
    axis.text.x = element_text(size=12, face="bold", colour = "black"),
    axis.text.y = element_text(size=12, face="bold", colour = "black"),
    strip.text.x = element_text(size = 14, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 14, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3)
  )


########################################## Fair Postprocess
Fair_data1 = data.frame(sort(mod_Part_1$fitted.values), cumsum(rep(1/length(y_A), length(y_A))))
colnames(Fair_data1) <- c("sort_y_1", "cumprob_1")

Fair_data2 = data.frame(sort(mod_Part_2$fitted.values), cumsum(fair_w))
colnames(Fair_data2) <- c("sort_y_2", "cumprob_2")
sp2 <- ggplot()+
  geom_line(data = Fair_data1, aes(x = sort_y_1, y = cumprob_1), size = 1.5, color = "aquamarine3")+
  geom_line(data = Fair_data2, aes(x = sort_y_2, y = cumprob_2), size = 1.5, color = "coral1") +
  geom_point() + labs(x = "h") + labs(y = "cdf of h")+
  ggtitle("Two-step") +
  theme(
    # LABLES APPEARANCE
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #panel.background = element_rect(fill = "transparent",colour = NA),
    #plot.background = element_rect(fill = "transparent",colour = NA),
    plot.title = element_text(hjust = 0.5, size=12, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=12, face="bold", colour = "black"),
    axis.title.y = element_text(size=12, face="bold", colour = "black"),
    axis.text.x = element_text(size=12, face="bold", colour = "black"),
    axis.text.y = element_text(size=12, face="bold", colour = "black"),
    strip.text.x = element_text(size = 14, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 14, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3)
  )

############################################# Fair Model fitting
Fair_data1 = data.frame(sort(mod_Part_1$fitted.values), cumsum(rep(1/length(y_A), length(y_A))))
colnames(Fair_data1) <- c("sort_y_1", "cumprob_1")

Fair_data2 = data.frame(sort(mod_Part_2$fitted.values), cumsum(fair_w_modfit))
colnames(Fair_data2) <- c("sort_y_2", "cumprob_2")

sp3 <- ggplot()+
  geom_line(data = Fair_data1, aes(x = sort_y_1, y = cumprob_1, color = "gender 1"), size = 1.5)+
  geom_line(data = Fair_data2, aes(x = sort_y_2, y = cumprob_2, color = "gender 2"), size = 1.5) +
  geom_point() + labs(x = "h") + labs(y = "cdf of h")+
  ggtitle("In model") +
  scale_colour_manual("",
                      breaks = c("gender 1", "gender 2"),
                      values = c("gender 2"= "coral1","gender 1"="aquamarine3")) +
  theme(
    # LABLES APPEARANCE
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #panel.background = element_rect(fill = "transparent",colour = NA),
    #plot.background = element_rect(fill = "transparent",colour = NA),
    plot.title = element_text(hjust = 0.5, size=12, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=12, face="bold", colour = "black"),
    axis.title.y = element_text(size=12, face="bold", colour = "black"),
    axis.text.x = element_text(size=12, face="bold", colour = "black"),
    axis.text.y = element_text(size=12, face="bold", colour = "black"),
    strip.text.x = element_text(size = 14, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 14, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3)
  )




ggarrange(sp1, sp2, sp3, ncol = 3, nrow = 1, widths = c(1, 1, 1.4))


# Regression plots --------------------------------------------------------
par(mfrow=c(1, 2))
plot(1:length(starting_mod_Part_1[-1]),starting_mod_Part_1[-1],
     ylim = c(-50, 50), pch = 15,
     xlab = "variable index", ylab = "Maximum likelihood estimate", main = "Two-step vs In model \n (Biological gender 1)")
points(1:length(Coeff_Part_1[-1]),Coeff_Part_1[-1],
      pch = 16, col = "red")
legend("bottomright", legend=c("Two-step", "In model"),
       cex=0.85, pch = 15:16, col = c("black", "red"))

plot(1:length(starting_mod_Part_2[-1]),starting_mod_Part_2[-1],
     ylim = c(-50, 50), pch = 15,
     xlab = "variable index", ylab = "Maximum likelihood estimate", main = "Two-step vs In model \n (Biological gender 2)")
points(1:length(Coeff_Part_2[-1]),Coeff_Part_2[-1],
       pch = 16, col = "red")
legend("bottomright", legend=c("Two-step", "In model"),
       cex=0.85, pch = 15:16, col = c("black", "red"))

