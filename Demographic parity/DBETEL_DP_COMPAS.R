# Dependencies ------------------------------------------------------------
library(DEoptim)

library(ggplot2)
library(egg)

library(dplyr)
library(hrbrthemes)

library(reshape2)

# Raw data (obtained from Kaggle) -------------------------------------------------------
# The raw dataset is attached in the data folder of the supplement
Raw_Data = read.csv("C:/Users/zovia/OneDrive/Desktop/compas-scores-raw.csv/compas-scores-raw.csv")
dim(Raw_Data)
colnames(Raw_Data)

# Subsetting, Pipe
Raw_Data = Raw_Data[Raw_Data$Ethnic_Code_Text %in% c("African-American", 
                                                     "African-Am",
                                                     "Caucasian"), ]
Raw_Data = Raw_Data[Raw_Data$DisplayText %in% c("Risk of Recidivism"),]
Raw_Data = Raw_Data[Raw_Data$AssessmentType %in% c("New"),]
# View(Raw_Data)
dim(Raw_Data)

# Data Preprocessing -----------------------------------------------------------

#### Response
y = Raw_Data$RawScore
hist(y, breaks = 20)

#### Protected variable
Raw_Data$Ethnic_Code_Text = 1*(Raw_Data$Ethnic_Code_Text %in% c("African-American", "African-Am"))    

#### Explanatory variables preprocessing
DOB = as.numeric(unlist(substr(Raw_Data$DateOfBirth, 7, 8)))
# max(DOB); min(DOB); 
x = (100 - DOB) + 13 # Age at 2013, Assumption?
hist(x)
Raw_Data$DateOfBirth = x

Raw_Data$Sex_Code_Text = as.numeric(Raw_Data$Sex_Code_Text=="Male")
Raw_Data$MaritalStatus = as.numeric(Raw_Data$MaritalStatus=="Married")
Raw_Data$Language = as.numeric(Raw_Data$Language=="English")
Raw_Data$LegalStatus = as.numeric(Raw_Data$LegalStatus=="Pretrial")
Raw_Data$CustodyStatus = as.numeric(Raw_Data$CustodyStatus %in% c("Jail Inmate", "Prison Inmate"))
# View(Raw_Data)

#### Our full dataset
Data = Raw_Data[ , c("RawScore",
                     "Sex_Code_Text", "Ethnic_Code_Text",
                     "DateOfBirth", "Language",
                     "LegalStatus", "CustodyStatus",
                     "MaritalStatus","RecSupervisionLevel" )]
colnames(Data) <- c("RawScore",
                    "Male_Yes_No", "Black_Yes_No",
                    "Age2013", "English_Yes_No",
                    "Pretrial_Yes_NO", "Inmate_Yes_No",
                    "Married_Yes_No","RecSupervisionLevel" )
View(Data)
dim(Data)


# Preliminaries -----------------------------------------------------------
#### Sub-populations
set.seed(3076)
subset_selected = sample(min(sum(Data$Black_Yes_No == 1), sum(Data$Black_Yes_No == 0)), 
       50, replace = FALSE, )
Part_Black = Data[Data$Black_Yes_No == 1, ]; Part_Black = Part_Black[subset_selected, ]
Part_NonBlack = Data[Data$Black_Yes_No == 0, ]; Part_NonBlack = Part_NonBlack[subset_selected, ]

# hist(Part_Black$RawScore, col=rgb(0,0,1,1/4), freq = F, 
#      ylim = c(0,1), xlim = c(min(Part_Black$RawScore)-1, max(Part_Black$RawScore)+1),
#      breaks = 20,
#      xlab = "Raw score",
#      main = "Recidivism score",
#      axes = FALSE)
# axis(1, pos = 0)
# axis(2, pos =-3)
# hist(Part_NonBlack$RawScore, breaks = 20, add = T, col=rgb(1,0,0,1/4), freq = F)
# legend("topright", c("African-americans", "Non African-americans") ,
#        col = c(rgb(0,0,1,1/4), rgb(1,0,0,1/4)), lty = 1, lwd = 5, cex = 1)

combined_data = data.frame(Score = c(Part_Black$RawScore, Part_NonBlack$RawScore),
                           Race = as.factor(c(rep("African American", length(Part_Black$RawScore)),
                                            rep("Non African American", length(Part_Black$RawScore))) )
                           )
colnames(combined_data) <- c("Raw_score", "Race")
# Histogram by group in ggplot2
ggplot(combined_data, aes(x = Raw_score, fill = Race)) + 
       geom_histogram(colour = "black",
       lwd = 0.75, alpha = 0.5,
       linetype = 1, bins = 10,
       position = "identity") +
       labs(x = "raw score") +
       ggtitle("Comparison of recidivism scores")



#### Linear regression
colnames(Data)
mod = lm(RawScore ~ ., data = Data[, -3])
summary(mod)
# hist(mod$fitted.values[which(Data$Black_Yes_No == 1)], col=rgb(0,0,1,1/4), freq = F,
#      xlab = "Raw Score", main = "Fitted: Black v/s Non-black", ylim = c(0,1), xlim = c(-3,3))
# hist(mod$fitted.values[which(Data$Black_Yes_No == 0)], add = T, col=rgb(1,0,0,1/4), freq = F)
# legend("topright", c("Afro-americans", "Non Afro-americans") ,
#        col = c(rgb(0,0,1,1/4), rgb(1,0,0,1/4)), lty = 1, lwd = 5, cex = 0.7)

# Group fairness by Post Processing ---------------------------------------------------------
#################################### Step 1: Fit LM (or the predictive model)
mod_Black = lm(RawScore ~ ., data = Part_Black[, -3])
summary(mod_Black)

mod_NonBlack = lm(RawScore ~ ., data = Part_NonBlack[, -3])
summary(mod_NonBlack)

# hist(mod_Black$fitted.values, col=rgb(0,0,1,1/4), freq = F,
#      xlab = "Raw Score", main = "Fitted: Afro-americans v/s Non Afro-americans", ylim = c(0,1), xlim = c(-3,3))
# hist(mod_NonBlack$fitted.values, add = T, col=rgb(1,0,0,1/4), freq = F)
# legend("topright", c("Afro-americans", "Non Afro-americans") ,
#        col = c(rgb(0,0,1,1/4), rgb(1,0,0,1/4)), lty = 1, lwd = 5, cex = 0.7)

################################## Step 2: Group fairness by reweighting
### Plotting outputs before reweighting
y_A = mod_NonBlack$fitted.values
y_B = mod_Black$fitted.values
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
objective_fn(w_example)

## Case 2
lambda = 1
set.seed(1234)
unnorm_w_example = runif(N2); w_example = unnorm_w_example/(sum(unnorm_w_example))
objective_fn(w_example)

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

## Output
## Vanilla
# Vanilla_data = data.frame(sort(y_A), cumsum(rep(1/length(y_A), length(y_A))),
#                           sort(y_B), cumsum(rep(1/length(y_B), length(y_B))))
# colnames(Vanilla_data) <- c("sort_y_A", "cumprob_A", "sort_y_B", "cumprob_B")
#saveRDS(Vanilla_data, "E:/temporary_betel/Diagram_Data/GF_Vanilla.rds")
# sp1 <- ggplot()+
#   geom_line(data = Vanilla_data, aes(x = sort_y_A, y = cumprob_A), color = "red")+
#   geom_line(data = Vanilla_data, aes(x = sort_y_B, y = cumprob_B), color = "black") +
#   geom_point() + labs(x = "h(x)") + labs(y = "cdf of h(x)")+
#   ggtitle("vanilla case") +
#   theme(
#     # LABLES APPEARANCE
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.background = element_rect(fill = "transparent",colour = NA),
#     plot.background = element_rect(fill = "transparent",colour = NA),
#     plot.title = element_text(hjust = 0.5, size=12, face= "bold", colour= "black" ),
#     axis.title.x = element_text(size=12, face="bold", colour = "black"),
#     axis.title.y = element_text(size=12, face="bold", colour = "black"),
#     axis.text.x = element_text(size=12, face="bold", colour = "black"),
#     axis.text.y = element_text(size=12, face="bold", colour = "black"),
#     strip.text.x = element_text(size = 14, face="bold", colour = "black" ),
#     strip.text.y = element_text(size = 14, face="bold", colour = "black"),
#     axis.line.x = element_line(color="black", size = 0.3),
#     axis.line.y = element_line(color="black", size = 0.3),
#     panel.border = element_rect(colour = "black", fill=NA, size=0.3)
#   )
# sp1

## Fair
# From saved data
# Fair_data = readRDS("E:/temporary_betel/Diagram_Data/GF_PostProcess.rds")
# fair_w = c(Fair_data[1 ,4], diff(Fair_data[ ,4]))
# From running
Fair_data = data.frame(sort(y_A), cumsum(rep(1/length(y_A), length(y_A))),
                      sort(y_B), cumsum(fair_w))
colnames(Fair_data) <- c("sort_y_A", "cumprob_A", "sort_y_B", "cumprob_B")
#saveRDS(Fair_data, "E:/temporary_betel/Diagram_Data/GF_PostProcess.rds")

# sp2 <- ggplot()+
#   geom_line(data = Fair_data, aes(x = sort_y_A, y = cumprob_A, color = "Non African-americans"))+
#   geom_line(data = Fair_data, aes(x = sort_y_B, y = cumprob_B, color = "African-americans")) +
#   geom_point() + labs(x = "h(x)") + labs(y = "cdf of h(x)")+
#   ggtitle(expression(paste(epsilon, "-fair case"))) +
#   scale_colour_manual("",
#                       breaks = c("Non African-americans", "African-americans"),
#                       values = c("Non African-americans"="red","African-americans"="black")) +
#   theme(
#     # LABLES APPEARANCE
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.background = element_rect(fill = "transparent",colour = NA),
#     plot.background = element_rect(fill = "transparent",colour = NA),
#     plot.title = element_text(hjust = 0.5, size=12, face= "bold", colour= "black" ),
#     axis.title.x = element_text(size=12, face="bold", colour = "black"),
#     axis.title.y = element_text(size=12, face="bold", colour = "black"),
#     axis.text.x = element_text(size=12, face="bold", colour = "black"),
#     axis.text.y = element_text(size=12, face="bold", colour = "black"),
#     strip.text.x = element_text(size = 14, face="bold", colour = "black" ),
#     strip.text.y = element_text(size = 14, face="bold", colour = "black"),
#     axis.line.x = element_line(color="black", size = 0.3),
#     axis.line.y = element_line(color="black", size = 0.3),
#     panel.border = element_rect(colour = "black", fill=NA, size=0.3)
#   )
# sp2
# ggarrange(sp1, sp2, ncol = 2, nrow = 1)


# Group fairness in model fitting -----------------------------------------------------
### Starting value ("good" starting values based on previous step)
# Weights
starting_w = fair_w ####
# Regression parameters
mod_Black = lm(RawScore ~ ., data = Part_Black[, -c(3,5)])
starting_mod_Black = coefficients(mod_Black)

mod_NonBlack = lm(RawScore ~ ., data = Part_NonBlack[, -c(3, 5)])
starting_mod_NonBlack = coefficients(mod_NonBlack)

### Objective funstion for fair model fitting
objective_fn <- function(param){
  # Parameters
  beta_A = as.vector(param[1:7]); beta_B = as.vector(param[8:14]); 
  w = as.vector(param[-c(1:14)])
  
  #### Calculation of Loss function
  y_A = beta_A%*%t(as.matrix(cbind(rep(1, nrow(Part_NonBlack)),Part_NonBlack[,-c(1,3,5)])))
  y_B = beta_B%*%t(as.matrix(cbind(rep(1, nrow(Part_Black)),Part_Black[,-c(1,3,5)])))
  
  l_A = sum((Part_NonBlack[, 1] - y_A)^2)/length(y_A)
  l_B = sum(w*((Part_Black[, 1] - y_B)^2))

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
  answer = 2*answer # discrete univariate case
  #### Objective function (to be minimised)
  return( l_A + l_B + (1-lambda)*((sum(w*log(w)))/log(N2)) + lambda*answer + 1e2*(round(sum(w), 3) - 1)^2)
  
}

### Testing the functions
## Case 1
lambda = 1
N2 = length(y_B)
w_example = rep(1/N2, N2)
#objective_fn(c(beta_A, beta_B, w_example))

## Case 2
lambda = 1
set.seed(1234)
unnorm_w_example = runif(N2); w_example = unnorm_w_example/(sum(unnorm_w_example))
#objective_fn(c(beta_A, beta_B, w_example))

### Optimization to carry out fair model fitting
maxIt <- 1500    
param_Init = c(starting_mod_NonBlack,
               starting_mod_Black,
               fair_w)
Init_box = cbind(0.9*param_Init, 1.1*param_Init)

Out= DEoptim(fn = objective_fn, 
             lower = apply(Init_box, 1, min), 
             upper = apply(Init_box, 1, max),
             control=list(NP = 2*(7 + 7 + N2), itermax = maxIt)
)
res = Out$optim$bestmem 

### Summary
Coeff_NonBlack = res[1:7]; Coeff_Black = res[8:14]; fair_w_modfit = res[-c(1:14)]
sum(fair_w_modfit) # Check
objective_fn(res)

Out_fn <- function(param){
  # Parameters
  beta_A = as.vector(param[1:7]); beta_B = as.vector(param[8:14]); 
  w = as.vector(param[-c(1:14)])
  
  #### Calculation of Loss function
  y_A = beta_A%*%t(as.matrix(cbind(rep(1, nrow(Part_NonBlack)),Part_NonBlack[,-c(1,3,5)])))
  y_B = beta_B%*%t(as.matrix(cbind(rep(1, nrow(Part_Black)),Part_Black[,-c(1,3,5)])))
  
  l_A = sum((Part_NonBlack[, 1] - y_A)^2)/length(y_A)
  l_B = sum(w*((Part_Black[, 1] - y_B)^2))
  
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
  answer = 2*answer # Discrete univariate case
  
  #### Parts of Objective function 
  return(list("y_A" = y_A, "y_B" = y_B, "w" = w,
              "l_A" = l_A, "l_B" = l_B, 
              "norm_entropy" = ((sum(w*log(w)))/log(N2)),
              "Wasserstein2" = answer))
}
Out_fn(res)


# Plotting ----------------------------------------------------------------
########################################## Vanilla (no fairness constraint)
Vanilla_data = data.frame(sort(mod_NonBlack$fitted.values), cumsum(rep(1/length(y_A), length(y_A))),
                          sort(mod_Black$fitted.values), cumsum(rep(1/length(y_B), length(y_B))))
colnames(Vanilla_data) <- c("sort_y_A", "cumprob_A", "sort_y_B", "cumprob_B")
#Vanilla_data = readRDS("E:/temporary_betel/Diagram_Data/Vanilla_data.rds")
sp1 <- ggplot()+
  geom_line(data = Vanilla_data, aes(x = sort_y_A, y = cumprob_A), size = 1.5, color = "aquamarine3")+
  geom_line(data = Vanilla_data, aes(x = sort_y_B, y = cumprob_B), size = 1.5, color = "coral1") +
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
#Fair_data = data.frame(sort(mod_NonBlack$fitted.values), cumsum(rep(1/length(y_A), length(y_A))),
#                       sort(mod_Black$fitted.values), cumsum(fair_w))
#colnames(Fair_data) <- c("sort_y_A", "cumprob_A", "sort_y_B", "cumprob_B")
Fair_data = readRDS("E:/temporary_betel/Diagram_Data/postprocess_data.rds")
sp2 <- ggplot()+
  geom_line(data = Fair_data, aes(x = sort_y_A, y = cumprob_A), size = 1.5, color = "aquamarine3")+
  geom_line(data = Fair_data, aes(x = sort_y_B, y = cumprob_B), size = 1.5, color = "coral1") +
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
Fair_data = data.frame(sort(y_A), cumsum(rep(1/length(y_A), length(y_A))),
                       sort(y_B), cumsum(fair_w_modfit))

colnames(Fair_data) <- c("sort_y_A", "cumprob_A", "sort_y_B", "cumprob_B")
#Fair_data <- readRDS("E:/temporary_betel/Diagram_Data/modelfit_data.rds")
sp3 <- ggplot()+
  geom_line(data = Fair_data, aes(x = sort_y_A, y = cumprob_A, color = "Non African-americans"), size = 1.5)+
  geom_line(data = Fair_data, aes(x = sort_y_B, y = cumprob_B, color = "African-americans"), size = 1.5) +
  geom_point() + labs(x = "h") + labs(y = "cdf of h")+
  ggtitle("In model") +
  scale_colour_manual("",
                      breaks = c("African-americans", "Non African-americans"),
                      values = c("African-americans"= "coral1","Non African-americans"="aquamarine3")) +
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

ggarrange(sp1, sp2, sp3, ncol = 3, nrow = 1)


# # Backup ------------------------------------------------------------------
# Vanilla_data = data.frame(sort(mod_NonBlack$fitted.values), cumsum(rep(1/length(y_A), length(y_A))),
#                           sort(mod_Black$fitted.values), cumsum(rep(1/length(y_B), length(y_B))))
# colnames(Vanilla_data) <- c("sort_y_A", "cumprob_A", "sort_y_B", "cumprob_B")
# #saveRDS(Vanilla_data, "E:/temporary_betel/Diagram_Data/Vanilla_data.rds")
# 
# ## Fair Postprocess
# Fair_data = data.frame(sort(mod_NonBlack$fitted.values), cumsum(rep(1/length(y_A), length(y_A))),
#                        sort(mod_Black$fitted.values), cumsum(fair_w))
# colnames(Fair_data) <- c("sort_y_A", "cumprob_A", "sort_y_B", "cumprob_B")
# #saveRDS(Fair_data, "E:/temporary_betel/Diagram_Data/postprocess_data.rds")
# 
# ## Fair Model fitting
# Fair_data = data.frame(sort(y_A), cumsum(rep(1/length(y_A), length(y_A))),
#                        sort(y_B), cumsum(fair_w_modfit))
# 
# colnames(Fair_data) <- c("sort_y_A", "cumprob_A", "sort_y_B", "cumprob_B")
# #saveRDS(Fair_data, "E:/temporary_betel/Diagram_Data/modelfit_data.rds")


# Regression plots --------------------------------------------------------

par(mfrow=c(1, 2))
plot(1:length(starting_mod_Black[-1]),starting_mod_Black[-1],
     ylim = c(-1, 1), pch = 15,
     xlab = "variable index", ylab = "Maximum likelihood estimate", main = "Two-step vs In model \n (African American)")
points(1:length(Coeff_Black[-1]),Coeff_Black[-1],
       pch = 16, col = "red")
legend("bottomright", legend=c("Two-step", "In model"),
       cex=0.85, pch = 15:16, col = c("black", "red"))

plot(1:length(starting_mod_NonBlack[-1]),starting_mod_NonBlack[-1],
     ylim = c(-1, 1), pch = 15,
     xlab = "variable index", ylab = "Maximum likelihood estimate", main = "Two-step vs In model \n (Non African American)")
points(1:length(Coeff_NonBlack[-1]),Coeff_NonBlack[-1],
       pch = 16, col = "red")
legend("bottomright", legend=c("Two-step", "In model"),
       cex=0.85, pch = 15:16, col = c("black", "red"))

