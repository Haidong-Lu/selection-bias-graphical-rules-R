library(tidyverse)
library(tidymodels)
library(simcausal)

# Set a seed for reproducibility
set.seed(12345)

# Create a simulated dataset
n_samples <- 1000000

data <- tibble(
  a = rbern(n_samples, .5),
  y = rbern(n_samples, .4),
  l = rbern(n_samples, .1 + .3*a + .5*y), #very few control healthy units selected so the effect is overestimated
  s = rbern(n_samples, .1 + .8*l)) #the larger l, the larger s

# Filter data to include only selected observations (S = 1)
selected_data <- data[data$s == 1, ]

# Truth
true_effect <- 0
model <- glm(y ~ a, data = data)
summary(model)

# Crude estimate in selected sample 
model_crude <- glm(y ~ a, data = selected_data)
summary(model_crude)

# Complete-case analysis adjusting for L in selected sample
model_ccl <- glm(y ~ a + l, data = selected_data)
summary(model_ccl)

# IPCW estimator using information on L in unselected sample as well
den.fit <- glm(s ~ l,
               data = data, family = binomial)

num.fit <- glm(s ~ 1,
               data = data, family = binomial)

den.pred <- predict(den.fit, type = "response")

num.pred <- predict(num.fit, type = "response")

selected_data_w <- cbind(data, num.pred, den.pred) %>% 
  filter(s == 1) %>% 
  mutate(w = num.pred/den.pred)

model_ipcw <- glm(y ~ a, weights=w, data=selected_data_w)
summary(model_ipcw)

# Yichi: try a different approach using only biased selected sample data to estimate weight
# issue: not observing S=0
den.fit2 <- glm(s ~ l,
                data = selected_data, family = binomial)

num.fit2 <- glm(s ~ 1,
                data = selected_data, family = binomial)

den.pred2 <- predict(den.fit2, type = "response")

num.pred2 <- predict(num.fit2, type = "response")

selected_data_w2 <- cbind(data, num.pred2, den.pred2) %>% 
  filter(s == 1) %>% 
  mutate(w2 = num.pred2/den.pred2)

model_ipcw2 <- glm(y ~ a, weights=w2, data=selected_data_w2)
summary(model_ipcw2)

# Yichi: try a non-stabilized weighting approach
# still unbiased
den.fit3 <- glm(s ~ l,
                data = data, family = binomial)

den.pred3 <- predict(den.fit3, type = "response")

selected_data_w3 <- cbind(data, den.pred3) %>% 
  filter(s == 1) %>% 
  mutate(w3 = 1/den.pred3)

model_ipcw3 <- glm(y ~ a, weights=w3, data=selected_data_w3)
summary(model_ipcw3)

# Yichi: make sure the scale does not change the estimation
# the same as the first model
den.fit4 <- glm(s ~ l,
                data = data, family = binomial)

num.fit4 <- glm(s ~ 1,
                data = data, family = binomial)

den.pred4 <- predict(den.fit4, type = "response")

num.pred4 <- predict(num.fit4, type = "response")

selected_data_w4 <- cbind(data, num.pred4, den.pred4) %>% 
  filter(s == 1) %>% 
  mutate(w4 = 79*num.pred4/den.pred4)

model_ipcw4 <- glm(y ~ a, weights=w4, data=selected_data_w4)
summary(model_ipcw4)

# Yichi: check the solution by multiplying a missing weight of receiving the observed treatment
# the same as the first model
den.fit5 <- glm(s ~ l,
                data = data, family = binomial)

num.fit5<- glm(s ~ 1,
               data = data, family = binomial)

den.pred5 <- predict(den.fit5, type = "response")

num.pred5 <- predict(num.fit5, type = "response")

selected_data_w5 <- cbind(data, num.pred5, den.pred5) %>% 
  filter(s == 1) %>% 
  mutate(w5 = num.pred5/den.pred5/(-0.2*a+0.5))

model_ipcw5 <- glm(y ~ a, weights=w5, data=selected_data_w5)
summary(model_ipcw5)

# Yichi: check the true IPW estimator
mean((data$a==1)*(data$s==1)*data$y/den.pred) - mean((data$a==0)*(data$s==1)*data$y/den.pred) #without prob of trt.
mean((data$a==1)*(data$s==1)*data$y/den.pred/0.5) - mean((data$a==0)*(data$s==1)*data$y/den.pred/0.5) #with prob of trt. (unbiased) HT
#the second one is unbiased, and did not use any info of Y in the unselected sample because of the indicator S==1
#but we do need the info of the total sample size N.
#new comment: this is different from the regression above using only the selected sample
#and it is the HT estimator, not the Hajek estimator with correspondence to the weighted regression
#so this is unbiased but with different results

#Yichi: check if we can do these in the selected sample (Hajek)
sum(data$y[(data$a==1)&(data$s==1)]/den.pred[(data$a==1)&(data$s==1)]/0.5)/sum(1/den.pred[(data$a==1)&(data$s==1)]/0.5) - 
  sum(data$y[(data$a==0)&(data$s==1)]/den.pred[(data$a==0)&(data$s==1)]/0.5)/sum(1/den.pred[(data$a==0)&(data$s==1)]/0.5) #with prob of trt. (unbiased) Hajek
#new comment: this results in the same conclusions as above, both use only the selected sample
#again this is the Hajek estimator, in theory corresponding to the weighted regression above

#Yichi: Correa estimator (3 approaches the same results - biased)
#APPROACH 1:Filter selected data to include observations w. L = 1 or L = 0
selected_data_l0 <- selected_data[selected_data$l == 0, ]
selected_data_l1 <- selected_data[selected_data$l == 1, ]
#Crude estimate in selected sample stratified by L
#p(l)
model_crude_l0 <- glm(y ~ a, data = selected_data_l0)
model_crude_l1 <- glm(y ~ a, data = selected_data_l1)
summary(model_crude_l0)
summary(model_crude_l1)
#Get the ATE = weighted average by the distribution of p(L) in the general pop.
model_crude_l0$coefficients[2]*(1-mean(data$l))+model_crude_l1$coefficients[2]*mean(data$l)
#APPROACH 2:Crude estimate in selected sample modified by L
model_crude_modified <- glm(y ~ a*l, data = selected_data)
summary(model_crude_modified)
model_crude_modified$coefficients[2]*(1-mean(data$l))+
  (model_crude_modified$coefficients[2]+model_crude_modified$coefficients[4])*mean(data$l)
#APPROACH 3:sample mean estimator
(mean(selected_data_l0$y[selected_data_l0$a == 1])-mean(selected_data_l0$y[selected_data_l0$a == 0]))*
  (1-mean(data$l)) +
  (mean(selected_data_l1$y[selected_data_l1$a == 1])-mean(selected_data_l1$y[selected_data_l1$a == 0]))*
  mean(data$l)
#APPROACH 4 (inspired):try the selected sample distribution of L
(mean(selected_data_l0$y[selected_data_l0$a == 1])-mean(selected_data_l0$y[selected_data_l0$a == 0]))*
  (1-mean(selected_data$l)) +
  (mean(selected_data_l1$y[selected_data_l1$a == 1])-mean(selected_data_l1$y[selected_data_l1$a == 0]))*
  mean(selected_data$l)
#APPROACH 5 (inspired):try the corresponding IPSW estimator by their results





#IPW estimator with variance estimation
#Input: data - contains four columns a, l, s and y (a and y will be used only for s == 1)
#       p_treat - the constant probability for each unit to be treated
#       estimator - the choice of IPW estimator to be estimated
IPW_inference <- function(data, p_treat = 0.5,
                          estimator = c("HT", "Hajek")){
  n_full <- dim(data)[1]
  
  #predicting selection
  s.model <- glm(s ~ l, data = data, family = binomial)
  s.prob.pred <- predict(s.model, type = "response")
  
  if (estimator == "HT"){
    #point estimate
    est.Y1 <- mean(data$a * data$s * data$y / s.prob.pred) / p_treat
    est.Y0 <- mean((1 - data$a) * data$s * data$y / s.prob.pred) / (1 - p_treat)
    est.effect <- est.Y1 - est.Y0
    
    #variance estimate - matrix components in sandwich estimator
    psi <- cbind(data$s - s.prob.pred, 
                 data$l * (data$s - s.prob.pred), 
                 (data$s * data$a * data$y) / (s.prob.pred * p_treat) - est.Y1,
                 (data$s * (1 - data$a) * data$y) / (s.prob.pred * (1 - p_treat)) - est.Y0)
    A_psi <- matrix(c(mean(s.prob.pred * (1 - s.prob.pred)), 
                      mean(data$l * s.prob.pred * (1 - s.prob.pred)), 
                      0, 
                      0,
                      mean(data$l * s.prob.pred * (1 - s.prob.pred)), 
                      mean((data$l)^2 * s.prob.pred * (1 - s.prob.pred)), 
                      0, 
                      0,
                      mean((data$s * data$a * data$y * (1 - s.prob.pred)) / (s.prob.pred * p_treat)), 
                      mean((data$s * data$a * data$y * data$l * (1 - s.prob.pred)) / (s.prob.pred * p_treat)),
                      1,
                      0,
                      mean((data$s * (1 - data$a) * data$y * (1 - s.prob.pred)) / (s.prob.pred * (1 - p_treat))),
                      mean((data$s * (1 - data$a) * data$y * data$l * (1 - s.prob.pred)) / (s.prob.pred * (1 - p_treat))),
                      0,
                      1),
                    4, 4, byrow = TRUE)
    A_psi_inv <- solve(A_psi)
    B_psi <- t(psi) %*% psi / dim(data)[1]
  } else if (estimator == "Hajek"){
    #point estimate
    est.Y1 <- mean(data$a * data$s * data$y / s.prob.pred) / mean(data$a * data$s / s.prob.pred)
    est.Y0 <- mean((1 - data$a) * data$s * data$y / s.prob.pred) / mean((1 - data$a) * data$s / s.prob.pred)
    est.effect <- est.Y1 - est.Y0
    
    #variance estimate - matrix components in sandwich estimator
    psi <- cbind(data$s - s.prob.pred, 
                 data$l * (data$s - s.prob.pred), 
                 (data$s * data$a * (data$y - est.Y1)) / s.prob.pred,
                 (data$s * (1 - data$a) * (data$y - est.Y0)) / s.prob.pred)
    A_psi <- matrix(c(mean(s.prob.pred * (1 - s.prob.pred)), 
                      mean(data$l * s.prob.pred * (1 - s.prob.pred)), 
                      0, 
                      0,
                      mean(data$l * s.prob.pred * (1 - s.prob.pred)), 
                      mean((data$l)^2 * s.prob.pred * (1 - s.prob.pred)), 
                      0, 
                      0,
                      mean((data$s * data$a * (data$y - est.Y1) * (1 - s.prob.pred)) / s.prob.pred), 
                      mean((data$s * data$a * (data$y - est.Y1) * data$l * (1 - s.prob.pred)) / s.prob.pred),
                      mean(data$s * data$a / s.prob.pred),
                      0,
                      mean((data$s * (1 - data$a) * (data$y - est.Y0) * (1 - s.prob.pred)) / s.prob.pred),
                      mean((data$s * (1 - data$a) * (data$y - est.Y0) * data$l * (1 - s.prob.pred)) / s.prob.pred),
                      0,
                      mean(data$s * (1 - data$a) / s.prob.pred)),
                    4, 4, byrow = TRUE)
    A_psi_inv <- solve(A_psi)
    B_psi <- t(psi) %*% psi / dim(data)[1]
  }
  cov_mat <- A_psi_inv %*% B_psi %*% t(A_psi_inv) / n_full
  est_var <- cov_mat[3, 3] + cov_mat[4, 4] - 2 * cov_mat[3, 4]
  return(list(point_est = est.effect,
              var_est = est_var,
              cov_mat = cov_mat))
}

HT <- IPW_inference(data = data, p_treat = 0.5, estimator = "HT")
HT
CI_HT <- c(HT$point_est - 1.96 * sqrt(HT$var_est),
           HT$point_est + 1.96 * sqrt(HT$var_est))
CI_HT
Hajek <- IPW_inference(data = data, p_treat = 0.5, estimator = "Hajek")
Hajek
CI_Hajek <- c(Hajek$point_est - 1.96 * sqrt(Hajek$var_est),
              Hajek$point_est + 1.96 * sqrt(Hajek$var_est))
CI_Hajek


#Test coverage: N = 100
set.seed(2024)
n_sim <- 10000
res_mat_n100 <- matrix(0, 3, 6)
rownames(res_mat_n100) <- c("HT", "Hajek", "Crude")
colnames(res_mat_n100) <- c("estimate", "95%CILB", "95%CIUB", "coverage rate", "rejection rate", "MSE")
for (i in 1:n_sim){
  n_samples <- 100
  
  data <- tibble(
    a = rbern(n_samples, .5),
    y = rbern(n_samples, .4),
    l = rbern(n_samples, .1 + .3*a + .5*y), #very few control healthy units selected so the effect is overestimated
    s = rbern(n_samples, .1 + .8*l)) #the larger l, the larger s)
  
  # IPCW estimator using information on L in unselected sample as well
  den.fit <- glm(s ~ l,
                 data = data, family = binomial)
  
  HT <- IPW_inference(data = data, p_treat = 0.5, estimator = "HT")
  res_mat_n100[1, 1] <- res_mat_n100[1, 1] + HT$point_est
  res_mat_n100[1, 2] <- res_mat_n100[1, 2] + HT$point_est - 1.96 * sqrt(HT$var_est)
  res_mat_n100[1, 3] <- res_mat_n100[1, 3] + HT$point_est + 1.96 * sqrt(HT$var_est)
  res_mat_n100[1, 4] <- res_mat_n100[1, 4] + 1 * (HT$point_est - 1.96 * sqrt(HT$var_est) <= true_effect ) * (HT$point_est + 1.96 * sqrt(HT$var_est) >= true_effect )
  res_mat_n100[1, 5] <- res_mat_n100[1, 5] + 1 * ((HT$point_est - 1.96 * sqrt(HT$var_est) > 0) | (HT$point_est + 1.96 * sqrt(HT$var_est) < 0))
  res_mat_n100[1, 6] <- res_mat_n100[1, 6] + (HT$point_est - true_effect)^2
    
  Hajek <- IPW_inference(data = data, p_treat = 0.5, estimator = "Hajek")
  res_mat_n100[2, 1] <- res_mat_n100[2, 1] + Hajek$point_est
  res_mat_n100[2, 2] <- res_mat_n100[2, 2] + Hajek$point_est - 1.96 * sqrt(Hajek$var_est)
  res_mat_n100[2, 3] <- res_mat_n100[2, 3] + Hajek$point_est + 1.96 * sqrt(Hajek$var_est)
  res_mat_n100[2, 4] <- res_mat_n100[2, 4] + 1 * (Hajek$point_est - 1.96 * sqrt(Hajek$var_est) <= true_effect) * (Hajek$point_est + 1.96 * sqrt(Hajek$var_est) >= true_effect )
  res_mat_n100[2, 5] <- res_mat_n100[2, 5] + 1 * ((Hajek$point_est - 1.96 * sqrt(Hajek$var_est) > 0) | (Hajek$point_est + 1.96 * sqrt(Hajek$var_est) < 0))
  res_mat_n100[2, 6] <- res_mat_n100[2, 6] + (Hajek$point_est - true_effect)^2
  
  # Filter data to include only selected observations (S = 1)
  selected_data <- data[data$s == 1, ]
  # Crude estimate in the selected sample 
  model_crude <- glm(y ~ a, data = selected_data)
  CI_crude <- confint(profile(model_crude))
  res_mat_n100[3, 1] <- res_mat_n100[3, 1] + model_crude$coefficients[2]
  res_mat_n100[3, 2] <- res_mat_n100[3, 2] + CI_crude[2,1]
  res_mat_n100[3, 3] <- res_mat_n100[3, 3] + CI_crude[2,2]
  res_mat_n100[3, 4] <- res_mat_n100[3, 4] + 1 * (CI_crude[2,1] <= true_effect) * (CI_crude[2,2] >= true_effect)
  res_mat_n100[3, 5] <- res_mat_n100[3, 5] + 1 * ((CI_crude[2,1] > 0) | (CI_crude[2,2] < 0))
  res_mat_n100[3, 6] <- res_mat_n100[3, 6] + (model_crude$coefficients[2] - true_effect)^2
  #summary(model_crude)
  print(i)
  
}
res_mat_n100 <- res_mat_n100/n_sim
res_mat_n100


#Test coverage: N = 1000
set.seed(2024)
n_sim <- 10000
res_mat_n1000 <- matrix(0, 3, 6)
rownames(res_mat_n1000) <- c("HT", "Hajek", "Crude")
colnames(res_mat_n1000) <- c("estimate", "95%CILB", "95%CIUB", "coverage rate", "rejection rate", "MSE")
for (i in 1:n_sim){
  n_samples <- 1000
  
  data <- tibble(
    a = rbern(n_samples, .5),
    y = rbern(n_samples, .4),
    l = rbern(n_samples, .1 + .3*a + .5*y), #very few control healthy units selected so the effect is overestimated
    s = rbern(n_samples, .1 + .8*l)) #the larger l, the larger s)
  
  # IPCW estimator using information on L in unselected sample as well
  den.fit <- glm(s ~ l,
                 data = data, family = binomial)
  
  HT <- IPW_inference(data = data, p_treat = 0.5, estimator = "HT")
  res_mat_n1000[1, 1] <- res_mat_n1000[1, 1] + HT$point_est
  res_mat_n1000[1, 2] <- res_mat_n1000[1, 2] + HT$point_est - 1.96 * sqrt(HT$var_est)
  res_mat_n1000[1, 3] <- res_mat_n1000[1, 3] + HT$point_est + 1.96 * sqrt(HT$var_est)
  res_mat_n1000[1, 4] <- res_mat_n1000[1, 4] + 1 * (HT$point_est - 1.96 * sqrt(HT$var_est) <= true_effect ) * (HT$point_est + 1.96 * sqrt(HT$var_est) >= true_effect )
  res_mat_n1000[1, 5] <- res_mat_n1000[1, 5] + 1 * ((HT$point_est - 1.96 * sqrt(HT$var_est) > 0) | (HT$point_est + 1.96 * sqrt(HT$var_est) < 0))
  res_mat_n1000[1, 6] <- res_mat_n1000[1, 6] + (HT$point_est - true_effect)^2
  
  Hajek <- IPW_inference(data = data, p_treat = 0.5, estimator = "Hajek")
  res_mat_n1000[2, 1] <- res_mat_n1000[2, 1] + Hajek$point_est
  res_mat_n1000[2, 2] <- res_mat_n1000[2, 2] + Hajek$point_est - 1.96 * sqrt(Hajek$var_est)
  res_mat_n1000[2, 3] <- res_mat_n1000[2, 3] + Hajek$point_est + 1.96 * sqrt(Hajek$var_est)
  res_mat_n1000[2, 4] <- res_mat_n1000[2, 4] + 1 * (Hajek$point_est - 1.96 * sqrt(Hajek$var_est) <= true_effect) * (Hajek$point_est + 1.96 * sqrt(Hajek$var_est) >= true_effect )
  res_mat_n1000[2, 5] <- res_mat_n1000[2, 5] + 1 * ((Hajek$point_est - 1.96 * sqrt(Hajek$var_est) > 0) | (Hajek$point_est + 1.96 * sqrt(Hajek$var_est) < 0))
  res_mat_n1000[2, 6] <- res_mat_n1000[2, 6] + (Hajek$point_est - true_effect)^2
  
  # Filter data to include only selected observations (S = 1)
  selected_data <- data[data$s == 1, ]
  # Crude estimate in the selected sample 
  model_crude <- glm(y ~ a, data = selected_data)
  CI_crude <- confint(profile(model_crude))
  res_mat_n1000[3, 1] <- res_mat_n1000[3, 1] + model_crude$coefficients[2]
  res_mat_n1000[3, 2] <- res_mat_n1000[3, 2] + CI_crude[2,1]
  res_mat_n1000[3, 3] <- res_mat_n1000[3, 3] + CI_crude[2,2]
  res_mat_n1000[3, 4] <- res_mat_n1000[3, 4] + 1 * (CI_crude[2,1] <= true_effect) * (CI_crude[2,2] >= true_effect)
  res_mat_n1000[3, 5] <- res_mat_n1000[3, 5] + 1 * ((CI_crude[2,1] > 0) | (CI_crude[2,2] < 0))
  res_mat_n1000[3, 6] <- res_mat_n1000[3, 6] + (model_crude$coefficients[2] - true_effect)^2
  #summary(model_crude)
  print(i)
  
}
res_mat_n1000 <- res_mat_n1000/n_sim
res_mat_n1000


#Test coverage: N = 10000
set.seed(2024)
n_sim <- 10000
res_mat_n10000 <- matrix(0, 3, 6)
rownames(res_mat_n10000) <- c("HT", "Hajek", "Crude")
colnames(res_mat_n10000) <- c("estimate", "95%CILB", "95%CIUB", "coverage rate", "rejection rate", "MSE")
for (i in 1:n_sim){
  n_samples <- 10000
  
  data <- tibble(
    a = rbern(n_samples, .5),
    y = rbern(n_samples, .4),
    l = rbern(n_samples, .1 + .3*a + .5*y), #very few control healthy units selected so the effect is overestimated
    s = rbern(n_samples, .1 + .8*l)) #the larger l, the larger s)
  
  # IPCW estimator using information on L in unselected sample as well
  den.fit <- glm(s ~ l,
                 data = data, family = binomial)
  
  HT <- IPW_inference(data = data, p_treat = 0.5, estimator = "HT")
  res_mat_n10000[1, 1] <- res_mat_n10000[1, 1] + HT$point_est
  res_mat_n10000[1, 2] <- res_mat_n10000[1, 2] + HT$point_est - 1.96 * sqrt(HT$var_est)
  res_mat_n10000[1, 3] <- res_mat_n10000[1, 3] + HT$point_est + 1.96 * sqrt(HT$var_est)
  res_mat_n10000[1, 4] <- res_mat_n10000[1, 4] + 1 * (HT$point_est - 1.96 * sqrt(HT$var_est) <= true_effect ) * (HT$point_est + 1.96 * sqrt(HT$var_est) >= true_effect )
  res_mat_n10000[1, 5] <- res_mat_n10000[1, 5] + 1 * ((HT$point_est - 1.96 * sqrt(HT$var_est) > 0) | (HT$point_est + 1.96 * sqrt(HT$var_est) < 0))
  res_mat_n10000[1, 6] <- res_mat_n10000[1, 6] + (HT$point_est - true_effect)^2
  
  Hajek <- IPW_inference(data = data, p_treat = 0.5, estimator = "Hajek")
  res_mat_n10000[2, 1] <- res_mat_n10000[2, 1] + Hajek$point_est
  res_mat_n10000[2, 2] <- res_mat_n10000[2, 2] + Hajek$point_est - 1.96 * sqrt(Hajek$var_est)
  res_mat_n10000[2, 3] <- res_mat_n10000[2, 3] + Hajek$point_est + 1.96 * sqrt(Hajek$var_est)
  res_mat_n10000[2, 4] <- res_mat_n10000[2, 4] + 1 * (Hajek$point_est - 1.96 * sqrt(Hajek$var_est) <= true_effect) * (Hajek$point_est + 1.96 * sqrt(Hajek$var_est) >= true_effect )
  res_mat_n10000[2, 5] <- res_mat_n10000[2, 5] + 1 * ((Hajek$point_est - 1.96 * sqrt(Hajek$var_est) > 0) | (Hajek$point_est + 1.96 * sqrt(Hajek$var_est) < 0))
  res_mat_n10000[2, 6] <- res_mat_n10000[2, 6] + (Hajek$point_est - true_effect)^2
  
  # Filter data to include only selected observations (S = 1)
  selected_data <- data[data$s == 1, ]
  # Crude estimate in the selected sample 
  model_crude <- glm(y ~ a, data = selected_data)
  CI_crude <- confint(profile(model_crude))
  res_mat_n10000[3, 1] <- res_mat_n10000[3, 1] + model_crude$coefficients[2]
  res_mat_n10000[3, 2] <- res_mat_n10000[3, 2] + CI_crude[2,1]
  res_mat_n10000[3, 3] <- res_mat_n10000[3, 3] + CI_crude[2,2]
  res_mat_n10000[3, 4] <- res_mat_n10000[3, 4] + 1 * (CI_crude[2,1] <= true_effect) * (CI_crude[2,2] >= true_effect)
  res_mat_n10000[3, 5] <- res_mat_n10000[3, 5] + 1 * ((CI_crude[2,1] > 0) | (CI_crude[2,2] < 0))
  res_mat_n10000[3, 6] <- res_mat_n10000[3, 6] + (model_crude$coefficients[2] - true_effect)^2
  #summary(model_crude)
  print(i)
  
}
res_mat_n10000 <- res_mat_n10000/n_sim
res_mat_n10000


res_df_c1 <- data.frame(rbind(res_mat_n100,
                              res_mat_n1000,
                              res_mat_n10000))
res_df_c1$estimator <- factor(rep(row.names(res_df_c1)[1:3], 3))
res_df_c1$sample_size <- rep(c("100", "1000", "10000"), each = 3)
row.names(res_df_c1) = c()
res_df_c1

library(ggplot2)
res_graph <- ggplot(res_df_c1, 
                    aes(x = sample_size, y = estimate, 
                        color = estimator)) + 
  geom_point(size = 3.0,
             position = position_dodge(width=0.2)) +
  geom_errorbar(aes(
    x = sample_size,
    ymin = X95.CILB,
    ymax = X95.CIUB,
    width = 0.005,
    color = estimator
  ),
  position = position_dodge(width=0.2)) +
  geom_hline(yintercept=0) +
  annotate("text", x="100", y=0, label="True ATE = 0", size=5, color="black", 
           hjust = 1.1, vjust = -0.9) +
  labs(
    title = "ATE and 95% confidence interval estimates (Case 1)",
    x = "sample size",
    y = "treatment effect"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 22),
        axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.position="none")
res_graph

MSE_graph <- ggplot(res_df_c1, 
                    aes(x = sample_size, y = MSE, 
                        color = estimator)) + 
  geom_point(size = 3.0,
             position = position_dodge(width=0.2)) +
  labs(
    title = "Mean squared error of the point estimate (Case 1)",
    x = "sample size",
    y = "mean squared error"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 22),
        axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.position="none")
MSE_graph

coverage_graph <- ggplot(res_df_c1, 
                         aes(x = sample_size, y = coverage.rate, 
                             color = estimator)) + 
  geom_point(size = 3.0,
             position = position_dodge(width=0.2)) +
  labs(
    title = "Coverage rate of the truth by the CI (Case 1)",
    x = "sample size",
    y = "coverage rate"
  ) +
  geom_hline(yintercept=0.95) +
  annotate("text", x="100", y=0.94, label="Confidence level = 0.95", size=5, color="black", 
           hjust = 0.6, vjust = -0.9) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 22),
        axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.position="none")
coverage_graph

error_graph <- ggplot(res_df_c1, 
                      aes(x = sample_size, y = rejection.rate, 
                          color = estimator)) + 
  geom_point(size = 3.0,
             position = position_dodge(width=0.2)) +
  labs(
    title = "Type I error rate (Case 1)",
    x = "sample size",
    y = "error rate"
  ) +
  geom_hline(yintercept=0.05) +
  annotate("text", x="100", y=-0.01, label="Significance level = 0.05", size=5, color="black", 
           hjust = 0.6, vjust = -0.7)  +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 22),
        axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20))
error_graph

#saveRDS(res_df_c1, "res_df_c1.rds")
#res_df_c1 <- readRDS(file = 'res_df_c1.rds')

library(ggpubr)
p <- ggarrange(res_graph, MSE_graph, coverage_graph, error_graph, ncol=2, nrow=2, 
          common.legend = TRUE, legend="bottom",
          font.label = list(size = 60))
#1600,1200

ggsave(
  "figure_simulation_results_case1.tiff",
  plot = p,
  width = 7,
  height = 7,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

library(gridExtra)
library(cowplot)
plot_grid(res_graph, MSE_graph, coverage_graph, error_graph, nrow = 2, ncol = 2)
#plot_grid(
#  plot_grid(res_graph, coverage_graph, nrow = 1, ncol = 2),
#  plot_grid(NULL, error_graph, NULL, nrow = 1, rel_widths = c(0.5, 1.5, 0.5)),
#  nrow = 2
#)
#1200,1000
#plot_grid(res_graph, coverage_graph, error_graph, align = "h", ncol = 2, rel_widths = c(5/16, 5/16, 6/16))
#grid.arrange(res_graph, coverage_graph, error_graph, ncol=3)

