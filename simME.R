## Load Libraries 
library(ggplot2)
library(dplyr)
library(MASS)
library(tidyr)
library(kableExtra)
library(pbapply)    
library(scales)
theme_set(theme_minimal())


## Data Generating Function
sim_df <- function(n, muX, sigmaX, alpha0, alphaX, alphaA,
                   beta0, betaY, betaX, beta_X2, betaA, betaXY, betaAY, sigma0, sigma1){
  ## Simulate X
  X <- mvrnorm(n, muX, sigmaX)
  
  
  ## Meassurement Error
  eps <- matrix(rnorm(n * ncol(X), mean = 0.2, sd = 0.1), nrow = n, ncol = ncol(X))
  X <- X + eps
  
  
  ## Simulate A
  A <- rbinom(n, 1, 0.5)
  
  ## Simulate Y depending on A and X
  lpY <- alpha0 + as.numeric(alphaX %*% t(X)) + alphaA *A
  p <- exp(lpY)/(1+exp(lpY))
  Y <- rbinom(n,1,p)
  
  ## Simulate M such that is normal coniditoned on X,Y,A
  eps1 <- rnorm(n,0,sigma1)
  eps0 <- rnorm(n,0,sigma0)
  
  M <- beta0 + betaY*Y + as.numeric(betaX %*% t(X)) + beta_X2*X[,1]^2 + 
    betaA*A + as.numeric(betaXY %*% t(X))* Y + betaAY*Y*A + eps1*Y + eps0*(1-Y)
  M <- as.numeric(M)
  df <- data.frame(Y = Y, A = A, M = M, X = X)
  par <- data.frame(n = n, muX = muX, sigmaX=sigmaX, alpha0=alpha0, 
                    alphaX=alphaX, alphaA=alphaA, beta0=beta0, betaY=betaY, 
                    betaX=betaX, beta_X2=beta_X2, betaA=betaA, betaXY=betaXY, betaAY=betaAY, 
                    sigma0=sigma0, sigma1=sigma1) 
  return(list(df, par))
}


## Set Parameters
muX = c(2,-1,1); sigmaX = diag(c(0.3, 0.3, 0.3), nrow = 3, ncol = 3)
alpha0 = -1 ; alphaX = c(-0.5,-1,1); beta_X2 = 1; alphaA = 1 ;
beta0 = 1; betaY = 5; betaX = c(2,3,2); betaA = 0.5; betaXY = c(3 ,3.5, 2); betaAY = 4
sigma0 <- 0.25; sigma1 <- 0.25


## Set n: 200, 500, 700, 1000, 2500, 5000
n <- 1000

## Empirical AUC function 
empiricalAUC_fast <- function(a, df) {
  idx <- df$A == a
  y   <- df$Y[idx]
  m   <- df$M[idx]
  n1  <- sum(y == 1L)
  n0  <- sum(y == 0L)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r   <- rank(m, ties.method = "average")
  s1  <- sum(r[y == 1L])
  (s1 - n1 * (n1 + 1) / 2) / (n1 * n0)
}


## True AUC's
trueAUC1 <- 0.9496
trueAUC0 <- 0.8613


####################### Correctly Specified Scenario ###########################

get_estimators <- function(a, df){
  ## 1) Fit a linear regression accross Y=1 and Y=0 separately. 
  df1 <- df %>% filter(Y==1)
  df0 <- df %>% filter(Y==0)
  
  modM1 <- lm(M ~ A + X.1 + X.2 + X.3 + I(X.1^2), data = df1)
  summary(modM1)
  
  modM0 <- lm(M ~ A + X.1 + X.2 + X.3 + I(X.1^2), data = df0)
  summary(modM0)
  
  ## For A=1
  new_all <- transform(df, A = a)         
  mu1 <- predict(modM1, new_all)           
  mu0 <- predict(modM0, new_all)          
  
  df$mu1 <- mu1
  df$mu0 <- mu0
  
  #Get M1XiXj
  m1XiXj <- pnorm(outer(df$mu1, df$mu0, "-")/sqrt(sigma(modM1)^2 + sigma(modM0)^2))
  diag(m1XiXj) <- NA_real_
  
  
  #2) Fit outcome model
  modY <- glm(Y ~ A + X.1 + X.2 + X.3, df, family = "binomial")
  
  ## Get M2XiXj
  df$PY1 <- predict(modY, new_all, type = "response")
  df$PY0 <- 1-df$PY1
  m2XiXj <- outer(df$PY1, df$PY0, "*")
  diag(m2XiXj) <- NA_real_
  
  
  #3) Fit propensity Score Model
  modA <- glm(A ~ 1, df, family = "binomial")
  p1 <- predict(modA, df, type = "response")
  df$PA <- if(a==1) p1 else 1-p1
  ipw <- outer(as.numeric(df$A==a)/df$PA, as.numeric(df$A==a)/df$PA, "*")
  diag(ipw) <- NA_real_
  
  #4) This is extra for the DR estimator 
  I_M   <- outer(df$M, df$M, FUN = ">")          # 1{Mi > Mj}
  I_10  <- outer(df$Y == 1, df$Y == 0, FUN = "&")# 1{Yi=1, Yj=0}
  I_ind <-  I_M & I_10
  diag(I_ind) <- NA_real_
  diag(I_10) <- NA_real_
  
  ## Standarized
  Est1 <-  sum(m1XiXj*I_10, na.rm = T)/sum(I_10, na.rm = T)
  
  ## Covariate Adjusted
  Est2 <- sum(m2XiXj*m1XiXj, na.rm = T)/sum(m2XiXj, na.rm = T)
  
  ## DR Estimator
  EstDR <- sum(m1XiXj*m2XiXj + ipw*(I_ind - m1XiXj*m2XiXj), na.rm = T)/sum(m2XiXj + ipw*(I_10 - m2XiXj),na.rm = T)
  
  ## Empirical 
  Emp <- empiricalAUC_fast(a, df)
  
  return(c(Estimator1 = Est1, Estimator2 = Est2, Robust = EstDR, Empirical = Emp))
}


set.seed(2)
N_sim <- 5000
EstimatorsA1 <- matrix(NA, nrow = N_sim, ncol = 4)
EstimatorsA0 <- matrix(NA, nrow = N_sim, ncol = 4)
colnames(EstimatorsA1) <- c("Estimator 1", "Estimator 2", "Robust", "Empirical")
colnames(EstimatorsA0) <- c("Estimator 1", "Estimator 2", "Robust", "Empirical")
for (i in 1:N_sim) {
  sim <- sim_df(n, muX, sigmaX, alpha0, alphaX, alphaA,
                beta0, betaY, betaX, beta_X2, betaA, betaXY, betaAY, sigma0, sigma1)
  df <- sim[[1]]
  par <- sim[[2]]
  EstimatorsA1[i, ] <- get_estimators(1, df)   
  EstimatorsA0[i, ] <- get_estimators(0, df) 
  #print(i)
}

print("forloop both correct done")

EstimatorsA1 <- as.data.frame(EstimatorsA1) %>% mutate(A = "A=1")
EstimatorsA0 <- as.data.frame(EstimatorsA0) %>% mutate(A = "A=0")
Estimators <- rbind(EstimatorsA0, EstimatorsA1)
Estimators <- Estimators %>%
  mutate(A = as.integer(sub("^A=", "", A))) %>%          # A=0 -> 0, etc.
  pivot_longer(
    cols = c( `Estimator 1`,`Estimator 2`, Robust, Empirical),
    names_to = "estimator",
    values_to = "estimate"
  )

print("both correct estimators done")

rel_eff_var1 <- Estimators %>%
  group_by(A, estimator) %>%
  summarise(var = var(estimate, na.rm = TRUE), .groups = "drop") %>%
  group_by(A) %>%
  mutate(var_unadj = var[estimator == "Empirical"],
         RE = var_unadj / var) %>%              # >1 ⇒ more efficient than Empirical
  filter(estimator != "Empirical") %>%
  dplyr::select(A, estimator, RE)
rel_eff_var1$scenario <- "both correct"

biasA1 <- EstimatorsA1[,1:4] - trueAUC1
biasA0 <- EstimatorsA0[,1:4] - trueAUC0


biasA1 <- as.data.frame(biasA1) %>% mutate(A = "A=1")
biasA0 <- as.data.frame(biasA0) %>% mutate(A = "A=0")
bias <- rbind(biasA0, biasA1) %>%
  mutate(A = as.integer(sub("^A=", "", A))) %>%         
  pivot_longer(
    cols = c( `Estimator 1`,`Estimator 2`, Robust, Empirical),
    names_to = "estimator",
    values_to = "bias"
  )

bias1 <- bias %>% 
  group_by(A, estimator) %>%
  summarize(Bias = round(mean(bias), 4))
bias1$scenario <- "both correct"


print("bothcorrect finished")


################ Correct AUC model but incorrect outcome model #################

get_estimators <- function(a, df){
  ## 1) Fit a linear regression accross Y=1 and Y=0 separately. 
  df1 <- df %>% filter(Y==1)
  df0 <- df %>% filter(Y==0)
  
  modM1 <- lm(M ~ A + X.1 + X.2 + X.3 + I(X.1^2), data = df1)
  summary(modM1)
  
  modM0 <- lm(M ~ A + X.1 + X.2 + X.3 + I(X.1^2), data = df0)
  summary(modM0)
  
  ## For A=1
  new_all <- transform(df, A = a)         
  mu1 <- predict(modM1, new_all)           
  mu0 <- predict(modM0, new_all)          
  
  df$mu1 <- mu1
  df$mu0 <- mu0
  
  #Get M1XiXj
  m1XiXj <- pnorm(outer(df$mu1, df$mu0, "-")/sqrt(sigma(modM1)^2 + sigma(modM0)^2))
  diag(m1XiXj) <- NA_real_
  
  
  #2) Fit outcome model
  modY <- glm(Y ~ A + X.1 + X.2, df, family = "binomial")
  
  ## Get M2XiXj
  df$PY1 <- predict(modY, new_all, type = "response")
  df$PY0 <- 1-df$PY1
  m2XiXj <- outer(df$PY1, df$PY0, "*")
  diag(m2XiXj) <- NA_real_
  
  
  #3) Fit propensity Score Model
  modA <- glm(A ~ 1, df, family = "binomial")
  p1 <- predict(modA, df, type = "response")
  df$PA <- if(a==1) p1 else 1-p1
  ipw <- outer(as.numeric(df$A==a)/df$PA, as.numeric(df$A==a)/df$PA, "*")
  diag(ipw) <- NA_real_
  
  #4) This is extra for the DR estimator 
  I_M   <- outer(df$M, df$M, FUN = ">")          # 1{Mi > Mj}
  I_10  <- outer(df$Y == 1, df$Y == 0, FUN = "&")# 1{Yi=1, Yj=0}
  I_ind <-  I_M & I_10
  diag(I_ind) <- NA_real_
  diag(I_10) <- NA_real_
  
  ## Standarized
  Est1 <-  sum(m1XiXj*I_10, na.rm = T)/sum(I_10, na.rm = T)
  
  ## Estimator 2
  Est2 <- sum(m2XiXj*m1XiXj, na.rm = T)/sum(m2XiXj, na.rm = T)
  
  ## DR Estimator
  EstDR <- sum(m1XiXj*m2XiXj + ipw*(I_ind - m1XiXj*m2XiXj), na.rm = T)/sum(m2XiXj + ipw*(I_10 - m2XiXj),na.rm = T)
  
  ## Empirical 
  Emp <- empiricalAUC_fast(a, df)
  
  return(c(Estimator1 = Est1, Estimator2 = Est2, Robust = EstDR, Empirical = Emp))
}

set.seed(2)
EstimatorsA1 <- matrix(NA, nrow = N_sim, ncol = 4)
EstimatorsA0 <- matrix(NA, nrow = N_sim, ncol = 4)
colnames(EstimatorsA1) <- c("Estimator 1", "Estimator 2", "Robust", "Empirical")
colnames(EstimatorsA0) <- c("Estimator 1", "Estimator 2", "Robust", "Empirical")
for (i in 1:N_sim) {
  #print(i)
  sim <- sim_df(n, muX, sigmaX, alpha0, alphaX, alphaA,
                beta0, betaY, betaX, beta_X2, betaA, betaXY, betaAY, sigma0, sigma1)
  df <- sim[[1]]
  par <- sim[[2]]
  EstimatorsA1[i, ] <- get_estimators(1, df)   
  EstimatorsA0[i, ] <- get_estimators(0, df) 
}

EstimatorsA1 <- as.data.frame(EstimatorsA1) %>% mutate(A = "A=1")
EstimatorsA0 <- as.data.frame(EstimatorsA0) %>% mutate(A = "A=0")
Estimators <- rbind(EstimatorsA0, EstimatorsA1)
Estimators <- Estimators %>%
  mutate(A = as.integer(sub("^A=", "", A))) %>%          # A=0 -> 0, etc.
  pivot_longer(
    cols = c( `Estimator 1`,`Estimator 2`, Robust, Empirical),
    names_to = "estimator",
    values_to = "estimate"
  )


rel_eff_var2 <- Estimators %>%
  group_by(A, estimator) %>%
  summarise(var = var(estimate, na.rm = TRUE), .groups = "drop") %>%
  group_by(A) %>%
  mutate(var_unadj = var[estimator == "Empirical"],
         RE = var_unadj / var) %>%              # >1 ⇒ more efficient than Empirical
  filter(estimator != "Empirical") %>%
  dplyr::select(A, estimator, RE)
rel_eff_var2$scenario <- "M1 correct"

biasA1 <- EstimatorsA1[,1:4] - trueAUC1
biasA0 <- EstimatorsA0[,1:4] - trueAUC0


biasA1 <- as.data.frame(biasA1) %>% mutate(A = "A=1")
biasA0 <- as.data.frame(biasA0) %>% mutate(A = "A=0")
bias <- rbind(biasA0, biasA1) %>%
  mutate(A = as.integer(sub("^A=", "", A))) %>%         
  pivot_longer(
    cols = c( `Estimator 1`,`Estimator 2`, Robust, Empirical),
    names_to = "estimator",
    values_to = "bias"
  )

bias2 <- bias %>% 
  group_by(A, estimator) %>%
  summarize(Bias = round(mean(bias), 4))
bias2$scenario <- "M1 correct"


################ Incorrect AUC model but correct outcome model #################

get_estimators <- function(a, df){
  ## 1) Fit a linear regression accross Y=1 and Y=0 separately. 
  df1 <- df %>% filter(Y==1)
  df0 <- df %>% filter(Y==0)
  
  modM1 <- lm(M ~ A + X.1 + X.2 + X.3, data = df1)
  summary(modM1)
  
  modM0 <- lm(M ~ A + X.1 + X.2 + X.3, data = df0)
  summary(modM0)
  
  ## For A=1
  new_all <- transform(df, A = a)         
  mu1 <- predict(modM1, new_all)           
  mu0 <- predict(modM0, new_all)          
  
  df$mu1 <- mu1
  df$mu0 <- mu0
  
  #Get M1XiXj
  m1XiXj <- pnorm(outer(df$mu1, df$mu0, "-")/sqrt(sigma(modM1)^2 + sigma(modM0)^2))
  diag(m1XiXj) <- NA_real_
  
  
  #2) Fit outcome model
  modY <- glm(Y ~ A + X.1 + X.2 + X.3, df, family = "binomial")
  
  ## Get M2XiXj
  df$PY1 <- predict(modY, new_all, type = "response")
  df$PY0 <- 1-df$PY1
  m2XiXj <- outer(df$PY1, df$PY0, "*")
  diag(m2XiXj) <- NA_real_
  
  
  #3) Fit propensity Score Model
  modA <- glm(A ~ 1, df, family = "binomial")
  p1 <- predict(modA, df, type = "response")
  df$PA <- if(a==1) p1 else 1-p1
  ipw <- outer(as.numeric(df$A==a)/df$PA, as.numeric(df$A==a)/df$PA, "*")
  diag(ipw) <- NA_real_
  
  #4) This is extra for the DR estimator 
  I_M   <- outer(df$M, df$M, FUN = ">")          # 1{Mi > Mj}
  I_10  <- outer(df$Y == 1, df$Y == 0, FUN = "&")# 1{Yi=1, Yj=0}
  I_ind <-  I_M & I_10
  diag(I_ind) <- NA_real_
  diag(I_10) <- NA_real_
  
  
  ## Standarized
  Est1 <-  sum(m1XiXj*I_10, na.rm = T)/sum(I_10, na.rm = T)
  
  ## Estimator 2
  Est2 <- sum(m2XiXj*m1XiXj, na.rm = T)/sum(m2XiXj, na.rm = T)
  
  ## DR Estimator
  EstDR <- sum(m1XiXj*m2XiXj + ipw*(I_ind - m1XiXj*m2XiXj), na.rm = T)/sum(m2XiXj + ipw*(I_10 - m2XiXj),na.rm = T)
  
  ## Empirical 
  Emp <- empiricalAUC_fast(a, df)
  
  return(c(Estimator1 = Est1, Estimator2 = Est2, Robust = EstDR, Empirical = Emp))
}

set.seed(2)
EstimatorsA1 <- matrix(NA, nrow = N_sim, ncol = 4)
EstimatorsA0 <- matrix(NA, nrow = N_sim, ncol = 4)
colnames(EstimatorsA1) <- c("Estimator 1", "Estimator 2", "Robust", "Empirical")
colnames(EstimatorsA0) <- c("Estimator 1", "Estimator 2", "Robust", "Empirical")
for (i in 1:N_sim) {
  sim <- sim_df(n, muX, sigmaX, alpha0, alphaX, alphaA,
                beta0, betaY, betaX, beta_X2, betaA, betaXY, betaAY, sigma0, sigma1)
  df <- sim[[1]]
  par <- sim[[2]]
  EstimatorsA1[i, ] <- get_estimators(1, df)   
  EstimatorsA0[i, ] <- get_estimators(0, df) 
  #print(i)
}

EstimatorsA1 <- as.data.frame(EstimatorsA1) %>% mutate(A = "A=1")
EstimatorsA0 <- as.data.frame(EstimatorsA0) %>% mutate(A = "A=0")
Estimators <- rbind(EstimatorsA0, EstimatorsA1)
Estimators <- Estimators %>%
  mutate(A = as.integer(sub("^A=", "", A))) %>%          # A=0 -> 0, etc.
  pivot_longer(
    cols = c( `Estimator 1`,`Estimator 2`, Robust, Empirical),
    names_to = "estimator",
    values_to = "estimate"
  )


rel_eff_var3 <- Estimators %>%
  group_by(A, estimator) %>%
  summarise(var = var(estimate, na.rm = TRUE), .groups = "drop") %>%
  group_by(A) %>%
  mutate(var_unadj = var[estimator == "Empirical"],
         RE = var_unadj / var) %>%              # >1 ⇒ more efficient than Empirical
  filter(estimator != "Empirical") %>%
  dplyr::select(A, estimator, RE)
rel_eff_var3$scenario <- "M2 correct"

biasA1 <- EstimatorsA1[,1:4] - trueAUC1
biasA0 <- EstimatorsA0[,1:4] - trueAUC0


biasA1 <- as.data.frame(biasA1) %>% mutate(A = "A=1")
biasA0 <- as.data.frame(biasA0) %>% mutate(A = "A=0")
bias <- rbind(biasA0, biasA1) %>%
  mutate(A = as.integer(sub("^A=", "", A))) %>%         
  pivot_longer(
    cols = c( `Estimator 1`,`Estimator 2`, Robust, Empirical),
    names_to = "estimator",
    values_to = "bias"
  )

bias3 <- bias %>% 
  group_by(A, estimator) %>%
  summarize(Bias = round(mean(bias), 4))
bias3$scenario <- "M2 correct"


################# Both AUC model and outcome model icnorrect ###################

get_estimators <- function(a, df){
  ## 1) Fit a linear regression accross Y=1 and Y=0 separately. 
  df1 <- df %>% filter(Y==1)
  df0 <- df %>% filter(Y==0)
  
  modM1 <- lm(M ~ A + X.1 + X.2 + X.3, data = df1)
  summary(modM1)
  
  modM0 <- lm(M ~ A + X.1 + X.2 + X.3, data = df0)
  summary(modM0)
  
  ## For A=1
  new_all <- transform(df, A = a)         
  mu1 <- predict(modM1, new_all)           
  mu0 <- predict(modM0, new_all)          
  
  df$mu1 <- mu1
  df$mu0 <- mu0
  
  #Get M1XiXj
  m1XiXj <- pnorm(outer(df$mu1, df$mu0, "-")/sqrt(sigma(modM1)^2 + sigma(modM0)^2))
  diag(m1XiXj) <- NA_real_
  
  
  #2) Fit outcome model
  modY <- glm(Y ~ A + X.1 + X.2, df, family = "binomial")
  
  ## Get M2XiXj
  df$PY1 <- predict(modY, new_all, type = "response")
  df$PY0 <- 1-df$PY1
  m2XiXj <- outer(df$PY1, df$PY0, "*")
  diag(m2XiXj) <- NA_real_
  
  
  #3) Fit propensity Score Model
  modA <- glm(A ~ 1, df, family = "binomial")
  p1 <- predict(modA, df, type = "response")
  df$PA <- if(a==1) p1 else 1-p1
  ipw <- outer(as.numeric(df$A==a)/df$PA, as.numeric(df$A==a)/df$PA, "*")
  diag(ipw) <- NA_real_
  
  #4) This is extra for the DR estimator 
  I_M   <- outer(df$M, df$M, FUN = ">")          # 1{Mi > Mj}
  I_10  <- outer(df$Y == 1, df$Y == 0, FUN = "&")# 1{Yi=1, Yj=0}
  I_ind <-  I_M & I_10
  diag(I_ind) <- NA_real_
  diag(I_10) <- NA_real_
  
  ## Standarized
  Est1 <-  sum(m1XiXj*I_10, na.rm = T)/sum(I_10, na.rm = T)
  
  ## Estimator 2
  Est2 <- sum(m2XiXj*m1XiXj, na.rm = T)/sum(m2XiXj, na.rm = T)
  
  ## DR Estimator
  EstDR <- sum(m1XiXj*m2XiXj + ipw*(I_ind - m1XiXj*m2XiXj), na.rm = T)/sum(m2XiXj + ipw*(I_10 - m2XiXj),na.rm = T)
  
  ## Empirical 
  Emp <- empiricalAUC_fast(a, df)
  
  return(c(Estimator1 = Est1, Estimator2 = Est2, Robust = EstDR, Empirical = Emp))
}


set.seed(2)
EstimatorsA1 <- matrix(NA, nrow = N_sim, ncol = 4)
EstimatorsA0 <- matrix(NA, nrow = N_sim, ncol = 4)
colnames(EstimatorsA1) <- c("Estimator 1", "Estimator 2", "Robust", "Empirical")
colnames(EstimatorsA0) <- c("Estimator 1", "Estimator 2", "Robust", "Empirical")
for (i in 1:N_sim) {
  sim <- sim_df(n, muX, sigmaX, alpha0, alphaX, alphaA,
                beta0, betaY, betaX, beta_X2, betaA, betaXY, betaAY, sigma0, sigma1)
  df <- sim[[1]]
  par <- sim[[2]]
  EstimatorsA1[i, ] <- get_estimators(1, df)   
  EstimatorsA0[i, ] <- get_estimators(0, df) 
  #print(i)
}

EstimatorsA1 <- as.data.frame(EstimatorsA1) %>% mutate(A = "A=1")
EstimatorsA0 <- as.data.frame(EstimatorsA0) %>% mutate(A = "A=0")
Estimators <- rbind(EstimatorsA0, EstimatorsA1)
Estimators <- Estimators %>%
  mutate(A = as.integer(sub("^A=", "", A))) %>%          # A=0 -> 0, etc.
  pivot_longer(
    cols = c( `Estimator 1`,`Estimator 2`, Robust, Empirical),
    names_to = "estimator",
    values_to = "estimate"
  )


rel_eff_var4 <- Estimators %>%
  group_by(A, estimator) %>%
  summarise(var = var(estimate, na.rm = TRUE), .groups = "drop") %>%
  group_by(A) %>%
  mutate(var_unadj = var[estimator == "Empirical"],
         RE = var_unadj / var) %>%              # >1 ⇒ more efficient than Empirical
  filter(estimator != "Empirical") %>%
  dplyr::select(A, estimator, RE)
rel_eff_var4$scenario <- "both incorrect"

biasA1 <- EstimatorsA1[,1:4] - trueAUC1
biasA0 <- EstimatorsA0[,1:4] - trueAUC0


biasA1 <- as.data.frame(biasA1) %>% mutate(A = "A=1")
biasA0 <- as.data.frame(biasA0) %>% mutate(A = "A=0")
bias <- rbind(biasA0, biasA1) %>%
  mutate(A = as.integer(sub("^A=", "", A))) %>%         
  pivot_longer(
    cols = c( `Estimator 1`,`Estimator 2`, Robust, Empirical),
    names_to = "estimator",
    values_to = "bias"
  )

bias4 <- bias %>% 
  group_by(A, estimator) %>%
  summarize(Bias = round(mean(bias), 4))
bias4$scenario <- "both incorrect"



## Join RESULTS:
rel_eff_var1$n <- n
rel_eff_var2$n <- n
rel_eff_var3$n <- n
rel_eff_var4$n <- n
rel_eff_var_n1 <- rbind(rel_eff_var1, rel_eff_var2, rel_eff_var3, rel_eff_var4)

bias1$n <- n
bias2$n <- n
bias3$n <- n
bias4$n <- n
bias_n1 <- rbind(bias1, bias2, bias3, bias4)

getwd()
write.csv(rel_eff_var_n1, "/home/mcolonva/AUC_HS/relative_efficienty4.csv")
write.csv(bias_n1, "/home/mcolonva/AUC_HS/bias4.csv")


