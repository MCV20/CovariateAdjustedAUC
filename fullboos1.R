library(haven)
library(dplyr)
library(nnet)


df <- read_sas("/home/Projects/LungCancerAUC/cmnanal.sas7bdat", NULL)

df <- df %>% 
  transmute(
    A   = as.numeric(rndgroup == 1),
    Y   =ifelse(!is.na(cancyr) & cancyr <= 6, 1L, 0L),
    age = age,
    race      = case_when(
      ethnic == 1 ~ "Hispanic",
      race == 1 ~ "White",
      race == 2 ~ "Black",
      race == 3 ~ "Asian",
      race == 4 ~ "American Indian or Alaskan Native",
      race == 5  ~ "Native Hawaiian or Pacific Islander",
      .default = NA_character_
    ),
    education      = recode(
      as.character(educat),
      `1` = 1,
      `2` = 1,
      `3` = 2,
      `4` = 3,
      `5` = 4,
      `6` = 5,
      `7` = 6,
      .default = NA_real_
    ),
    BMI = 703 * weight / height^2,
    COPD = diagcopd,
    PHC = ifelse(cancblad == 1 | cancbrea == 1 | canccerv == 1 | canccolo == 1 |
                   cancesop == 1 | canckidn == 1 | canclung == 1 | cancnasa == 1 |
                   cancoral == 1 | cancpanc == 1 | cancphar == 1 | cancstom == 1 |
                   cancthyr == 1 | canctran == 1, 1, 0),
    FHC = ifelse(fambrother == 1 | famchild == 1 | famfather == 1 | fammother == 1 | famsister == 1, 1, 0),
    packyears = pkyr,
    cigsmok   = cigsmok,
    smokeday = smokeday,
    smokeyr = smokeyr,
    smoke_inconsistent = (cigsmok == 1 & !is.na(age_quit)),
    age_quit_clean = if_else(cigsmok == 1, NA_real_, age_quit),
    quit_time = case_when(
      cigsmok == 0 & !is.na(age_quit_clean) ~ pmax(0, age - age_quit_clean),
      cigsmok == 1 ~ 0,         
      TRUE ~ NA_real_             
    )
  ) %>%
  select(!age_quit_clean)

df <- df[complete.cases(df),]


## Take Coefficients from model (PAPER)
beta <- c(`(Intercept)` = -4.532506, 
          age = 0.0778868,
          raceBlack = 0.3944778,
          raceHispanic = -0.7434744, 
          raceAsian = -0.466585,
          `raceAmerican Indian or Alaskan Native` = 0,
          `raceNative Hawaiian or Pacific Islander`= 1.027152,
          education = -0.0812744,
          BMI =  -0.0274194,
          COPD = 0.3553063,
          PHC = 0.4589971,
          FHC = 0.587185,
          cigsmok = 0.2597431,
          smokeday = -1.822606,
          smokeyr = 0.0317321, 
          quit_time = -0.0308572)

df$race <- factor(
  df$race,
  levels = c(
    "White",
    "Black",
    "Hispanic",
    "Asian",
    "American Indian or Alaskan Native",
    "Native Hawaiian or Pacific Islander"
  )
)
df$smokeday <- (df$smokeday/10)^-1 - 0.4021541613
df$age <- df$age - 62
df$education <- df$education - 4
df$BMI <- df$BMI - 27
df$smokeyr <- df$smokeyr - 27
df$quit_time <- df$quit_time - 10

## Get Marker
X <- model.matrix( ~ age + race + education + BMI + COPD + PHC + FHC +
                     cigsmok + smokeday + smokeyr + quit_time, data = df )
logit_M <- X %*% beta
df$M <- logit_M


get_estimators_fast <- function(df, a) {
  ## Get indices for Y=1 and Y=0
  idx1 <- df$Y == 1; n1 <- sum(idx1)
  idx0 <- df$Y == 0; n0 <- sum(idx0)
  
  ## Get data for cases and controls
  df1 <- df[df$Y == 1, , drop = FALSE]
  df0 <- df[df$Y == 0, , drop = FALSE]
  
  
  
  ## Normal linear models for M | X, Y, A
  modM1 <- lm(M ~ A + age + education + BMI + COPD + PHC + FHC +
                cigsmok + smokeday + smokeyr + quit_time + packyears,
              data = df1)
  
  modM0 <- lm(M ~ A + age + education + BMI + COPD + PHC + FHC +
                cigsmok + smokeday + smokeyr + quit_time + packyears,
              data = df0)
  
  ## Evaluate everything at intervention arm A = a
  new_all <- df
  new_all$A <- a
  
  ## Get mu's
  mu1 <- predict(modM1, newdata = new_all, xlev = list(race = race_lvls))
  mu0 <- predict(modM0, newdata = new_all, xlev = list(race = race_lvls))
  
  ## Get s for Mi-Mj
  s1     <- summary(modM1)$sigma
  s0     <- summary(modM0)$sigma
  s <- sqrt(s1^2 + s0^2)
  
  ## Get the case and control mu's 
  mu1_cases <- mu1[idx1]
  mu0_controls <- mu0[idx0]
  
  ## Get standarized estimator
  num1 <- 0
  for (i in 1:sum(idx1)) {
    num1 <- num1 + sum(pnorm(mu1_cases[i] - mu0_controls, 0, s))
  }
  
  Est1 <- num1/(n1*n0)
  
  ## Now for covariate adjusted estimator we need:
  ## Here we decompose the sum as the sum accross all combinations - sum of equal i = j
  
  ## First get outcome model and probabilitys of case and control
  modY <- glm(Y ~ A + age + education + BMI + COPD + PHC + FHC +
                cigsmok + smokeday + smokeyr + quit_time + packyears,
              data = df,
              family = binomial())
  
  PY1 <- predict(modY, newdata = new_all, type = "response",
                 xlev = list(race = race_lvls))
  PY0 <- 1 - PY1
  df$PY1 <- PY1
  df$PY0 <- PY0
  
  ## Numerator m1(Xi,Xj)m2(Xi,Xj)
  
  ## sum_ij
  num11 <- 0
  for (i in 1:length(mu1)){
    num11 <- num11 + sum(pnorm(mu1[i] - mu0, 0, s) * PY1[i] * PY0)
  }
  
  ## sum_ii
  num12 <- sum( pnorm(mu1- mu0, 0, s) * PY1 * PY0)
  
  ## Numerator
  num1 <- num11 - num12 #this numerator will be used in the DR estimation
  
  
  ## Denominator m2(Xi,Xj)   
  #sumij - sumii
  den1 <- sum(PY1)*sum(PY0) - sum(PY1*PY0) ## This denominator will be used in the DR estimation
  
  ## Get Covariate Adjuste Estimator
  Est2 <- num1/den1
  
  
  ## Now get the robust estimator:
  #sum(m1XiXj*m2XiXj + ipw*(I_ind - m1XiXj*m2XiXj), na.rm = T) /
  #sum(m2XiXj + ipw*(I_10 - m2XiXj),na.rm = T)
  # We have sum(m1XiXj*m2XiXj ) and sum(m2XiXj) so we focus on the augmentation terms
  
  
  ## Get numerator for IPW part
  
  ## Get Propensity Score
  modA <- glm(A ~ 1, data = df, family = binomial())
  p1   <- predict(modA, type = "response")
  PA <- if (a == 1) p1 else 1 - p1
  
  ## Get indices for Y=1 and Y=0 on treatment arm A=a
  idx_case_a  <- which(df$A == a & df$Y == 1)
  idx_ctrl_a  <- which(df$A == a & df$Y == 0)
  
  
  ## The IPW term is zero when i=j so we focus on sum_ij
  ## Since we have indicators for Y=1 Y=0 on A=1 we use those indices
  numIPW <- 0
  denIPW <- 0
  for (i in 1:length(idx_case_a)) {
    numIPW <- numIPW +  sum(1/PA[idx_case_a[i]] * 1/PA[idx_ctrl_a] * as.numeric(df$M[idx_case_a[i]] > df$M[idx_ctrl_a]))
    denIPW <- denIPW +  sum(1/PA[idx_case_a[i]] * 1/PA[idx_ctrl_a])
    
  }
  
  ## Now AUgmentation Subctraction term
  ## Get numerator for substraction term
  idxA <-  which(df$A == a)
  
  numA1 <- 0
  denA1 <- 0
  for (i in 1:length(mu1[idxA])){
    numA1 <- numA1 + sum(pnorm(mu1[idxA][i] - mu0[idxA], 0, s) * PY1[idxA][i] * PY0[idxA] * 1/PA[idxA][i] * 1/PA[idxA])
    denA1 <- denA1 + sum(PY1[idxA][i] * PY0[idxA] * 1/PA[idxA][i] * 1/PA[idxA])
  }
  
  ## sum_ii
  numA2 <- sum( pnorm(mu1[idxA]- mu0[idxA], 0, s) * PY1[idxA] * PY0[idxA] * 1/PA[idxA] * 1/PA[idxA])
  denA2 <- sum(PY1[idxA] * PY0[idxA] * 1/PA[idxA] * 1/PA[idxA])
  
  
  
  numR <- num1 + numIPW - (numA1 - numA2)
  denR <- den1 + denIPW - (denA1 - denA2)
  
  
  EstR <- numR/denR
  
  
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
  
  Emp <- empiricalAUC_fast(a, df)
  
  return(c(
    Emp = Emp,
    Est1      = Est1,
    Est2      = Est2,
    EstDR     = EstR
  ))
}



## Bootstrap
n <- nrow(df)
B <- 1000
estimators1 <- matrix(NA_real_, ncol = 4, nrow = B,
                      dimnames = list(NULL, c("Emp", "Est1", "Est2", "EstDR")))
set.seed(1)
for (i in 1:B) {
  cat("Bootstrap", i, "of", B, "\n")
  ids <- sample(n, replace = TRUE)
  df_sample <- df[ids, , drop = FALSE]
  table(df_sample$race)
  estimators1[i, ] <- get_estimators_fast(df_sample, 1)
}

write.csv(estimators1, "/home/mcolonva/Boostrap/fullBoostrap1.csv", row.names = FALSE)

