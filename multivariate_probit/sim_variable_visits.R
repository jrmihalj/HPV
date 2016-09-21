# generating data from a multivariate probit model
library(clusterGeneration)

d <- 2 # number of dimensions (strains)
n_patients <- 250
n_visits_patient <- ceiling(-abs(rnorm(n_patients,0,2.5))) + 10 # most have all 10, but some less
n <- sum(n_visits_patient)
patient <- rep(1:n_patients, times=n_visits_patient)
visit <- NULL
for(i in 1:n_patients){
  temp <- c(1:n_visits_patient[i])
  visit <- c(visit, temp)
}
n_visit_max <- 10

# Calculate random effect variances, and make sure they sum to 1
variance <- array(dim = c(d, 2))
variance[, 2] <- runif(d, .1, .9)
variance[, 1] <- 1 - variance[, 2]

# construct correlation matrices
Rho_patient <- genPositiveDefMat(d, covMethod = 'onion',
                                 eta = 2, rangeVar = c(1, 1))$Sigma
Rho_visit <- genPositiveDefMat(d, covMethod = 'onion',
                               eta = 2, rangeVar = c(1, 1))$Sigma

# construct covariance matrices
sdev <- sqrt(variance)
Sigma_patient <- diag(sdev[, 1]) %*% Rho_patient %*% diag(sdev[, 1])
Sigma_visit <- diag(sdev[, 2]) %*% Rho_visit %*% diag(sdev[, 2])

L_1 <- chol(Sigma_visit)
L_2 <- chol(Sigma_patient)

# draw random effect values
eps1 <- array(rnorm(n * d), dim = c(n, d)) %*% L_1
eps2 <- array(rnorm(n_patients * d), dim = c(n_patients, d)) %*% L_2

# add up the nested random effects:
randefs <- eps1 + eps2[patient, ] 

# create a fixed effect of previous co-habitating species' occurrences:
betas_phi <- array(rnorm(d*d, 0, 1), dim=c(d,d))
betas_gam <- array(rnorm(d*d, 0, 1), dim=c(d,d))
diag(betas_gam) <- 0

# create a fixed effect of the time between visits (tbv) in units of days
# This covariate is strain specific, and can affect either phi or gamma
betas_tbv_phi <- rnorm(d, 0, 0.5)
betas_tbv_gam <- rnorm(d, 0, 0.5)

# generate the tbv covariate values
# I'll assume most at 14 days, some up to ~30
tbv <- floor(abs(rnorm(n * d, 0, 5))) + 14
# center, scale, and put into a matrix
tbv <- matrix(scale(tbv), nrow=n, ncol=d)
#tbv <- (tbv)^(4/7)
#hist(tbv)

# create species-specific intercepts for y_star
# assume these are low, because that is reflected in our data
alphas <- rnorm(d, -2, 1)   # transformed to unit normal scale, for the probit

# estimate latent variable, y*, and the observed occurrence, y:
y_star <- array(0, dim=c(n,d))
y <- array(0, dim=c(n,d))

# create the fixed effects
mu_all <- array(0, dim=c(n,d))

for(i in 1:n){
  
  if(visit[i] == 1){# For the initial visit, use only random effects
    
    y_star[i, ] <- randefs[i, ]
    mu_all[i, ] <- 0
    
  }else{
    for(j in 1:d){
      
      if(y[i - 1, j] == 1){ # If the previous observation for strain j in patient i was PRESENT, then persistence:
        mu_all[i,j] <- (y[i - 1, ] %*% betas_phi)[j] + tbv[i,j] * betas_tbv_phi[j]
      }else{ # If the previous observation for strain j in patient i was ABSENT, then colonization:
        mu_all[i,j] <- (y[i - 1, ] %*% betas_gam)[j] + tbv[i,j] * betas_tbv_gam[j]
      }
      
    }
    
    y_star[i, ] <- alphas + mu_all[i, ] + randefs[i, ]
    
  }
  
  y[i, ] <- ifelse(y_star[i, ] > 0, 1, 0)
  
}

# check overall prevalences
colSums(y) / nrow(y)

# # check correlations
# library(tidyverse)
# hist(y_star)
# df_test <- data.frame(y_star, patient, visit)
# 
# sub <- df_test %>%
#   filter(patient %in% sample(1:n_patients, size=10)) %>%
#   ggplot(aes(x=X2, y=X3))+
#   geom_point(shape=19, aes(color=factor(patient))) +
#   facet_wrap(~patient, nrow=2)
# print(sub)
# 
# cor(df_test[,1:2])
