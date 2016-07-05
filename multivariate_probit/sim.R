# generating data from a multivariate probit model
library(clusterGeneration)
d <- 3 # number of dimensions (strains)

n_visits <- 10
n_patients <- 150
n <- n_patients * n_visits
patient <- rep(1:n_patients, times = n_visits)
visit <- rep(1:n_visits, each=n_patients)

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

# draw randome effect values
eps1 <- array(rnorm(n * d), dim = c(n, d)) %*% L_1
eps2 <- array(rnorm(n_patients * d), dim = c(n_patients, d)) %*% L_2

# create a fixed effect of previous co-habitating species' occurrences:
betas_phi <- array(rnorm(d*d, 0, 1), dim=c(d,d))
betas_gam <- array(rnorm(d*d, 0, 1), dim=c(d,d))
diag(betas_gam) <- 0

# add up the nested random effects:
randefs <- eps1 + eps2[patient, ] 

# create species-specific intercepts for y_star
base_prob <- runif(d, 0, 1) # base-line probability of occurrence
alphas <- qnorm(base_prob)   # transformed to unit normal scale, for the probit

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
      
      if(y[i - n_patients, j] == 1){
        mu_all[i,j] <- (y[i - n_patients, ] %*% betas_phi)[j]
      }else{
        mu_all[i,j] <- (y[i - n_patients, ] %*% betas_gam)[j]
      }

    }
    
    y_star[i, ] <- alphas + mu_all[i, ] + randefs[i, ]
    
  }
  
  y[i, ] <- ifelse(y_star[i, ] > 0, 1, 0)
  
}

# check overall prevalences
colSums(y) / nrow(y)


