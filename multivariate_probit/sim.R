# generating data from a multivariate probit model
# with a varying number of visits among patients
library(clusterGeneration)


## high level data parameters --------------------------------------------------
n_strains <- 2
n_patients <- 250
n_visit_max <- 10


# generate the number of visits per patient, which is between 1 and 10
n_visits_patient <- rbinom(n_patients, n_visit_max, prob = .8)
if (min(n_visits_patient) == 0) {
  zeros <- n_visits_patient == 0
  n_visits_patient[zeros] <- 1
}
n <- sum(n_visits_patient)


# generate indices for patients and visits
patient <- rep(1:n_patients, times = n_visits_patient)
visit <- NULL
for (i in 1:n_patients) {
  visit <- c(visit, 1:n_visits_patient[i])
}
stopifnot(length(patient) == length(visit))
stopifnot(length(patient) == n)


## Generate parameter values -------------------------------------------------
# random effect variances (must sum to 1 bc. probit!)
variance <- array(dim = c(n_strains, 2))
variance[, 2] <- runif(n_strains, .1, .9)
variance[, 1] <- 1 - variance[, 2]

# correlation matrices
Rho_patient <- genPositiveDefMat(n_strains, covMethod = 'onion',
                                 eta = 2, rangeVar = c(1, 1))$Sigma
Rho_visit <- genPositiveDefMat(n_strains, covMethod = 'onion',
                               eta = 2, rangeVar = c(1, 1))$Sigma

# covariance matrices
sdev <- sqrt(variance)
Sigma_patient <- diag(sdev[, 1]) %*% Rho_patient %*% diag(sdev[, 1])
Sigma_visit <- diag(sdev[, 2]) %*% Rho_visit %*% diag(sdev[, 2])

L_1 <- chol(Sigma_visit)
L_2 <- chol(Sigma_patient)

# draw random effect values
eps1 <- array(rnorm(n * n_strains), dim = c(n, n_strains)) %*% L_1
eps2 <- array(rnorm(n_patients * n_strains), dim = c(n_patients, n_strains)) %*% L_2

# add up the nested random effects by matching patients to observations
randefs <- eps1 + eps2[patient, ]

# create a fixed effect of previous co-habitating strain occurrences:
betas_phi <- matrix(rnorm(n_strains**2, 0, 1), nrow = n_strains) # persistence
betas_gam <- matrix(rnorm(n_strains**2, 0, 1), nrow = n_strains) # colonization
diag(betas_gam) <- 0 # 

# create a fixed effect of the time between visits (tbv) in units of days
# This covariate is strain specific, and can affect either phi or gamma
betas_tbv_phi <- rnorm(n_strains, 0, 0.5)
betas_tbv_gam <- rnorm(n_strains, 0, 0.5)

# generate the tbv covariate values
# I'll assume most at 14 days, some up to ~30
tbv <- floor(abs(rnorm(n, 0, 5))) + 14
# center, scale
tbv <- as.vector(scale(tbv))

# create species-specific intercepts for y_star
# assume these are low, because that is reflected in our data
alphas <- rnorm(n_strains, -2, 1)   # transformed to unit normal scale, for the probit




## Generate observation histories ----------------------------------------------
y_star <- array(dim = c(n, n_strains)) # latent (continuous)
y <- array(dim = c(n, n_strains)) # observed (0 or 1)

# fixed effects array
mu_all <- array(dim = c(n, n_strains))

for (i in 1:n) {
  if (visit[i] == 1) {# For the initial visit, use random effects + intercept
    y_star[i, ] <- randefs[i, ] + alphas
    mu_all[i, ] <- 0
  } else {
    for (j in 1:n_strains) {
      if (y[i - 1, j] == 1) { # strain j in patient i was present: persistence
        mu_all[i,j] <- (y[i - 1, ] %*% betas_phi)[j] + 
          tbv[i] * betas_tbv_phi[j]
      } else {# strain j in patient i was absent: colonization
        mu_all[i,j] <- (y[i - 1, ] %*% betas_gam)[j] + 
          tbv[i] * betas_tbv_gam[j]
      }
    }
    y_star[i, ] <- alphas + mu_all[i, ] + randefs[i, ]
  }
  y[i, ] <- ifelse(y_star[i, ] > 0, 1, 0)
}

# check overall prevalences
apply(y, 2, mean)

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
