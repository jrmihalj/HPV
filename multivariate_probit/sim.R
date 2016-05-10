# generating data from a multivariate probit model
library(clusterGeneration)
d <- 3 # number of dimensions (strains)

n_visits <- 15
n_patients <- 100
n <- n_patients * n_visits
patient <- rep(1:n_patients, times = n_visits)

nrep <- 1000
vars <- array(dim = c(d, d, nrep))

for (i in 1:nrep) {
  # we need variances for each strain in each of the two random effect levels
  # and they need to sum to one. Each row in `variance` corresponds to a strain
  # and each column is a ranef level
  variance <- array(dim = c(d, 2))
  variance[, 2] <- runif(d, .1, .9)
  variance[, 1] <- 1 - variance[, 2]
  stopifnot(all(apply(variance, 1, sum) == 1))

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

  eps1 <- array(rnorm(n * d), dim = c(n, d)) %*% L_1
  eps2 <- array(rnorm(n_patients * d), dim = c(n_patients, d)) %*% L_2
  patient <- rep(1:n_patients, times = n_visits)

  y_star <- eps1 + eps2[patient, ]

  vars[, , i] <- var(y_star)
}

# visualize the variances and correlations
par(mfrow = c(d, d))
for (i in 1:d) {
  for (j in 1:d) {
    hist(vars[i, j, ], breaks = 50,
         main = paste0('Cov(', i, ', ', j, ')'))
    if (i == j) {
      abline(v = 1, col = 2, lwd = 3, lty = 2)
    }
  }
}

pairs(y_star, col = patient)

# convert the last realization to a bernoulli r.v.
y <- ifelse(y_star > 0, 1, 0)
