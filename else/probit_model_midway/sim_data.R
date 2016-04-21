# generative model from sebastian gonzalez paper
library(clusterGeneration)
library(smoothmest)

time <- rep(1:n_timesteps, each = n_site)
site <- rep(1:n_site, n_timesteps)

# Patient-level corr matrix
# Rp <- genPositiveDefMat(dim=m, #Number of columns/rows
#                         covMethod = 'onion', 
#                         rangeVar = c(1, 1), # Range of variances
#                         eta=2)$Sigma

Rp <- read.csv(Rp_filename)[,-1]
# observation-level ranefs
ep <- matrix(nrow = n_site, ncol = m)
for (i in 1:n_site){
  ep[i, ] <- rnorm(m) %*% t(chol(Rp))
}

# observation-level corr matrix
# Patient-level corr matrix
Ro <- read.csv(Ro_filename)[,-1]

# observation-level ranefs
eo <- matrix(nrow = n, ncol = m)
for (i in 1:n){
  eo[i, ] <- rnorm(m) %*% t(chol(Ro))
}

# Add the ranef
e_all <- matrix(nrow = n, ncol = m)
for(i in 1:n){
  e_all[i, ] <- ep[site[i], ] + eo[i, ]
}

# Normalize
e_all <- (e_all - mean(e_all)) / sd(e_all)

# fixed effects design matrix, parameters, and linear predictor
# (note currently this is an intercept-only model, but X could be expanded)
k <- 1
X <- matrix(1, nrow = n, ncol = k)
beta <- matrix(rnorm(m * k), nrow = k)

# begin by simulating the state at time = 1
mu <- X %*% beta
z <- mu + qlogis(pnorm(e_all))
y <- matrix(nrow = n, ncol = m)
y[time == 1, ] <- ifelse(z[time == 1, ] > 0, 1, 0)

# then define species "interaction" parameters, and subsequent time series
beta_intxn <- matrix(rdoublex(m ^ 2, 0, 1), nrow = m)
for (i in 2:n_timesteps) {
  z[time == i, ] <- z[time == i, ] + y[time == i - 1, ] %*% beta_intxn
  y[time == i, ] <- ifelse(z[time == i, ] > 0, 1, 0)
}


# I'm going to assume that we begin at t = 2, assuming that we know the 
# state at time 1, so that we can use it as a covariate
# we need to construct the design matrix 
# (this was ambiguous in the original paper, but this seems defensible)
y_use <- y[time > 1, ]
X_use <- as.matrix(X[time > 1, ])
X_intxn <- y[time < n_timesteps, ]
