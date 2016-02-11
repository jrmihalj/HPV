# generative model from sebastian gonzalez paper
library(clusterGeneration)
library(smoothmest)


# define parameters -------------------------------------
m <- 3                    # species
n_timesteps <- 11          # visits/repeat observations
n_site <- 100               # locations/patients
n <- n_timesteps * n_site # number of observations
time <- rep(1:n_timesteps, each = n_site)
site <- rep(1:n_site, n_timesteps)

# observation-level corr matrix
R <- genPositiveDefMat(dim=m, #Number of columns/rows
                       covMethod = 'onion', 
                       rangeVar = c(1, 1), # Range of variances
                       eta=2)$Sigma

# observation-level ranefs
e <- matrix(nrow = n, ncol = m)
for (i in 1:n){
  e[i, ] <- rnorm(m) %*% t(chol(R))
}

# fixed effects design matrix, parameters, and linear predictor
# (note currently this is an intercept-only model, but X could be expanded)
k <- 1
X <- matrix(1, nrow = n, ncol = k)
beta <- matrix(rnorm(m * k), nrow = k)

# begin by simulating the state at time = 1
mu <- X %*% beta
z <- mu + qlogis(pnorm(e))
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
