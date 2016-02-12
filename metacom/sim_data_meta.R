# generative model, modified from Dorazio et al 2010

library(clusterGeneration)
library(smoothmest)

# Required  custom functions:
Logit <- function(x){
  log(x) - log(1 - x)
}

AntiLogit <- function(x){
  exp(x) / (exp(x) + 1)
}

# define parameters -------------------------------------
m <- 2                   # species
n_timesteps <- 11          # visits/repeat observations
n_site <- 150               # locations/patients
n <- n_timesteps * n_site # number of observations
time <- rep(1:n_timesteps, each = n_site)
site <- rep(1:n_site, n_timesteps)

# Patient-level corr matrix
Rp <- genPositiveDefMat(dim=m * 2, #Number of columns/rows (phi and gamma for each species)
                        covMethod = 'onion', 
                        rangeVar = c(1, 1), # Range of variances
                        eta=2)$Sigma

# observation-level ranefs
ep <- matrix(nrow = n_site, ncol = m * 2)
for (i in 1:n_site){
  ep[i, ] <- rnorm(m * 2) %*% t(chol(Rp))
}

# observation-level corr matrix
Ro <- genPositiveDefMat(dim=m *2 , #Number of columns/rows
                        covMethod = 'onion', 
                        rangeVar = c(1, 1), # Range of variances
                        eta=2)$Sigma

# observation-level ranefs
eo <- matrix(nrow = n, ncol = m * 2)
for (i in 1:n){
  eo[i, ] <- rnorm(m * 2) %*% t(chol(Ro))
}

# Add the ranef
e_all <- matrix(nrow = n, ncol = m * 2)
for(i in 1:n){
  e_all[i, ] <- ep[site[i], ] + eo[i, ]
}

# Normalize
e_all <- (e_all - mean(e_all)) / sd(e_all)

#######
# PHI #
#######
# fixed effects design matrix, parameters, and linear predictor
# (note currently this is an intercept-only model, but X could be expanded)
k_phi <- 1
X_phi <- matrix(1, nrow = n, ncol = k_phi)
beta_phi <- matrix(rnorm(m * k_phi), nrow = k_phi)

mu_phi <- X_phi %*% beta_phi
z_phi <- mu_phi + qlogis(pnorm(e_all[,1:m]))

#########
# GAMMA #
#########
# fixed effects design matrix, parameters, and linear predictor
# (note currently this is an intercept-only model, but X could be expanded)
k_gam <- 1
X_gam <- matrix(1, nrow = n, ncol = k_gam)
beta_gam <- matrix(rnorm(m * k_gam), nrow = k_gam)

mu_gam <- X_gam %*% beta_gam
z_gam <- mu_gam + qlogis(pnorm(e_all[,(m+1):ncol(e_all)]))

#######
# PSI #
#######

z_tot <- matrix(NA, nrow=n, ncol=m)
y <- matrix(NA, nrow=n, ncol=m)

z_tot[time==1, ] = rnorm(n_site * m) # .5 prob of initial occurrence...
y[time == 1, ] <- ifelse(z_tot[time == 1, ] > 0, 1, 0)

for (i in 2:n_timesteps) {
  z_tot[time == i, ] <- y[time == i - 1, ] * z_phi[time == i - 1, ] + (1 - y[time == i - 1, ]) * z_gam[time == i - 1, ]
  y[time == i, ] <- ifelse(z_tot[time == i, ] > 0, 1, 0)
}
