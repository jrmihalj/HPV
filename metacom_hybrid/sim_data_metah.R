# This is a hybrid model of Dorazio et al 2010 and Sebastian-Gonzalez et al. 2010
# Longitudinal effects are divided into effects on colonization and persistence
# There are fixed effects of each species at t-1 on each other species' col + persist
# There is also a correlated and nested random effect at time t.


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
n_timesteps <- 20          # visits/repeat observations
n_site <- 200               # locations/patients
n <- n_timesteps * n_site # number of observations
time <- rep(1:n_timesteps, each = n_site)
site <- rep(1:n_site, n_timesteps)

# Patient-level corr matrix
Rp <- genPositiveDefMat(dim=m, #Number of columns/rows
                        covMethod = 'onion', 
                        rangeVar = c(1, 1), # Range of variances
                        eta=2)$Sigma

#Rp <- matrix(c(1,0,0,1), ncol=m, byrow=T)

# observation-level ranefs
ep <- matrix(nrow = n_site, ncol = m)
for (i in 1:n_site){
  ep[i, ] <- rnorm(m) %*% t(chol(Rp))
}

# observation-level corr matrix
# Ro <- genPositiveDefMat(dim=m, #Number of columns/rows
#                         covMethod = 'onion', 
#                         rangeVar = c(1, 1), # Range of variances
#                         eta=2)$Sigma

Ro <- matrix(c(1,0,0,1), ncol=m, byrow=T)

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

#######
# PHI #
#######
# fixed effects design matrix, parameters, and linear predictor
# (note currently this is an intercept-only model, but X could be expanded)
k_phi <- 1
X_phi <- matrix(1, nrow = n, ncol = k_phi)
#beta_phi <- matrix(rnorm(m * k_phi), nrow = k_phi)
beta_phi <- matrix(rep(0,m*k_phi), nrow=k_phi)

# species-specific fixed effects on phi (used below)
#beta_phi_intxn <- matrix(rdoublex(m ^ 2, 0, 1), nrow = m)
beta_phi_intxn <- matrix(rep(0,m^2), ncol=m)
mu_phi <- X_phi %*% beta_phi

#########
# GAMMA #
#########
# fixed effects design matrix, parameters, and linear predictor
# (note currently this is an intercept-only model, but X could be expanded)
k_gam <- 1
X_gam <- matrix(1, nrow = n, ncol = k_gam)
#beta_gam <- matrix(rnorm(m * k_gam), nrow = k_gam)
beta_gam <- matrix(rep(0,m*k_phi), nrow=k_phi)


# species-specific fixed effects on gam (used below)
# Doesn't make sense to have a intra-specific effect on colonization 
# (because the species wasn't there in t-1), so these vals = 0
#beta_gam_intxn <- matrix(rdoublex(m ^ 2, 0, 1), nrow = m)
beta_gam_intxn <- matrix(rep(0,m^2), ncol=m)
diag(beta_gam_intxn) <- 0

mu_gam <- X_gam %*% beta_gam

##############
# Occurrence #
##############

z_tot <- matrix(NA, nrow=n, ncol=m)
y <- matrix(NA, nrow=n, ncol=m)

# first observation:
z_tot[time==1, ] =  qlogis(pnorm(e_all[time==1, ])) #only based on correlations
y[time == 1, ] <- ifelse(z_tot[time == 1, ] > 0, 1, 0)

# subsequent observations:
for (i in 2:n_timesteps) {
                        #Persistence:
  z_tot[time == i, ] <- (y[time == i - 1, ]) * (mu_phi[time == i - 1, ] + y[time == i - 1, ] %*% beta_phi_intxn) +
                        #Colonization:
                        (1 - y[time == i - 1, ]) * (mu_gam[time == i - 1, ] + y[time == i - 1, ] %*% beta_gam_intxn) +
                        #Correlated and nested random effects:
                        qlogis(pnorm(e_all[time==i, ]))
  
  y[time == i, ] <- ifelse(z_tot[time == i, ] > 0, 1, 0)
}
