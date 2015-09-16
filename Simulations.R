##############################
# Simulations for 2 strains
##############################

# HEADER #

source("sim_func.R")
require(mnormt) #needed for multivariate normal distrib
require(lme4) #needed for GLMMs

# Required  custom functions:
Logit <- function(x){
  log(x) - log(1 - x)
}

AntiLogit <- function(x){
  exp(x) / (exp(x) + 1)
}

##############################
# Test (null) Function:
##############################
null <- sim_func()
head(null)
#null$Y

##########################################################################################
##########################################################################################

# Question 1: 
# Can hidden within-host correlations bias the estimation of 
# among-host covariate effects for pathogen strains?

# Simulation 1a - Similar strains, corr in gam
# Methods:
# - Assign no correlation between strains among hosts in any probabilities
# - Assign strong correlation in within-host gam (i.e. strong priority effects)
# - Assume no effect of across-time covariates


nsims <- 100
ints <- mat.or.vec(nr=nsims,nc=2)
slops <- mat.or.vec(nr=nsims,nc=2)
pres <- mat.or.vec(nr=nsims,nc=2)

for(i in 1:nsims){
  Sim <- sim_func(n.pat = 100, rwg = 0.8, 
                  bpat1g = 3, bpat2g = 3)
  
  Sim_null <- sim_func(n.pat = 100, rwg = 0, 
                       bpat1g = 3, bpat2g = 3)
  
  Simnull.mod <- glmer(Y ~ X.pat.vals + (1|Patient), family= binomial, data = subset(Sim_null, Strain==1))
  Sim.mod <- glmer(Y ~ X.pat.vals + (1|Patient), family= binomial, data = subset(Sim, Strain==1))
  
  ints[i,1] <- Simnull.mod@beta[1]
  ints[i,2] <- Sim.mod@beta[1]
  
  slops[i,1] <- Simnull.mod@beta[2]
  slops[i,2] <- Sim.mod@beta[2]
  
  # Pres counts
  pres[i,1] <- length(which(subset(Sim_null, Strain==1)$Y == 1))
  pres[i,2] <- length(which(subset(Sim, Strain==1)$Y == 1))
  
}
# White = null, Gray = alternative
hist(ints[,1], breaks=20)
hist(ints[,2], breaks=20, add=T, col="gray")

hist(slops[,1], breaks=20)
hist(slops[,2], breaks=20, add=T, col="gray")

hist(pres[,1], breaks=20)
hist(pres[,2], breaks=20, add=T, col="gray")


# Observations:
# - No obvious effect on pres counts, intercepts or slopes. 

##############################

# Simulation 1b - Different strains, corr in gam
# Methods:
# - Assign no correlation between strains among hosts in any probabilities
# - Assign strong correlation in within-host gam (i.e. strong facilitatory priority effects)
# - Assume no effect of across-time covariates


nsims <- 100
ints <- mat.or.vec(nr=nsims,nc=2)
slops <- mat.or.vec(nr=nsims,nc=2)
pres <- mat.or.vec(nr=nsims,nc=2)

for(i in 1:nsims){
  Sim <- sim_func(n.pat = 100, rwg = 0.8, 
                  bpat1g = 0, bpat2g = -5)
  
  Sim_null <- sim_func(n.pat = 100, rwg = 0, 
                       bpat1g = 0, bpat2g = -5)
  
  Simnull.mod <- glmer(Y ~ X.pat.vals + (1|Patient), family= binomial, data = subset(Sim_null, Strain==1))
  Sim.mod <- glmer(Y ~ X.pat.vals + (1|Patient), family= binomial, data = subset(Sim, Strain==1))
  
  ints[i,1] <- Simnull.mod@beta[1]
  ints[i,2] <- Sim.mod@beta[1]
  
  slops[i,1] <- Simnull.mod@beta[2]
  slops[i,2] <- Sim.mod@beta[2]
  
  # Pres counts
  pres[i,1] <- length(which(subset(Sim_null, Strain==1)$Y == 1))
  pres[i,2] <- length(which(subset(Sim, Strain==1)$Y == 1))
  
}

# White = null, Gray = alternative
hist(ints[,2], breaks=20)
hist(ints[,1], breaks=20, add=T, col="gray")

hist(slops[,1], breaks=20)
hist(slops[,2], breaks=20, add=T, col="gray")

hist(pres[,1], breaks=20)
hist(pres[,2], breaks=20, add=T, col="gray")


# Observations:
# -


##############################

# Simulation 1c - Similar strains, corr in phi
# Methods:
# - Assign no correlation between strains among hosts in any probabilities
# - Assign strong positive correlation in within-host phi (i.e. strong facilitatory priority effects)
# - Assume no effect of across-time covariates

nsims <- 100
ints <- mat.or.vec(nr=nsims,nc=2)
slops <- mat.or.vec(nr=nsims,nc=2)
pres <- mat.or.vec(nr=nsims,nc=2)

for(i in 1:nsims){
  Sim <- sim_func(n.pat = 100, rwp = 0.8, 
                  bpat1p = 3, bpat2p = 3)
  
  Sim_null <- sim_func(n.pat = 100, rwp = 0, 
                       bpat1p = 3, bpat2p = 3)
  
  Simnull.mod <- glmer(Y ~ X.pat.vals + (1|Patient), family= binomial, data = subset(Sim_null, Strain==1))
  Sim.mod <- glmer(Y ~ X.pat.vals + (1|Patient), family= binomial, data = subset(Sim, Strain==1))
  
  ints[i,1] <- Simnull.mod@beta[1]
  ints[i,2] <- Sim.mod@beta[1]
  
  slops[i,1] <- Simnull.mod@beta[2]
  slops[i,2] <- Sim.mod@beta[2]
  
  # Pres counts
  pres[i,1] <- length(which(subset(Sim_null, Strain==1)$Y == 1))
  pres[i,2] <- length(which(subset(Sim, Strain==1)$Y == 1))
  
}
# White = null, Gray = alternative
hist(ints[,1], breaks=20)
hist(ints[,2], breaks=20, add=T, col="gray")

hist(slops[,1], breaks=20)
hist(slops[,2], breaks=20, add=T, col="gray")

hist(pres[,1], breaks=20)
hist(pres[,2], breaks=20, add=T, col="gray")

##########################################################################################
##########################################################################################

# Question 2: 
# Can among-host correlations bias the estimation of 
# among-host covariate effects for pathogen strains?

# Simulation 2a - Similar strains, corr in gam
# Methods:
# - Assign no correlation between strains among hosts in any probabilities
# - Assign strong correlation in among-host gam (i.e. strong priority effects)
# - Assume no effect of across-time covariates


nsims <- 100
ints <- mat.or.vec(nr=nsims,nc=2)
slops <- mat.or.vec(nr=nsims,nc=2)
pres <- mat.or.vec(nr=nsims,nc=2)

for(i in 1:nsims){
  Sim <- sim_func(n.pat = 100, rag = 0.8, 
                  bpat1g = 3, bpat2g = 3)
  
  Sim_null <- sim_func(n.pat = 100, rag = 0, 
                       bpat1g = 3, bpat2g = 3)
  
  Simnull.mod <- glmer(Y ~ X.pat.vals + (1|Patient), family= binomial, data = subset(Sim_null, Strain==1))
  Sim.mod <- glmer(Y ~ X.pat.vals + (1|Patient), family= binomial, data = subset(Sim, Strain==1))
  
  ints[i,1] <- Simnull.mod@beta[1]
  ints[i,2] <- Sim.mod@beta[1]
  
  slops[i,1] <- Simnull.mod@beta[2]
  slops[i,2] <- Sim.mod@beta[2]
  
  # Pres counts
  pres[i,1] <- length(which(subset(Sim_null, Strain==1)$Y == 1))
  pres[i,2] <- length(which(subset(Sim, Strain==1)$Y == 1))
  
}
# White = null, Gray = alternative
hist(ints[,1], breaks=20)
hist(ints[,2], breaks=20, add=T, col="gray")

hist(slops[,1], breaks=20)
hist(slops[,2], breaks=20, add=T, col="gray")

hist(pres[,1], breaks=20)
hist(pres[,2], breaks=20, add=T, col="gray")


# Observations:
# - No obvious effect on pres counts, intercepts or slopes. 

