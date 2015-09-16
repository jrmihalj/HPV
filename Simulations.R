##############################
# Simulations for 2 strains
##############################

# HEADER #

source("sim_func.R")
require(mnormt) #needed for multivariate normal distrib

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
# - Assign strong positive correlation in within-host gam (i.e. strong facilitatory priority effects)
# - Assume no effect of across-time covariates


Sim1a <- sim_func(n.pat = 100, rwp = 0.8, 
                  bpat1p = 0, bpat2p = 0, 
                  bpat1g = 5, bpat2g = 5,
                  btime1p = 0, btime2p = 0,
                  btime1g = 0, btime2g = 0)
Sim1a_null <- sim_func(n.pat = 100, rwp = 0, 
                       bpat1p = 0, bpat2p = 0, 
                       bpat1g = 5, bpat2g = 5,
                       btime1p = 0, btime2p = 0,
                       btime1g = 0, btime2g = 0)

# Now run a regression to see if we can recover the effect of among-patient covariate
plot(Y+runif(length(Y), 0,0.35) ~ X.pat.vals, data = subset(Sim1a_null, Strain==1))
plot(Y+runif(length(Y), 0,0.35) ~ X.pat.vals, data = subset(Sim1a, Strain==1))

# Pres counts
length(which(subset(Sim1a_null, Strain==1)$Y == 1))
length(which(subset(Sim1a, Strain==1)$Y == 1))

# Patient covariate effects
library(lme4)
sim1anull.mod <- glmer(Y ~ X.pat.vals + (1|Patient), family= binomial, data = subset(Sim1a_null, Strain==1))
sim1anull.mod@beta[1] # Intercept
sim1anull.mod@beta[2] # Slope

sim1a.mod <- glmer(Y ~ X.pat.vals + (1|Patient), family= binomial, data = subset(Sim1a, Strain==1))
sim1a.mod@beta[1] # Intercept
sim1a.mod@beta[2] # Slope

##############################

# Simulation 1b - Different strains, corr in gam
# Methods:
# - Assign no correlation between strains among hosts in any probabilities
# - Assign strong positive correlation in within-host gam (i.e. strong facilitatory priority effects)
# - Assume no effect of across-time covariates


Sim1b <- sim_func(n.pat = 100, rwp = 0.8, 
                  bpat1p = 0, bpat2p = 0, 
                  bpat1g = 5, bpat2g = -5,
                  btime1p = 0, btime2p = 0,
                  btime1g = 0, btime2g = 0)
Sim1b_null <- sim_func(n.pat = 100, rwp = 0, 
                       bpat1p = 0, bpat2p = 0, 
                       bpat1g = 5, bpat2g = -5,
                       btime1p = 0, btime2p = 0,
                       btime1g = 0, btime2g = 0)

# Now run a regression to see if we can recover the effect of among-patient covariate
plot(Y+runif(length(Y), 0,0.35) ~ X.pat.vals, data = subset(Sim1b_null, Strain==1))
plot(Y+runif(length(Y), 0,0.35) ~ X.pat.vals, data = subset(Sim1b, Strain==1))

# Pres counts
length(which(subset(Sim1b_null, Strain==1)$Y == 1))
length(which(subset(Sim1b, Strain==1)$Y == 1))

# Patient covariate effects
library(lme4)
sim1bnull.mod <- glmer(Y ~ X.pat.vals + (1|Patient), family= binomial, data = subset(Sim1b_null, Strain==1))
sim1bnull.mod@beta[1] # Intercept
sim1bnull.mod@beta[2] # Slope

sim1b.mod <- glmer(Y ~ X.pat.vals + (1|Patient), family= binomial, data = subset(Sim1b, Strain==1))
sim1b.mod@beta[1] # Intercept
sim1b.mod@beta[2] # Slope


