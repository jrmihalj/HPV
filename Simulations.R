##############################
# Simulations for 2 strains
##############################

# HEADER #

source("sim_func.R")
source("occur_func.R")
source("null_func.R")
require(mnormt) #needed for multivariate normal distrib
require(tidyr)
require(dplyr)

# Required  custom functions:
Logit <- function(x){
  log(x) - log(1 - x)
}

AntiLogit <- function(x){
  exp(x) / (exp(x) + 1)
}

##############################
# Test the functions:
##############################
test <- sim_func()
head(test)

occur_func(test)
null_func(test, 100)

##########################################################################################
##########################################################################################

# Question 1: 
# How do correlations affect coinfection?

# Simulation 1a - Similar/different strains, within-host corr in phi and gamma
# Methods:
# - Assign strong within-host correlation
# - Assign no among-host correlations
# - Assume no effect of covariates


n.sims <- 50
out <- list()

for(i in 1:n.sims){
  
  Sim <- sim_func(n.pat = 100, rag = -.9, rap = 0, rwg = 0, rwp = 0, 
                  bpat1p = 0, bpat2p = 0)
  
  # Run null model analysis:
  out[[i]] <- null_func(Sim, 1000)
  
}

# Get the co.dif in a usable format...
coinf <- sapply(out, "[", 4)
coinf <- as_data_frame(coinf) # from dplyr
coinf <- as.vector(t(coinf))
hist(coinf, breaks=20)

# Observations:
# - Within-host correlations in phi or gamma very rarely affect the number of coinfections.
# - This makes sense because correlations in psi (occurrence probability), which emerge from strong correlations in phi or gamma,
#   do not change the AVERAGE occurrence probabilities.
# - What these correlations would change is the frequency of observing the SAME outcome (pres/abs) for both strains within visits, 
#   not just co-occurrence.

# - However, correlations in phi ACROSS PATIENTS can affect coinfections. 
#   - Positive correlation can increase coinfection rate, negative correlation tends to decrease rates.
#   - This is because the correlation at the patient-level changes the average psi for that patient (therefore for all of that
#     patient's visits). This allows for within-patient co-occurrence to increase on average, even though overall psi for each strain
#     does not change. 

# - Correlations in gamma ACROSS PATIENTS do not seem to affect coinfections...
#   - Need to think about why

##########################################################################################
##########################################################################################

# Question 2: 
# How do among-host correlations affect coinfection?

# Simulation 1a - Similar/different strains, corr in phi (persistence)
# Methods:
# - Assign no among-host correlations
# - Assign strong within-host correlation in phi (i.e. strong competition/facilitation)
# - Assume no effect of covariates


n.sims <- 1
out <- list()

for(i in 1:n.sims){
  
  Sim <- sim_func(n.pat = 100, rag = 0.9, rap = 0, rwg = 0.9, rwp = 0, 
                  bpat1p = 0, bpat2p = 0)
  
  # Run null model analysis:
  out[[i]] <- null_func(Sim, 1000)
  
}


