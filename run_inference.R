#################################################
# Run the inference model on simulated data 
# Sylvia Ranjeva
# 9.2015
#################################################
source("sim_func.R")
source("occur_func.R")
source("null_func.R")
require(mnormt) #needed for multivariate normal distrib
require(tidyr)
require(dplyr)
library(rstan)

# Required  custom functions:
Logit <- function(x){
  log(x) - log(1 - x)
}

AntiLogit <- function(x){
  exp(x) / (exp(x) + 1)
}

#### Set params and run simulation
n.pat = 100
rag = .9
rap = 0
rwg = .9
rwp = 0
bpat1p = 0
bpat2p = 0

Sim <- sim_func(n.pat = 100, rag = 0.9, rap = 0, rwg = 0.9, rwp = 0, 
                bpat1p = 0, bpat2p = 0)

### Generate variables for stan model input:
Visit = rep(1:(nrow(Sim)/ length(unique(Sim$Strain))), each=length(unique(Sim$Strain)))

modelInput <- list(
  n_strains = length(unique(Sim$Strain)),
  n_obs = nrow(Sim),
  n_patients = length(unique(Sim$Patient)),
  n_visits_total = length(unique(Visit)),
  Y = Sim$Y,
  strain = Sim$Strain,
  patient = Sim$Patient,
  visit_pat = Sim$Visit.pat,
  Visit = Visit,
  X_time = Sim$X.time.vals,
  X_pat = Sim$X.pat.vals
)

for(j in 1:length(modelInput)){
  modelInput[[j]][is.na(modelInput[[j]])] <- 0
  #modelInput[[j]] <- as.data.frame(modelInput[[j]])
  assign(names(modelInput)[j], modelInput[[j]])
}

input <- c(names(modelInput))
nChains =1
iter = 1000
modelFile <- "metacommunity_model.stan"

fit <- stan(
  file = modelFile, 
  data = input,
  save_dso = TRUE,
  iter = iter,
  chains = nChains)

summary <- summary(fit)$summary
  
