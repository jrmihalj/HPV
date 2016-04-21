#################################################
# Run the inference model on simulated data 
# Sylvia Ranjeva
# 2.2016
#################################################
source("sim_func.R")
require(mnormt) #needed for multivariate normal distrib
require(tidyr)
require(dplyr)
library(rstan)
library(scales)

Sim <- sim_func()
Occ <- Sim$Occ
df <- data.frame(pat = Occ$Patient, visit = Occ$Visit, strain = Occ$Strain, Y = Occ$Y_m)

l_df <- split(df, df$visit )
l_df <- lapply(l_df, function(x) { x["visit"] <- NULL; x })
l_df <- lapply(l_df, reshape, idvar = "pat", timevar = "strain", direction = "wide")
l_df <- lapply(l_df, function(x) { x["pat"] <- NULL; x })

X <- array(NA, dim=c(Sim$pars$n.vis, Sim$pars$n.pat, Sim$pars$n.strains))

for(i in 1:Sim$pars$n.vis){
  X[i,,] <- as.matrix(l_df[[i]])
}
 
#X <- aperm(X, c(3,1,2))
### Generate variables for stan model input:
modelInput <- list(
  n_strains = length(unique(Occ$Strain)),
  n_obs = nrow(Occ),
  n_patients = length(unique(Occ$Patient)),
  n_visits = max(Occ$Visit),
  Y = Occ$Y_m,
  strain = Occ$Strain,
  patient = Occ$Patient,
  visit_pat = Sim$Occ$Visit,
  X = X
)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

nChains <- 2
burn <- 100
n.thin <- 2
n.samp <- 1000
iter <- burn + (n.thin*n.samp)
modelFile <- "longitudinal_model.stan"

m_init <- NULL
m_init <- stan(modelFile, data=modelInput, iter=10, chains=1, pars='lp__')

fit <- NULL
fit <- stan(
  data = modelInput,
  file = modelFile,
  fit = m_init,
  warmup = burn,
  thin = n.thin,
  iter = iter,
  chains = nChains)

stan_trace(fit, 'lp__')
stan_trace(fit, "betas")

post <- extract(fit)

quartz()
stan_plot(fit, "cor_patient")

quartz()
par(mfrow=c(modelInput$n_strains, modelInput$n_strains), bty="n")
for(i in 1:modelInput$n_strains){
  for(j in 1:modelInput$n_strains){
    
    hist(post$betas[,i,j], main="", xlab=paste("Strains ", i, ",", j, sep=""))
    abline(v=Sim$pars$betas[i,j])
    
  }
}



summary <- summary(fit)$summary
library(ggmcmc)
ggd <- ggs(fit)
ggs_caterpillar(ggd, 'alpha_patient')
ggs_caterpillar(ggd, 'alpha_time')
