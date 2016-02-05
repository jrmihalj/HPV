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
df <- data.frame(pat = Occ$Patient, strain = Occ$Strain, Y = Occ$Y_m)

### Generate variables for stan model input:
modelInput <- list(
  n_strains = length(unique(Occ$Strain)),
  n_obs = nrow(Occ),
  n_patients = length(unique(Occ$Patient)),
  n_visits_total = max(Occ$Visit),
  Y = Occ$Y_m,
  strain = Occ$Strain,
  patient = Occ$Patient,
  visit_pat = Sim$Occ$Visit,
)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

nChains <- 3
iter <- 1200
modelFile <- "metacommunity_model_across_pt_alt.stan"

fit <- stan(
  fit = m_init, 
  data = modelInput,
  iter = iter,
  chains = nChains, 
  cores = nChains)

traceplot(fit, 'lp__')
traceplot(fit)

post <- extract(fit)

quartz()
br <- seq(0, .1 + max(post$tau_patient), .1)
alph <- .4
par(mfrow=c(modelInput$n_strains * 2, modelInput$n_strains * 2), bty='n')
for (i in 1:(modelInput$n_strains*2)){
  for (j in 1:(modelInput$n_strains*2)){
    if (i == j){
      hist(post$tau_patient[, i], breaks=br, main='', yaxt='n', col='grey', 
           ylab='')
      abline(v=sqrt(Sim$pars$Sig.across[i, j]), col='red', lty=2, lwd=2)
    } else if (i < j) {
        plot(density(post$cor_patient[, i, j]), main='', yaxt='n', xlab='', 
             xlim=c(-1, 1), ylab='')
        abline(v=Sim$pars$Rho_across[i, j])
      } else if (i > j) {
        i_m <- apply(post$alpha_patient[, , i], 2, median)
        j_m <- apply(post$alpha_patient[, , j], 2, median)
        plot(i_m, j_m, col=alpha(1, alph), cex=.4)
      }
  }
}


quartz()
br <- seq(0, .1 + max(post$tau_time), .1)
alph <- .2
par(mfrow=c(modelInput$n_strains * 2, modelInput$n_strains * 2), bty='n')
for (i in 1:(modelInput$n_strains*2)){
  for (j in 1:(modelInput$n_strains*2)){
    if (i == j){
      hist(post$tau_time[, i], breaks=br, main='', yaxt='n', col='grey', 
           ylab='')
      abline(v=sqrt(Sim$pars$Sig.within[i, j]), col='red', lty=2, lwd=2)
    } else if (i < j) {
      plot(density(post$cor_time[, i, j]), main='', yaxt='n', xlab='', 
           xlim=c(-1, 1), ylab='')
      abline(v=Sim$pars$Rho_within[i, j])
    } else if (i > j) {
      i_m <- apply(post$alpha_time[, , i], 2, median)
      j_m <- apply(post$alpha_time[, , j], 2, median)
      plot(i_m, j_m, col=alpha(1, alph), cex=.4)
    }
  }
}

summary <- summary(fit)$summary
library(ggmcmc)
ggd <- ggs(fit)
ggs_caterpillar(ggd, 'alpha_patient')
ggs_caterpillar(ggd, 'alpha_time')
