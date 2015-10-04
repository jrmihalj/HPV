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
library(scales)

# Required  custom functions:
Logit <- function(x){
  log(x) - log(1 - x)
}

AntiLogit <- function(x){
  exp(x) / (exp(x) + 1)
}

#### Set params and run simulation
n.pat = 200
bpat1p = 0
bpat2p = 0

Sim <- sim_func(n.pat = n.pat, n.vis=10, manual=FALSE)

Sim$Occ %>%
  ggplot(aes(x=phi, y=gam, color=factor(Patient))) + 
  facet_wrap(~Strain) + 
  geom_point()

if (n.pat < 300){
  Sim$Occ %>%
    ggplot(aes(x=Visit.Pat, y=psi)) + 
    geom_line(aes(color=factor(Strain))) + 
    facet_wrap(~Patient)
}

### Generate variables for stan model input:
Visit = rep(1:(nrow(Sim$Occ)/ length(unique(Sim$Occ$Strain))), 
            each=length(unique(Sim$Occ$Strain)))

modelInput <- list(
  n_strains = length(unique(Sim$Occ$Strain)),
  n_obs = nrow(Sim$Occ),
  n_patients = length(unique(Sim$Occ$Patient)),
  n_visits_total = length(unique(Visit)),
  Y = Sim$Occ$Y,
  strain = Sim$Occ$Strain,
  patient = Sim$Occ$Patient,
  visit_pat = Sim$Occ$Visit.Pat,
  Visit = Visit,
  X_time = Sim$Occ$X.time.vals,
  X_pat = Sim$Occ$X.pat.vals
)

nChains <- 3
iter <- 1000
modelFile <- "metacommunity_model.stan"

watch <- c('phi_mean', 'gam_mean', 'psi_mean', 
           #'beta_pat_mean', 'beta_time_mean', 
           #'beta_pat_sd', 'beta_time_sd', 
           #'beta_pat', 'beta_time', 
           'tau_patient', 'alpha_patient', #'tau_time', 
           'cor_patient', 'psi'#, #'cor_time'
           )
if (!('m_init' %in% ls())){
  m_init <- stan(modelFile, data=modelInput, iter=10, chains=1, pars='lp__')
}
fit <- stan(
  fit = m_init, 
  data = modelInput,
  iter = iter,
  chains = nChains, 
  cores = nChains, 
  pars = watch)

traceplot(fit, 'lp__')
traceplot(fit)

traceplot(fit, 'cor_patient')
Sim$pars$Rho_across
traceplot(fit, 'tau_patient')
sqrt(diag(Sim$pars$Sig.across))

library(ggmcmc)
ggd <- ggs(fit)
ggs_caterpillar(ggd, 'alpha_patient')

post <- extract(fit)

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
             xlim=c(-1, 1))
        abline(v=Sim$pars$Rho_across[i, j])
      } else if (i > j) {
        i_m <- apply(post$alpha_patient[, , i], 2, median)
        j_m <- apply(post$alpha_patient[, , j], 2, median)
        plot(i_m, j_m, col=alpha(1, alph), cex=.4)
      }
  }
}


summary <- summary(fit)$summary
