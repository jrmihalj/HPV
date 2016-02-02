# fitting the simpler model
source('simpler/sim_d.R')
d <- sim_d(n.pat = 200,
           n.vis = 10, 
           ## eta controls the degree of correlation:
           etaA = .1, # Among-patient eta
           nstrain = 2)
library(ggplot2)
d$occ %>%
  ggplot(aes(x=visit, y=z, col=factor(strain))) + 
  geom_line() + 
  facet_wrap(~patient)

d$occ %>%
  ggplot(aes(x=visit, y=psi, col=factor(strain))) + 
  geom_line() + 
  facet_wrap(~patient)

library(rstan)
stan_d <- list(n_patients = dim(d$z)[1], 
               n_vis = dim(d$z)[2], 
               n_strains = dim(d$z)[3], 
               z=d$z)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

m_init <- stan('simpler/mod.stan', data=stan_d, iter=10, chains=1)
m_fit <- stan(fit = m_init, data=stan_d, iter=600, chains=3, #cores=3, 
              pars=c('cor_patient', 'tau_patient', 'alpha_patient'))

traceplot(m_fit, pars=c('cor_patient', 'tau_patient'))

post <- extract(m_fit)

library(scales)
br <- seq(0, .1 + max(post$tau_patient), .1)
alph <- .4
par(mfrow=c(stan_d$n_strains * 2, stan_d$n_strains * 2), bty='n', 
    mar=c(2, 2, 1, 1))
for (i in 1:(stan_d$n_strains*2)){
  for (j in 1:(stan_d$n_strains*2)){
    if (i == j){
      hist(post$tau_patient[, i], breaks=br, main='', yaxt='n', col='grey', 
           ylab='')
      abline(v=sqrt(d$pars$Sig.across[i, j]), col='red', lty=3, lwd=2)
      abline(v=sd(d$pars$alpha_patient[, j]), col='blue', lty=2, lwd=2)
    } else if (i < j) {
      plot(density(post$cor_patient[, i, j]), main='', yaxt='n', xlab='', 
           xlim=c(-1, 1), ylab='')
      abline(v=d$pars$Rho_across[i, j], col='red', lty=3, lwd=2)
      abline(v=cor(d$pars$alpha_patient[, i], d$pars$alpha_patient[, j]), 
             col='blue', lty=2, lwd=2)
    } else if (i > j) {
      i_m <- apply(post$alpha_patient[, , i], 2, median)
      j_m <- apply(post$alpha_patient[, , j], 2, median)
      plot(i_m, j_m, col=alpha(1, alph), cex=.4)
    }
  }
}
# sample quantities in dashed blue, population params in red

ap <- melt(post$alpha_patient, varnames = c('iter', 'patient', 'column')) %>%
  group_by(patient, column) %>%
  summarize(med = median(value))

true_ap <- melt(d$pars$alpha_patient, varnames=c('patient', 'column'))

apd <- full_join(ap, true_ap)
ggplot(apd, aes(x=value, y=med)) + 
  geom_point() + 
  geom_abline(yintercept=0, slope=1, linetype='dashed') + 
  xlab('True value') + 
  ylab('Posterior median') + 
  facet_wrap(~column)

