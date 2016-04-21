# fitting the simpler model
library(rstan)
library(ggplot2)
#library(reshape)
source('sim_data_meta.R')

stan_d <- list(n_patients = n_site, 
               n_obs = n, 
               n_strains = m,
               patient = site,
               visit = time,
               z = y)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

m_init <- NULL
m_init <- stan('longitudinal_model.stan', data=stan_d, iter=10, chains=1)
m_fit <- stan(fit = m_init, data=stan_d, warmup=400, iter=1400, chains=2, #cores=3, 
              pars=c('cor_patient', 'cor_time', 'psi'))

traceplot(m_fit, pars=c('lp__'))

post <- extract(m_fit)

# evaluate recovery of R_p
traceplot(m_fit, pars = 'cor_patient', inc_warmup = TRUE)
Rp_df <- data.frame(Rp = c(Rp), y = ((m*2)^2):1)
plot(m_fit, pars = 'cor_patient') + 
  geom_point(data = Rp_df, aes(x = Rp, y = y), 
             size = 4, col = 'blue') + 
  ggtitle('Recovery of patient-level correlations')

# evaluate recovery of R_o
traceplot(m_fit, pars = 'cor_time', inc_warmup = TRUE)
Ro_df <- data.frame(Ro = c(Ro), y = ((m*2)^2):1)
plot(m_fit, pars = 'cor_time') + 
  geom_point(data = Ro_df, aes(x = Ro, y = y), 
             size = 4, col = 'blue') + 
  ggtitle('Recovery of observation-level correlations')

# evaluate recovery of z_tot (psi)
mean_psi <- apply(post$psi, c(2,3), mean)

plot(mean_psi[,2] ~ z_tot[,2], xlab="data", ylab="model")
abline(b=1, a=0)
