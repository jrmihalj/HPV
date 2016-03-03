# fitting the simpler model
library(rstan)
library(ggplot2)
#library(reshape)
source('sim_data_metah.R')

stan_d <- list(n_unit = n_site, 
               n_time = n_timesteps,
               n = n, 
               m = m,
               unit = site,
               time = time,
               y = y,
               y_mat = y)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

m_init <- NULL
m_init <- stan('metah.stan', data=stan_d, iter=10, chains=1)

m_fit <- stan(fit = m_init, data=stan_d, warmup=100, iter=1100, chains=2, #cores=3, 
              pars=c('R_p', 'R_o', 'z', "beta_phi", "beta_gam"))

traceplot(m_fit, pars=c('lp__'))

post <- extract(m_fit)

# evaluate recovery of beta_phi
traceplot(m_fit, pars = 'beta_phi', inc_warmup = TRUE)
beta_phi_df <- data.frame(beta = c(t(beta_phi_intxn)), y = length(beta_phi_intxn):1)
plot(m_fit, pars = 'beta_phi') + 
  geom_point(data = beta_phi_df, aes(x = beta, y = y), 
             size = 4, col = 'blue') + 
  ggtitle('Recovery of species interaction t - 1 to t effects')

# evaluate recovery of beta_gam
traceplot(m_fit, pars = 'beta_gam', inc_warmup = TRUE)
beta_gam_df <- data.frame(beta = c(t(beta_gam_intxn)), y = length(beta_gam_intxn):1)
plot(m_fit, pars = 'beta_gam') + 
  geom_point(data = beta_gam_df, aes(x = beta, y = y), 
             size = 4, col = 'blue') + 
  ggtitle('Recovery of species interaction t - 1 to t effects')


# evaluate recovery of R_p
traceplot(m_fit, pars = 'R_p', inc_warmup = TRUE)
Rp_df <- data.frame(Rp = c(Rp), y = (m^2):1)
plot(m_fit, pars = 'R_p') + 
  geom_point(data = Rp_df, aes(x = Rp, y = y), 
             size = 4, col = 'blue') + 
  ggtitle('Recovery of patient-level correlations')

# evaluate recovery of R_o
traceplot(m_fit, pars = 'R_o', inc_warmup = TRUE)
Ro_df <- data.frame(Ro = c(Ro), y = (m^2):1)
plot(m_fit, pars = 'R_o') + 
  geom_point(data = Ro_df, aes(x = Ro, y = y), 
             size = 4, col = 'blue') + 
  ggtitle('Recovery of observation-level correlations')


# evaluate recovery of z_tot
mean_z <- apply(post$z, c(2,3), mean)

plot(mean_z[,1] ~ z_tot[,1], xlab="data", ylab="model")
abline(b=1, a=0)
