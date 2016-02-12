source('probit/sim_data.R')
library(rstan)
rstan_options(auto_write = TRUE)
stan_d <- list(n = nrow(y_use), 
               m = m, 
               k = k, 
               n_unit = n_site,
               unit = site[1:nrow(y_use)],
               X = X_use, 
               X_intxn = X_intxn,
               y = y_use)
str(stan_d)

watch <- c('beta', 'beta_intxn', 'R_p', 'R_o')

m_fit <- NULL
m_fit <- stan('probit/probit.stan', data = stan_d, cores = 2, chains = 2, 
              iter = 1000, pars = watch, init_r = .1)

traceplot(m_fit, 'lp__')

# evaluate recovery of beta
traceplot(m_fit, pars = 'beta', inc_warmup = TRUE)
beta_df <- data.frame(beta = c(beta), y = length(beta):1)
plot(m_fit, pars = 'beta') + 
  geom_point(data = beta_df, aes(x = beta, y = y), 
             size = 4, col = 'blue') + 
  ggtitle('Recovery of fixed effects')

# evaluate recovery of beta_intxn
traceplot(m_fit, pars = 'beta_intxn', inc_warmup = TRUE)
beta_intxn_df <- data.frame(beta = c(t(beta_intxn)), y = length(beta_intxn):1)
plot(m_fit, pars = 'beta_intxn') + 
  geom_point(data = beta_intxn_df, aes(x = beta, y = y), 
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
plot(m_fit, pars = 'R_p') + 
  geom_point(data = Ro_df, aes(x = Ro, y = y), 
             size = 4, col = 'blue') + 
  ggtitle('Recovery of observation-level correlations')
