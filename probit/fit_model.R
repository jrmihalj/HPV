source('probit/sim_data.R')
library(rstan)
rstan_options(auto_write = TRUE)
stan_d <- list(n = nrow(y_use), 
               m = m, 
               k = k, 
               X = X_use, 
               X_intxn = X_intxn,
               y = y_use)
str(stan_d)

watch <- c('beta', 'beta_intxn', 'R')

m_fit <- stan('probit/probit.stan', data = stan_d, cores = 2, chains = 2, 
              iter = 800, pars = watch, init_r = .1)

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

# evaluate recovery of R
traceplot(m_fit, pars = 'R', inc_warmup = TRUE)
R_df <- data.frame(R = c(R), y = (m^2):1)
plot(m_fit, pars = 'R') + 
  geom_point(data = R_df, aes(x = R, y = y), 
             size = 4, col = 'blue') + 
  ggtitle('Recovery of observation-level correlations')
