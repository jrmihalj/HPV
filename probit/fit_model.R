source('probit/sim_data.R')
library(rstan)
rstan_options(auto_write = TRUE)
stan_d <- list(n = n, 
               m = m, 
               k = k, 
               X = X, 
               y = y)

watch <- c('beta', 'R')

m_fit <- stan('probit/probit.stan', data = stan_d, cores = 2, chains = 2, 
              iter = 800, pars = watch, init_r = .1)

# evaluate recovery of beta
beta_df <- data.frame(beta = c(beta), y = length(beta):1)
plot(m_fit, pars = 'beta') + 
  geom_point(data = beta_df, aes(x = beta, y = y), size = 4, col = 'blue')

# evaluate recovery of R
R_df <- data.frame(R = c(R), y = (m^2):1)
plot(m_fit, pars = 'R') + 
  geom_point(data = R_df, aes(x = R, y = y), size = 4, col = 'blue')

traceplot(m_fit, 'beta')
traceplot(m_fit, 'R')

