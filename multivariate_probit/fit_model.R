source('multivariate_probit/sim.R')

library(rstan)
stan_d <- list(n = n, d = d, y = y, eta = 2)
m_fit <- stan('multivariate_probit/model.stan',
              data = stan_d, cores = 2, chains = 2, iter = 600,
              pars = c('Rho', 'Rho_prior'))
m_fit
traceplot(m_fit)

# evaluate parameter recovery
plot(density(extract(m_fit)$Rho[, 1, 2]), xlim = c(-1, 1),
     main = 'Posterior density: correlation')
lines(density(extract(m_fit)$Rho_prior[, 1, 2]), lty = 2)
abline(v = Rho[1, 2], col = 'red')
abline(v = cor(y_star[, 1], y_star[, 2]), col = 'blue')
legend('topright',
       col = c('red', 'blue', 'black'),
       legend = c('Population', 'Sample', 'Prior'),
       lty = c(1, 1, 2),
       bty = 'n')
