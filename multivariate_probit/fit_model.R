source('multivariate_probit/sim.R')

library(rstan)
stan_d <- list(n = n, d = d, y = y, eta = 2,
               patient = patient, n_patient = n_patients,
               dir_prior = c(1, 1))
m_fit <- stan('multivariate_probit/twolevel.stan',
              data = stan_d, cores = 2, chains = 2, iter = 1000,
              pars = c('Rho_patient', 'Rho_visit',
                       'sd_visit', 'sd_patient', 'var_mat'))
m_fit
traceplot(m_fit, pars = c('var_mat'))

traceplot(m_fit, pars = 'Rho_patient', inc_warmup = TRUE) +
  ylim(-1, 1)
traceplot(m_fit, pars = 'Rho_visit', inc_warmup = TRUE) +
  ylim(-1, 1)

post <- extract(m_fit)

# evaluate parameter recovery of among patient correlations
par(mfrow = c(d, d), mar = c(4, 3, 2, 2))
br <- seq(-1, 1, .05)
for (i in 1:d) {
  for (j in 1:d) {
    if (i == j) {
      plot.new()
    } else {
      hist(post$Rho_patient[, i, j], xlim = c(-1, 1), breaks = br)
      abline(v = Rho_patient[i, j], col = 2, lwd = 2)
    }
  }
}

# check recovery of visit-level correlations
for (i in 1:d) {
  for (j in 1:d) {
    if (i == j) {
      plot.new()
    } else {
      hist(post$Rho_visit[, i, j], xlim = c(-1, 1), breaks = br)
      abline(v = Rho_visit[i, j], col = 2, lwd = 2)
    }
  }
}

# check recovery of variances
par(mfrow = c(2, d))
rev_i <- 2:1
for (i in 1:2) {
  for (j in 1:d) {
    hist(post$var_mat[, j, i], breaks = 30, xlim = range(post$var_mat))
    abline(v = variance[j, rev_i[i]], col = 2, lwd = 3)
  }
}

