source('multivariate_probit/sim.R')

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


stan_d <- list(n = n, d = d, y = y, eta = 2,
               patient = patient, n_patient = n_patients,
               n_visit = n_visits, visit = visit,
               dir_prior = c(1, 1))

m_fit <- stan('multivariate_probit/twolevel.stan',
              data = stan_d, chains = 2, iter = 1000, #warmup = 1000,
              pars = c('Rho_patient', 'Rho_visit',
                       'sd_visit', 'sd_patient', 'var_mat',
                       'betas'))

# check convergence:
summary(m_fit)$summary[,"Rhat"]

#check traceplots:
traceplot(m_fit, pars = 'lp__')

traceplot(m_fit, pars = 'var_mat')

traceplot(m_fit, pars = 'Rho_patient') +
  ylim(-1, 1)
traceplot(m_fit, pars = 'Rho_visit') +
  ylim(-1, 1)

traceplot(m_fit, pars = 'betas')

# extract posterior draws
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

# check recovery of betas:
par(mfrow=c(d,d))
for (i in 1:d) {
  for (j in 1:d) {
    hist(post$betas[, i, j], main="", xlab=paste("beta [", i, "," , j, "]"))
    abline(v = betas[i, j], col = 2, lwd = 2)
  }
}

