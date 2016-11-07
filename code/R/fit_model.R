# Fit the model to simulated data ------------------------------------------
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source('code/R/sim.R')
stan_d <- list(n = n, 
               n_strains = n_strains, 
               y = y,
               tbv = tbv,
               eta = 2,
               patient = patient, n_patient = n_patients,
               n_visit_max = n_visit_max, visit = visit,
               dir_prior = c(1, 1))


inits_f <- function(){
  list(betas_phi = array(rnorm(n_strains*n_strains,0,1), dim=c(n_strains,n_strains)),
       betas_gam = array(rnorm(n_strains*n_strains,0,1), dim=c(n_strains,n_strains)),
       betas_tbv_phi = rnorm(n_strains, 0, 1),
       betas_tbv_gam = rnorm(n_strains, 0, 1),
       alphas = rnorm(n_strains, -2, 1),
       e_patient = array(rnorm(n_patients*n_strains,0,1), dim=c(n_patients, n_strains)),
       abs_ystar = array(abs(rnorm(n*n_strains,0,2)), dim=c(n,n_strains))
       )
}

params <- c('Rho_patient', 'Rho_visit',
            'sd_visit', 'sd_patient', 'var_mat',
            'betas_phi', 'betas_gam', 
            'betas_tbv_phi', 'betas_tbv_gam', 'alphas')

test <- stan('code/stan/twolevel.stan',
              data = stan_d, chains = 1, iter = 2,
              init = inits_f,
              pars = params)

start_time <- Sys.time()
m_fit <- stan(fit = test,
              data = stan_d,
              init = inits_f,
              chains = 3, iter = 2000, warmup = 1000, #thin=3,
              pars = params,
              control=list(max_treedepth=13))
end_time <- Sys.time()

time_taken <- start_time - end_time
print(time_taken)
# check convergence:
summary(m_fit)$summary[,"Rhat"]

#check traceplots:
traceplot(m_fit, pars = 'lp__')

traceplot(m_fit, pars = 'var_mat')

traceplot(m_fit, pars = 'Rho_patient') +
  ylim(-1, 1)
traceplot(m_fit, pars = 'Rho_visit') +
  ylim(-1, 1)

traceplot(m_fit, pars = 'betas_phi')
traceplot(m_fit, pars = 'betas_gam')
traceplot(m_fit, pars = 'betas_tbv_phi')
traceplot(m_fit, pars = 'betas_tbv_gam')
traceplot(m_fit, pars = 'alphas')

# extract posterior draws
post <- extract(m_fit)

# evaluate parameter recovery of among patient correlations
d = m_fit@par_dims$Rho_patient[1]
quartz(height=10, width=10)
par(mfrow = c(d, d), mar = c(4, 3, 2, 2))
br <- seq(-1, 1, .05)
for (i in 1:d) {
  for (j in 1:d) {
    if (i == j) {
      plot.new()
    } else {
      colr=NULL
      colr = ifelse(findInterval(0, quantile(post$Rho_patient[, i, j], probs = c(0.025, 0.975)))==1,
                    "gray", "blue")
      hist(post$Rho_patient[, i, j], xlim = c(-1, 1), breaks = br, col = colr,
           main=paste("Rho_patient[", i, ",", j, "]", sep=""), xlab="")
      abline(v = quantile(post$Rho_patient[, i, j], probs = c(0.025, 0.975)))
      #abline(v = Rho_patient[i, j], col = 2, lwd = 2)
    }
  }
}

# check recovery of visit-level correlations
quartz(height=10, width=10)
par(mfrow = c(d, d), mar = c(4, 3, 2, 2))
br <- seq(-1, 1, .05)
for (i in 1:d) {
  for (j in 1:d) {
    if (i == j) {
      plot.new()
    } else {
      colr=NULL
      colr = ifelse(findInterval(0, quantile(post$Rho_visit[, i, j], probs = c(0.025, 0.975)))==1,
             "gray", "blue")
      hist(post$Rho_visit[, i, j], xlim = c(-1, 1), breaks = br, col = colr,
           main=paste("Rho_visit[", i, ",", j, "]", sep=""), xlab="")
      abline(v = quantile(post$Rho_visit[, i, j], probs = c(0.025, 0.975)))
      #abline(v = Rho_visit[i, j], col = 2, lwd = 2)
    }
  }
}

# Check correlation between Rho_visit and Rho_pat
RP_med = apply(post$Rho_patient, c(2,3), median)
RV_med = apply(post$Rho_visit, c(2,3), median)

plot(as.vector(RP_med) ~ as.vector(RV_med))
# No correlation

# check recovery of variances
quartz(height=10, width=10)
par(mfrow = c(2, d))
rev_i <- 2:1
for (i in 1:2) {
  for (j in 1:d) {
    hist(post$var_mat[, j, i], breaks = 30, xlim = range(post$var_mat))
    #abline(v = variance[j, rev_i[i]], col = 2, lwd = 3)
  }
}

# check recovery of betas_phi:
quartz(height=10, width=10)
par(mfrow = c(d, d), mar = c(4, 3, 2, 2))
for (i in 1:d) {
  for (j in 1:d) {
    colr=NULL
    colr = ifelse(findInterval(0, quantile(post$betas_phi[, i, j], probs = c(0.025, 0.975)))==1,
                  "gray", "blue")
    hist(post$betas_phi[, i, j], main="", 
         #xlim = range(betas_phi) + c(-1,1), 
         xlim = c(-1,1), col = colr,
         xlab=paste("beta_phi [", i, "," , j, "]"))
    abline(v = quantile(post$betas_phi[, i, j], probs = c(0.025, 0.975)))
    #abline(v = betas_phi[i, j], col = 2, lwd = 2)
  }
}

# check recovery of betas_gam:
quartz(height=10, width=10)
par(mfrow = c(d, d), mar = c(4, 3, 2, 2))
for (i in 1:d) {
  for (j in 1:d) {
    colr=NULL
    colr = ifelse(findInterval(0, quantile(post$betas_gam[, i, j], probs = c(0.025, 0.975)))==1,
                  "gray", "blue")
    hist(post$betas_gam[, i, j], main="", 
         #xlim = range(betas_gam) + c(-1,1),
         xlim = c(-1,1), col = colr,
         xlab=paste("beta_gam [", i, "," , j, "]"))
    abline(v = quantile(post$betas_gam[, i, j], probs = c(0.025, 0.975)))
    #abline(v = betas_gam[i, j], col = 2, lwd = 2)
  }
}

# check recovery of betas_tbv_phi and betas_tbv_gam:
quartz(height=10, width=5)
par(mfrow=c(d, 2), mar = c(4, 3, 2, 2))
for(i in 1:d){
  colr=NULL
  colr = ifelse(findInterval(0, quantile(post$betas_tbv_phi[, i], probs = c(0.025, 0.975)))==1,
                "gray", "blue")
  hist(post$betas_tbv_phi[, i], main="", 
       #xlim = range(betas_tbv_phi) + c(-1,1),
       xlim = c(-0.5,0.5), col = colr,
       xlab=paste("beta_tbv_phi [", i, "]"))
  abline(v = quantile(post$betas_tbv_phi[, i], probs = c(0.025, 0.975)))
  colr=NULL
  colr = ifelse(findInterval(0, quantile(post$betas_tbv_gam[, i], probs = c(0.025, 0.975)))==1,
                "gray", "blue")
  hist(post$betas_tbv_gam[, i], main="", 
       #xlim = range(betas_tbv_gam) + c(-1,1),
       xlim = c(-0.2,0.2), col = colr,
       xlab=paste("beta_tbv_gam [", i, "]"))
  abline(v = quantile(post$betas_tbv_gam[, i], probs = c(0.025, 0.975)))
}


# check recovery of alphas:
quartz(height=10, width=5)
par(mfrow=c(d,2), mar = c(4, 3, 2, 2))
for (i in 1:d) {
  hist(post$alphas[, i], main="", 
       #xlim = range(alphas) + c(-1,1),
       xlab=paste("alphas [", i, "]"))
  hist(pnorm(post$alphas[, i]), main="", 
       #xlim = range(alphas) + c(-1,1),
       xlim=c(0,0.10),
       xlab=paste("alphas_prob [", i, "]"))
  #abline(v = alphas[i], col = 2, lwd = 2)
}


