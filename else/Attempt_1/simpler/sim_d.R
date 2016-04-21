# Function for simulations:
require(tidyr)
require(dplyr)
require(clusterGeneration)

sim_d <- function(n.pat = 100,
                     n.vis = 10, 
                     ## eta controls the degree of correlation:
                     etaA = 100, # Among-patient eta
                     etaW = 100, # Within-patient eta
                     nstrain = 2
){
  
  # generate ranef cov matrices
  Sig.across <- genPositiveDefMat(dim=nstrain*2, #Number of columns/rows
                                  covMethod = 'onion', 
                                  rangeVar = c(.5, 3), # Range of variances
                                  eta=etaA)$Sigma
  
  # generate correlation matrices
  Rho_across <- cov2cor(Sig.across)

  # draw the random effects:
  alpha_patient <- matrix(rnorm(n.pat * nstrain * 2), 
                          nrow=n.pat) %*% t(chol(Sig.across))

  # patient-visit indices
  visit <- rep(1:n.vis, times=n.pat)
  
  # Calculate phis and gammas:
  lphi <- NULL
  lgam <- NULL
  
  mu_phi <- 0
  mu_gamma <- 0
  
  phi <- array(mu_phi, dim=c(n.pat, n.vis, nstrain))
  gamma <- array(mu_gamma, dim=c(n.pat, n.vis, nstrain))
  for (i in 1:n.pat) { 
    for (j in 1:n.vis){
      for (k in 1:nstrain){
        phi[i, j, k] <- alpha_patient[i, k]
        gamma[i, j, k] <- alpha_patient[i, k + nstrain]
      }
    }
  }
  phi <- plogis(phi)
  gamma <- plogis(gamma)
  
  # initial occupancy probability
  mu_psi0 <- 0
  sd_psi0 <- 1
  psi0 <- rnorm(n.pat * nstrain, mu_psi0, sd_psi0)
  
  # occurrence states
  psi <- array(dim = c(n.pat, n.vis, nstrain))
  z <- array(dim = c(n.pat, n.vis, nstrain))
  
  # t=1 
  psi[, 1, ] <- plogis(psi0)
  z[, 1, ] <- rbinom(n.pat * nstrain, 1, psi[, 1, ])
  
  # subsequent timesteps
  for (i in 1:n.pat) {
    for (t in 2:n.vis) {
      for (k in 1:nstrain) {
        psi[i, t, k] <- z[i, t - 1, k] * phi[i, t-1, k] + 
          (1 - z[i, t -1, k]) * gamma[i, t - 1, k]
        z[i, t, k] <- rbinom(1, 1, psi[i, t, k])
      }
    }
  }
  
  require(reshape2)
  mz <- reshape2::melt(z, varnames=c('patient', 'visit', 'strain'))
  names(mz)[names(mz) == 'value'] <- 'z'
  mpsi <- melt(psi, varnames=c('patient', 'visit', 'strain'), value.name='psi')
  names(mpsi)[names(mpsi) == 'value'] <- 'psi'
  occ <- full_join(mz, mpsi)

  # Create a data.frame
  pars <- list(Sig.across=Sig.across, #Sig.within=Sig.within, 
               Rho_across=Rho_across, #Rho_within=Rho_within, 
               alpha_patient = alpha_patient)
  return(list(occ=occ, pars=pars, z=z))
}