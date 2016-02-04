# Function for simulations:

require(mnormt) 
require(reshape2)
require(clusterGeneration)
require(smoothmest)

# Required  custom functions:
Logit <- function(x){
  log(x) - log(1 - x)
}

AntiLogit <- function(x){
  exp(x) / (exp(x) + 1)
}

sim_func <- function(n.pat = 5,
                     n.vis = 10, 
                     init.psi = 0.5, #Initial occur. prob (0-1)
                     lap.lambda = 0.75, #lambda of the double exponential (laplace) distr.
                     eta = .4, #For positive definite matrix
                     n.strains = 2
  ){
  
  n.strains <- n.strains #number of HPV strains

  # Fixed t-1 Species Co-occurrence Covariates (Storage)
  X <- array(dim=c(n.pat,n.vis-1,n.strains))
  
  # Draw the beta values from a laplace distribution:
  betas <- rdoublex(n.strains^2, 0, lap.lambda)
  betas <- matrix(betas, nrow=n.strains)
  
  #Covariance matrix (across patients) (co-occurrences at t, not t-1)
  Sig.across <- genPositiveDefMat(dim=n.strains, #Number of columns/rows
                                  covMethod = 'onion', 
                                  rangeVar = c(.01, .5), # Range of variances
                                  eta=eta)$Sigma
  
  #Derive correlation matrix
  Rho.across <- cov2cor(Sig.across)

  # Draw the random effects:
  eij <- NULL #Indexed eij[Patient,Strain]
  eij <- rmnorm(n.pat, mean = rep(0, n.strains), varcov=Sig.across) # patient-level
  
  cov.eff <- array(dim=c(n.pat,n.vis,n.strains))
  lpsi <- array(dim=c(n.pat,n.vis,n.strains))
  psi <- array(dim=c(n.pat,n.vis,n.strains))
  Y <- array(dim=c(n.pat,n.vis,n.strains)) # Y[Patient, Strain, Observation]
  
  for(i in 1:n.pat){
    
    for(j in 1:n.vis){
      
      for(k in 1:n.strains){
        
        if(j == 1){ # For first observation, no cov. effects
          cov.eff[i,j,k] <- 0
        }else{
          cov.eff[i,j,k] <- sum(betas[k,] * X[i,j-1,])
        }
        
        lpsi[i,j,k] <- Logit(init.psi) + cov.eff[i,j,k] + eij[i,k]
        psi[i,j,k] <- AntiLogit(lpsi[i,j,k])
        
        Y[i,j,k] <- rbinom(1,1,psi[i,j,k])
        if(j<n.vis) X[i,j,k] <- Y[i,j,k]
        
      } # End k
    } # End j
  } # End i
  
  Patient <- melt(Y)[,1]
  Visit <- melt(Y)[,2]
  Strain <- melt(Y)[,3]
  Y <- melt(Y)$value
  psi <- melt(psi)$value
  
  # Create a data.frame
  Occ <- data.frame(Patient, Visit, Strain, psi, Y)
  pars <- list(eij = eij,
               X = X,
               betas = betas,
               Sig.across=Sig.across, #Sig.within=Sig.within, 
               Rho.across=Rho.across)#, Rho_within=Rho_within)
  return(list(Occ=Occ, pars=pars))
}