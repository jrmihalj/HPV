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

sim_func <- function(n.pat = 20,
                     n.vis = 10, 
                     n.strains = 5,
                     init.psi = 0.5, #Initial occur. prob (0-1)
                     lap.lambda = 0.75, #lambda of the double exponential (laplace) distr.
                     eta = .4 #For positive definite matrix
  ){
  
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
  alpha.ij <- NULL #Indexed alpha.ij[Patient,Strain]
  alpha.ij <- rmnorm(n.pat, mean = rep(0, n.strains), varcov=Sig.across) # patient-level
  
  # Storage for looping
  cov.eff <- array(dim=c(n.pat,n.vis,n.strains)) #Summed covariate effects
  lpsi <- array(dim=c(n.pat,n.vis,n.strains)) #psi on logit scale
  psi <- array(dim=c(n.pat,n.vis,n.strains)) #psi in probability space
  Y <- array(dim=c(n.pat,n.vis,n.strains)) #Y[Patient, Visit, Strain]
  
  for(i in 1:n.pat){ # Patient
    
    for(j in 1:n.vis){ # Visit/Observation
      
      for(k in 1:n.strains){ # Strain/Species
        
        if(j == 1){ # For first observation, no cov. effects
          cov.eff[i,j,k] <- 0
        }else{
          cov.eff[i,j,k] <- sum(betas[k,] * X[i,j-1,])
        }
        
        lpsi[i,j,k] <- Logit(init.psi) + cov.eff[i,j,k] + alpha.ij[i,k]
        psi[i,j,k] <- AntiLogit(lpsi[i,j,k])
        
        Y[i,j,k] <- rbinom(1,1,psi[i,j,k])
        if(j<n.vis) X[i,j,k] <- Y[i,j,k]
        
      } # End k
    } # End j
  } # End i
  
  Patient <- melt(Y)[,1]
  Visit <- melt(Y)[,2]
  Strain <- melt(Y)[,3]
  Y_m <- melt(Y)$value
  psi_m <- melt(psi)$value
  
  # Create a data.frame
  Occ <- data.frame(Patient, Visit, Strain, psi_m, Y_m)
  pars <- list(alpha.ij = alpha.ij,
               X = X,
               Y = Y,
               psi = psi,
               betas = betas,
               Sig.across=Sig.across,
               Rho.across=Rho.across,
               n.pat = n.pat,
               n.vis = n.vis, 
               n.strains = n.strains,
               init.psi = init.psi,
               lap.lambda = lap.lambda,
               eta = eta)
  return(list(Occ=Occ, pars=pars))
}

##################
# Visualize
test <- sim_func()

test$pars$Rho.across
test$pars$betas

#Look at among-patient correlations (across all time samples)
cor(as.vector(test$pars$psi[,,2]),as.vector(test$pars$psi[,,3]))
plot(test$pars$psi[,,2]~test$pars$psi[,,3])

#Look at time effects (e.g. the effect of presence of strain 4 in t-1 on strain 5 in t)
# Choose 2 strains (look at effect of s1 in t-1 on s2 in t):
s1 <- 1 #t-1
s2 <- 3 #t 
temp.Y <- mat.or.vec(nr=test$pars$n.pat*(test$pars$n.vis-1), 3) # temp.Y[t-1, t, patient]
temp.Y[,1] <- as.vector(test$pars$Y[, 1:(test$pars$n.vis-1), s1]) #Choose a strain
temp.Y[,2] <- as.vector(test$pars$psi[, 2:test$pars$n.vis, s2]) #Choose a strain
temp.Y[,3] <- rep(c(1:test$pars$n.pat), times=(test$pars$n.vis-1)) #Choose a strain

plot(temp.Y[,2]~temp.Y[,1], xlab=paste("Strain",s1,"t-1"), ylab=paste("Strain",s2, "psi"))




