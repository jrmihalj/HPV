# Function for simulations:

require(mnormt) 
require(tidyr)
require(dplyr)
require(clusterGeneration)

# Required  custom functions:
Logit <- function(x){
  log(x) - log(1 - x)
}

AntiLogit <- function(x){
  exp(x) / (exp(x) + 1)
}

sim_func <- function(n.pat = 100,
                     n.vis = 10, 
                     # Within and among host covariate effects for phi and gamma:
                     bpat1g = 0, 
                     bpat2g = 0,
                     btime1g = 0,
                     btime2g = 0,
                     bpat1p = 0,
                     bpat2p = 0,
                     btime1p = 0,
                     btime2p = 0,
                     # Correlations:
                     manual = 1,
                     # If 'manual' == T, then specify correlations
                     ## Among-patients:
                     raG1G2 = 0,
                     raP1P2 = 0, 
                     raP1G1 = 0,
                     raP2G2 = 0,
                     raP1G2 = 0,
                     raP2G1 = 0,
                     ## Within-patients:
                     rwG1G2 = 0,
                     rwP1P2 = 0, 
                     rwP1G1 = 0,
                     rwP2G2 = 0,
                     rwP1G2 = 0,
                     rwP2G1 = 0,
                     ## sd across patients for each strain (phi and gamma):
                     ## Assume all sd equal
                     sa = 1,
                     ## sd within patients for each strain (phi and gamma):
                     ## Assume all sd equal
                     sw = .4,
                     # If 'manual' == FALSE, generate pos. definite covariance matrix automatically
                     ## eta controls the degree of correlation:
                     etaA = 2, # Among-patient eta
                     etaW = 2, # Within-patient eta
                     #Global probabilities:
                     globphi = .5, 
                     globgam = .5, 
                     globpsi = .5
  ){
  
  n.strains <- 2 #number of HPV strains
  n.patients <- n.pat #number of patients tested
  n.visits <- rep(n.vis, n.patients) #number of visits per patient
  n.obs.perpat <- n.strains*n.visits #number of "tests" per patient
  n.obs <- sum(n.obs.perpat) #total number of "tests" (One test per strain per patient per visit)
  
  # Indices:
  Strain <- rep(c(1,2), times=n.obs/n.strains) #Strain index
  Patient <- NULL
  for(i in 1:n.patients){
    Patient <- c(Patient, rep(i, times=n.obs.perpat[i]))
  }
  # Overall visit:
  Visit <- rep(1:sum(n.visits), each=n.strains) 
  # Visit number per patient:
  Visit.Pat <- NULL
  for(i in 1:n.patients){
    Visit.Pat <- c(Visit.Pat, rep(1:n.visits[i], each=n.strains))
  }
  
  # Covariates (centered to mean and scaled to 1 sd):
  q.pat <- 1 #number of covariates measured at the patient-level (but not over time)
  q.time <- 1 #number of covariates measured for each patient at each visit
  
  X.pat.vals <- NULL
  for(i in 1:n.patients){
    X.pat.vals <- c(X.pat.vals, rep(rnorm(1,0,1), times=n.obs.perpat[i]))
  }
  
  X.time.vals <- NULL
  for(i in 1:n.patients){
    X.time.vals <- c(X.time.vals, rep(rnorm(n.visits[i],0,1), each=n.strains))
  }
  
  # Covariate responses:
  beta.pat.phi <- c(bpat1p, bpat2p)
  beta.time.phi <- c(btime1p, btime2p)
  beta.pat.gam <- c(bpat1g, bpat2g)
  beta.time.gam <- c(btime1g, btime2g)
  
  #Assign the global mean values
  mean.phi <- Logit(globphi)
  mean.gam <- Logit(globgam)
  mean.psi <- Logit(globpsi)
  
  #Correlations:
  if(manual == 1){
    
    ## Among-patients:
    
    #Covariance matrix:
    #Phi1, Phi2, Gamma1, Gamma2 = columns
    Sig.across <- matrix(
      c(sa^2, raP1P2*sa^2, raP1G1*sa^2, raP1G2*sa^2,
        raP1P2*sa^2, sa^2, raP2G1*sa^2, raP2G2*sa^2,
        raP1G1*sa^2, raP2G1*sa^2, sa^2, raG1G2*sa^2,
        raP1G2*sa^2, raP2G2*sa^2, raG1G2*sa^2, sa^2),
      ncol=4
    )
    
    ## Within-patients:
    
    #Covariance matrix:
    #Phi1, Phi2, Gamma1, Gamma2 = columns
    
    Sig.within <- matrix(
      c(sw^2, rwP1P2*sw^2, rwP1G1*sw^2, rwP1G2*sw^2,
        rwP1P2*sw^2, sw^2, rwP2G1*sw^2, rwP2G2*sw^2,
        rwP1G1*sw^2, rwP2G1*sw^2, sw^2, rwG1G2*sw^2,
        rwP1G2*sw^2, rwP2G2*sw^2, rwG1G2*sw^2, sw^2),
      ncol=4
    )
    
  }else{
    
    Sig.across <- genPositiveDefMat(dim=n.strains*2, #Number of columns/rows
                                    covMethod = 'onion', 
                                    rangeVar = c(.5, 1), # Range of variances
                                    eta=etaA)$Sigma
    
    Sig.within <- genPositiveDefMat(dim=n.strains*2, #Number of columns/rows
                                    covMethod = 'onion', 
                                    rangeVar = c(.05, .3), 
                                    eta=etaW)$Sigma
    
  }
  
  
  # Draw the random effects:
  # Phi1, Phi2, Gamma1, Gamma2 = columns
  b.0.k <- rmnorm(n.patients, mean = rep(0, n.strains*2), varcov=Sig.across) # patient-level
  b.0.t <- rmnorm(sum(n.visits), mean = rep(0, n.strains*2), varcov=Sig.within) # visit-level
  
  # Calculate phis and gammas:
  lphi <- NULL
  lgam <- NULL
  
  # Formula from our model:
  # (phi/gamma) <- global_mean + intercept_among + covariate_among + 
  #                     intercept_within + covariate_within
  for(j in 1:n.obs){
    lphi[j] <- mean.phi + b.0.k[Patient[j],Strain[j]] + beta.pat.phi[Strain[j]] * X.pat.vals[j] + 
      b.0.t[Visit[j],Strain[j]] + beta.time.phi[Strain[j]] * X.time.vals[j]
    lgam[j] <- mean.gam + b.0.k[Patient[j],(n.strains + Strain[j])] + beta.pat.gam[Strain[j]] * X.pat.vals[j] + 
      b.0.t[Visit[j],(n.strains + Strain[j])] + beta.time.gam[Strain[j]] * X.time.vals[j]
  }
  
  
  # First we need the initial occurrence probability for each strain for each patient.
  # We assume this probability is not influenced by covariates, but this assumption
  # could be relaxed in the future. Instead we assume the initial probabilities are
  # generated randomly from a global occurrence probability.
  
  # Calculate AntiLogit for phi and gam
  phi <- AntiLogit(lphi)
  gam <- AntiLogit(lgam)
  
  psi <- NULL
  for(j in 1:n.obs){
    
    if(Visit.Pat[j]==1){ # If it's the patient's first visit
      psi[j] <- AntiLogit(rnorm(1,mean.psi,.2))
    }else{
      psi[j] <- phi[j-n.strains] * psi[j-n.strains] + gam[j-n.strains] * (1 - psi[j-n.strains]) 
    }
    
  }
  
  Y <- NULL
  for(z in 1:n.obs){
    Y[z] <- rbinom(1, 1, psi[z])
  }
  
  
  
  # Create a data.frame
  Occ <- data.frame(Y, psi, phi, gam, Strain, Patient, Visit.Pat, X.time.vals, X.pat.vals)
  
  return(Occ)
}