# Test how parameters affect observed correlations in occupancy probability

n.pat = 200;
n.vis = 10; 
# Within and among host covariate effects for phi and gamma:
bpat1g = 0; 
bpat2g = 0;
btime1g = 0;
btime2g = 0;
bpat1p = 0;
bpat2p = 0;
btime1p = 0;
btime2p = 0;
rag = -0.9; #rho.across.gamma
rap = -0.9; #rho.across.phi
rwg = 0.9; #rho.within.gamma
rwp = 0.9; #rho.within.phi
# sd across patients for each strain (phi and gamma):
sap1 = .6; 
sap2 = .6;
sag1 = .6;
sag2 = .6;
# sd within patients for each strain (phi and gamma):
swp1 = .2; 
swp2 = .2;
swg1 = .2;
swg2 = .2;
#Global probabilities:
globphi = .5; 
globgam = .5; 
globpsi = .5 


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

#Across-patient std.dev. in persistence for each strain:
sig.across.phi.1 <- sap1
sig.across.phi.2 <- sap2

#Across-patient std.dev. in colonization for each strain:
sig.across.gam.1 <- sag1
sig.across.gam.2 <- sag2

#Within-patient std.dev. in persistence for each strain:
sig.within.phi.1 <- swp1
sig.within.phi.2 <- swp2

#within-patient std.dev. in colonization for each strain:
sig.within.gam.1 <- swg1
sig.within.gam.2 <- swg2

#assign the global mean values
mean.phi <- Logit(globphi)
mean.gam <- Logit(globgam)
mean.psi <- Logit(globpsi)

i <- 1

#Correlation:
if(i == 1) {rho.across.phi <- rap} else {rho.across.phi <- 0}

#Covariance matrix:
Sig.across.phi <- matrix(
  c(sig.across.phi.1^2, rho.across.phi*sig.across.phi.1*sig.across.phi.2,
    rho.across.phi*sig.across.phi.1*sig.across.phi.2, sig.across.phi.2^2),
  ncol=2
)

#Correlation:
if(i == 1) {rho.across.gam <- rag} else {rho.across.gam <- 0}

#Covariance matrix:
Sig.across.gam <- matrix(
  c(sig.across.gam.1^2, rho.across.gam*sig.across.gam.1*sig.across.gam.2,
    rho.across.gam*sig.across.gam.1*sig.across.gam.2, sig.across.gam.2^2),
  ncol=2
)

#Correlation:
if(i == 1) {rho.within.phi <- rwp} else {rho.within.phi <- 0}

#Covariance matrix:
Sig.within.phi <- matrix(
  c(sig.within.phi.1^2, rho.within.phi*sig.within.phi.1*sig.within.phi.2,
    rho.within.phi*sig.within.phi.1*sig.within.phi.2, sig.within.phi.2^2),
  ncol=2
)

#Correlation:
if(i == 1) {rho.within.gam <- rwg} else {rho.within.gam <- 0}

#Covariance matrix:
Sig.within.gam <- matrix(
  c(sig.within.gam.1^2, rho.within.gam*sig.within.gam.1*sig.within.gam.2,
    rho.within.gam*sig.within.gam.1*sig.within.gam.2, sig.within.gam.2^2),
  ncol=2
)

# Draw the random effects:
b.0.phi <- rmnorm(n.patients, mean=rep(0,n.strains), varcov=Sig.across.phi) 
b.0.gam <- rmnorm(n.patients, mean=rep(0,n.strains), varcov=Sig.across.gam)

b.1t.phi <- rmnorm(sum(n.visits), mean=rep(0,n.strains), varcov=Sig.within.phi)
b.1t.gam <- rmnorm(sum(n.visits), mean=rep(0,n.strains), varcov=Sig.within.gam)

# Calculate phis and gammas:
phi <- NULL
gam <- NULL

# Formula from our model:
# (phi/gamma) <- global_mean + intercept_among + covariate_among + 
#                     intercept_within + covariate_within
for(j in 1:n.obs){
  phi[j] <- mean.phi + b.0.phi[Patient[j],Strain[j]] + beta.pat.phi[Strain[j]] * X.pat.vals[j] + 
    b.1t.phi[Visit[j],Strain[j]] + beta.time.phi[Strain[j]] * X.time.vals[j]
  gam[j] <- mean.gam + b.0.gam[Patient[j],Strain[j]] + beta.pat.gam[Strain[j]] * X.pat.vals[j] + 
    b.1t.gam[Visit[j],Strain[j]] + beta.time.gam[Strain[j]] * X.time.vals[j]
}


# First we need the initial occurrence probability for each strain for each patient.
# We assume this probability is not influenced by covariates, but this assumption
# could be relaxed in the future. Instead we assume the initial probabilities are
# generated randomly from a global occurrence probability.

# Calculate AntiLogit for phi and gam
lphi <- AntiLogit(phi)
lgam <- AntiLogit(gam)

psi <- NULL
for(j in 1:n.obs){
  
  if(Visit.Pat[j]==1){ # If it's the patient's first visit
    psi[j] <- AntiLogit(rnorm(1,mean.psi,.2))
  }else{
    psi[j] <- lphi[j-n.strains] * psi[j-n.strains] + lgam[j-n.strains] * (1 - psi[j-n.strains]) 
  }
  
}

Y_Alt <- NULL
for(i in 1:n.obs){
  Y_Alt[i] <- rbinom(1, 1, psi[i])
}



# Create a data.frame
Occ <- data.frame(Y_Alt, lphi, lgam, psi, Strain, Patient, Visit.Pat, X.time.vals, X.pat.vals)

t1 <- subset(Occ, Strain == 1)
t2 <- subset(Occ, Strain == 2)

plot(t1$lphi ~ t2$lphi)
plot(t1$psi ~ t2$psi)

# comb <- cbind(t2$psi, t1$psi, t1[,6])
# head(comb)
# colnames(comb) <- c("psi2", "psi1", "Pat")
# comb <- data.frame(comb)
# 
# library(ggplot2)
# ggplot(comb, aes(x=psi2, y=psi1, color=factor(Pat)))+
#   geom_point()+
#   theme_classic()

# Weird GLM
m1 <- glm(t1$Y_Alt ~ t2$Y_Alt, family="binomial")
summary(m1)

