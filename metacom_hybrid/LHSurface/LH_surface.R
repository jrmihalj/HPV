# Exploring the Likelihood Surface of betas and correlations
#################################################################
#################################################################

# First, set wd and source some data:
setwd("./metacom_hybrid/LHSurface")
source("../sim_data_metah.R")

# Create a quick function to generate the likelihood:
# For now, have the correlations known
LH_calc <- function(betas, fix.beta=F, fix.rho=T){
  
  #######  #######  #######
  # PHI and GAMMA betas #
  #######  #######  #######
  
  if(fix.beta==T){
    bp_temp <- beta_phi_intxn;
    bg_temp <- beta_gam_intxn;
  }else{
    bp_temp <- matrix(c(betas[1], betas[2],
                        betas[3], betas[4]), ncol=2, byrow=T)
    bg_temp <- matrix(c(0, betas[5],
                        betas[6], 0), ncol=2, byrow=T)
  }
  
  if(fix.rho==F){ # Then we're simulating different correlations each time
    Rp_temp <- genPositiveDefMat(dim=m, #Number of columns/Ro_tempws
                            covMethod = 'onion', 
                            rangeVar = c(1, 1), # Range of variances
                            eta=2)$Sigma
    
    # observation-level ranefs
    ep_temp <- matrix(nrow = n_site, ncol = m)
    for (i in 1:n_site){
      ep_temp[i, ] <- rnorm(m) %*% t(chol(Rp_temp))
    }
    
    # observation-level corr matrix
    # Ro_temp <- genPositiveDefMat(dim=m, #Number of columns/rows
    #                         covMethod = 'onion', 
    #                         rangeVar = c(1, 1), # Range of variances
    #                         eta=2)$Sigma
    
    Ro_temp <- matrix(c(1,0,0,1), ncol=m, byrow=T)
    
    # observation-level ranefs
    eo_temp <- matrix(nrow = n, ncol = m)
    for (i in 1:n){
      eo_temp[i, ] <- rnorm(m) %*% t(chol(Ro_temp))
    }
    
    # Add the ranef
    e_all_temp <- matrix(nrow = n, ncol = m)
    for(i in 1:n){
      e_all_temp[i, ] <- ep_temp[site[i], ] + eo_temp[i, ]
    }
    
    # Normalize
    e_all_temp <- (e_all_temp - mean(e_all_temp)) / sd(e_all_temp)
    
  }else{
    
    e_all_temp <- e_all
    
  }

  ##############
  # Occurrence #
  ##############
  
  z_temp <- matrix(NA, nrow=n, ncol=m)
  y_temp <- matrix(NA, nrow=n, ncol=m)
  
  # first observation:
  z_temp[time==1, ] =  qlogis(pnorm(e_all_temp[time==1, ])) #only based on correlations
  y_temp[time == 1, ] <- ifelse(z_tot[time == 1, ] > 0, 1, 0)
  
  # subsequent observations:
  for (i in 2:n_timesteps) {
    #Persistence:
    z_temp[time == i, ] <- (y_temp[time == i - 1, ]) * (y_temp[time == i - 1, ] %*% bp_temp) +
      #Colonization:
      (1 - y_temp[time == i - 1, ]) * (y_temp[time == i - 1, ] %*% bg_temp) +
      #Correlated and nested random effects:
      qlogis(pnorm(e_all_temp[time==i, ]))
    
    y_temp[time == i, ] <- ifelse(z_temp[time == i, ] > 0, 1, 0)
  }
  
  prob_temp <- NULL; LH_all <- NULL; LH <- NULL;
  
  prob_temp <- pnorm(z_temp)
  LH_all <-  dbinom(y, 1, prob=prob_temp, log=TRUE)
  #Fix Inf/-Inf values (e.g. y=1, prob=0)
  LH_all <- ifelse(is.infinite(LH_all)==T, -10, LH_all) # Inf means very low probability
  LH <- sum(LH_all)
  
  if(fix.rho==F){
    cat(paste(LH,Rp_temp[1,2],Ro_temp[1,2], sep=","),paste("\n"),file="LHSurfaceRho.csv", append=T, sep="")
  }
  
  return(LH)
}

# Test function:
LH_calc(rep(0,6), fix.beta = F, fix.rho = T)

#################################################################
#################################################################
true_b <- c(beta_phi_intxn[1,],beta_phi_intxn[2,],beta_gam_intxn[1,2],beta_gam_intxn[2,1])
#lower <- round(min(true_b)-0.1, 1); upper <- round(max(true_b)+0.1, 1);
lower <- -5; upper <- 5;
val <- lower; 
b_temp <- rep(lower, 6);
n.sweeps <- 10;
sweep <- 1;

# Conduct a grid search over each beta parameter:
cat("LH,b1,b2,b3,b4,b5,b6\n", file="LHSurfaceBeta.csv")
while(sweep <= n.sweeps){
  starter <- NULL
  starter <- sample.int(6, 1)
  this <- starter
  while(this < (6 + starter)){#number of params
    PC <- 1 + this %% 6
    val <- lower
    while(val <= upper){
      b_temp <- runif(6, lower, upper)
      b_temp[PC] <- val
      
      LH_temp <- NULL
      LH_temp <- LH_calc(b_temp, fix.beta = F, fix.rho = T)
      
      cat(paste(LH_temp),",",paste(b_temp,collapse=",",sep=","),paste("\n"), 
          file="LHSurfaceBeta.csv", append=T, sep="")
      val <- val + runif(1,0.01,.1)
    }
    
    this <- this + 1
  }
  sweep <- sweep + 1
}

#################################################################
#################################################################
library(readr)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)

# Read in the data now
LHdf <- read_csv(file="LHSurfaceBeta.csv", col_names=T)

# Remove NA vals
LHdf <- filter(LHdf, is.infinite(LH)==F)

# Transform likelihood (now smaller values are better)
LHdf$LogLH <- log(-2*LHdf$LH)

# Plot it
these <- c(2:7)
all_plots <- list()
counter <- 1
param_names <- c("bp11", "bp12", "bp21", "bp22", "bg12", "bg21")
#quartz(height=6, width=6)
for(i in these){
  for(j in these){
    if( (i == j) | (i > j)){
    
      p <- nullGrob()

      if(i==j){ 
        p <- textGrob(param_names[i-1])
      }
      
      all_plots[[counter]] <- p
      counter <- counter + 1
      
    } else {
      all_plots[[counter]] <- ggplot(LHdf)+
                                aes_string(x=colnames(LHdf)[j], y=colnames(LHdf)[i])+
                                theme_classic()+
                                geom_point(aes(color=LogLH), size=.2, alpha=.5)+
                                labs(x="",y="", title="")+
                                scale_y_continuous(limits=c(lower,upper), breaks=c(lower,0,upper))+
                                scale_x_continuous(limits=c(lower,upper), breaks=c(lower,0,upper))+
                                scale_color_gradient("LHood", low="blue", high="yellow")+
                                annotate("point", x=true_b[j-1], y=true_b[i-1], color="red")+
                                guides(color=F)
      counter <- counter + 1
    }
  }
}

# Plot on one page
quartz(height=10,width=10)
plot_it <- arrangeGrob(grobs=all_plots, ncol=6)
grid.arrange(plot_it)

#################################################################
#################################################################
#################################################################
#################################################################

# Now look at the Rho's 

library(readr)
library(dplyr)
library(ggplot2)

n.samples <- 2000
cat("LH,Rp,Ro\n", file="LHSurfaceRho.csv", sep="")

i = 1;
while(i <= n.samples){
  LH_calc(fix.beta = T, fix.rho = F)
  i <- 1 + i
}

# Pull in the data:
LH_Rho <- read_csv(file="LHSurfaceRho.csv")

quartz(height=5, width=8)
par(mfrow=c(1,2),
    mai=c(0.5,0.5,0,0),
    omi=c(0,.25,0,0))
plot(LH_Rho$LH ~ LH_Rho$Rp, axes="F", ylim=c((min(LH_Rho$LH)-20), (max(LH_Rho$LH)-20)))
axis(1, at=c(-1, 0, 1), mgp=c(3,.4,0))
axis(2, at=c(round((min(LH_Rho$LH)-20)), round(max(LH_Rho$LH)-20)), mgp=c(3,.4,0))
par(xpd=F)
abline(v=Rp[1,2], col="red", lwd=1.2)
mtext("Rp",1,outer=F,line=1.2)

par(mai=c(0.5,0,0,0.2), xpd=T)
plot(LH_Rho$LH ~ LH_Rho$Ro, axes="F", ylim=c((min(LH_Rho$LH)-20), (max(LH_Rho$LH)-20)))
axis(1, at=c(-1, 0, 1), mgp=c(3,.4,0))
par(xpd=F)
abline(v=Ro[1,2], col="red", lwd=1.2)
mtext("Ro",1,outer=F,line=1.2)
mtext("Likelihood", side=2, outer=T, line=0)

quartz(height=5, width=5)
ggplot(LH_Rho, aes(Ro, Rp))+
  geom_point(aes(color=LH), size=1)+
  theme_classic()+
  scale_y_continuous(limits=c(-1,1), breaks=c(-1,0,1))+
  scale_x_continuous(limits=c(-1,1), breaks=c(-1,0,1))+
  scale_color_gradient("log(-2*LH)", low="yellow", high="blue")+
  annotate("point", x=Ro[1,2], y=Rp[1,2], color="red", size=4)+
  labs(x="Observation-level Rho", y="Patient-level Rho")+
  guides(color=F)
  
  