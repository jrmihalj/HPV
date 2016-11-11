# Function calculates waic for occupancy model 
# input mcmc.list with posterior sims and data.frame for model input
# output list: log posterior predictive density, p_WAIC2, WAIC, mean likelihood
calc_waic_ZIM <- function(posts, data){
  with(data,{
    
    ######################
    # STRUCTURE THE DATA #
    ######################
    Y <- c(off1,off2,off3)
    
    mus <- rbind(t(posts$mu1), t(posts$mu2), t(posts$mu3))
    
    if(type %in% c("ZINB_C", "ZHNB_C")){
      mu0s <- rbind(t(posts$mu1_0), t(posts$mu2_0), t(posts$mu3_0))
    }
    
    if(MF){ #Dif phi for MF in experiment 1 and 2
      phis <- array(0, dim = c(dim(posts$phi)[1], 3, 2))
      for(e in 1:2){
        for(s in 1:2){
          phis[,s,e] <- posts$phi[,s,e]
        }
      }
      phi_line <- matrix(posts$phi_line, ncol=1)
      phis[,3,1] <- phi_line # not sex specific
      phis[,3,2] <- phi_line
      
      this_phi <- rep(1:3, c(N_obs1,N_obs2,N_obs3))
      this_sex <- c(data$Sex, data$Sex2, rep(1, N_obs3))
      
    }else{
      phis <- rbind(t(posts$phi[,1]), t(posts$phi[,2]), t(posts$phi[,3]))
      this_phi <- rep(1:3, c(N_obs1,N_obs2,N_obs3))
    }
    
    
    N_all <- N_obs1 + N_obs2 + N_obs3
    Nsamples <- dim(posts[[1]])[1]
    
    L <- array(dim=c(N_all, Nsamples))
    
    L_bar <- rep(NA, N_all)
    var_LL <- rep(NA, N_all)
    
    ##########################
    # LIKELIHOOD CALCULATION #
    ##########################
    for (i in 1:N_all){
      
      
      if(type=="NB"){
        if(MF){
          L[i, ] <- dnbinom(rep(Y[i], Nsamples), size=phis[, this_phi[i], this_sex[i]], mu=exp(mus[i,]), log=TRUE)
        }else{
          L[i, ] <- dnbinom(rep(Y[i], Nsamples), size=phis[this_phi[i], ], mu=exp(mus[i,]), log=TRUE)
        }
      }
      if(type=="ZINB_NC"){
        if(MF){
          L[i, ] <- dzinb(rep(Y[i], Nsamples), k=phis[, this_phi[i], this_sex[i]], lambda=exp(mus[i, ]), omega=plogis(posts$theta), log=TRUE)
        }else{
          L[i, ] <- dzinb(rep(Y[i], Nsamples), k=phis[this_phi[i], ], lambda=exp(mus[i, ]), omega=plogis(posts$theta), log=TRUE)
        }
      }
      if(type=="ZINB_C"){
        if(MF){
          L[i, ] <- dzinb(rep(Y[i], Nsamples), k=phis[, this_phi[i], this_sex[i]], lambda=exp(mus[i, ]), omega=plogis(mu0s[i, ]), log=TRUE)
        }else{
          L[i, ] <- dzinb(rep(Y[i], Nsamples), k=phis[this_phi[i], ], lambda=exp(mus[i, ]), omega=plogis(mu0s[i, ]), log=TRUE) 
        }
      }
      
      if(type=="ZHNB_NC"){
        
        if(Y[i]==0){
          L[i, ] <- dbinom(rep(Y[i], Nsamples), 1, prob=1-posts$theta, log=TRUE)
        }else{#Truncated neg binom distribution
          if(MF){
            L[i, ] <- dposnegbin(rep(Y[i], Nsamples), size=phis[, this_phi[i], this_sex[i]], munb=exp(mus[i, ]), log=TRUE)
          }else{
            #L[i, ] <- dnbinom(rep(Y[i], Nsamples), size=phis[this_phi[i], ], mu=exp(mus[i, ]), log=TRUE)
            L[i, ] <- dposnegbin(rep(Y[i], Nsamples), size=phis[this_phi[i], ], munb=exp(mus[i, ]), log=TRUE)
          }
          
        }
        
      }
      
      if(type=="ZHNB_C"){
        
        if(Y[i]==0){
          L[i, ] <- dbinom(rep(Y[i], Nsamples), 1, prob=1-plogis(mu0s[i, ]), log=TRUE)
        }else{#Truncated neg binom distribution
          if(MF){
            L[i, ] <- dposnegbin(rep(Y[i], Nsamples), size=phis[, this_phi[i], this_sex[i]], munb=exp(mus[i, ]), log=TRUE)
          }else{
            #L[i, ] <- dnbinom(rep(Y[i], Nsamples), size=phis[this_phi[i], ], mu=exp(mus[i, ]), log=TRUE)
            L[i, ] <- dposnegbin(rep(Y[i], Nsamples), size=phis[this_phi[i], ], munb=exp(mus[i, ]), log=TRUE) 
          }
          
        }
        
      }
      
      if(type=="PNBH_NC"){
        
        if(Y[i]<=p_hurdle){
          L[i, ] <- dpois(rep(Y[i], Nsamples), lambda = posts$lambda, log=TRUE)
        }else{#Truncated neg binom distribution
          #LogLike = loglike(y) - log(1 -like11)?
          L[i, ] <- dnbinom(rep(Y[i], Nsamples), size=phis[this_phi[i], ], mu=exp(mus[i, ]), log=TRUE) -
                      log( 1 - dnbinom(rep((p_hurdle+1), Nsamples), size=phis[this_phi[i], ], mu=exp(mus[i, ]), log=F))
        }
        
      }
      
    }

    ##################
    # SUMMARIZE WAIC #
    ##################
    
    for(i in 1:N_all){
      L_bar[i] <- mean(exp(L[i, ]))
      var_LL[i] <- var(L[i, ])
    }  
    
    if(remove){
      these = which(log10(L_bar) < -5)
      L_bar[these] = NA
      var_LL[these] = NA
    }
    
    lppd <- sum(log(L_bar), na.rm=T)
    p_WAIC <- sum(var_LL, na.rm=T)
    WAIC <- -2 * (lppd - p_WAIC)
    return(list(lppd=lppd, p_WAIC=p_WAIC, WAIC=WAIC, L_bar=L_bar))
  })
}