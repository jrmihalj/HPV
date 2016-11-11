# Function calculates waic for occupancy model 
# input mcmc.list with posterior sims and data.frame for model input
# output list: log posterior predictive density, p_WAIC2, WAIC, mean likelihood
calc_waic_ZIM <- function(posts, data){
  with(data,{
    
    ######################
    # STRUCTURE THE DATA #
    ######################
    
    
    
    N_all <- N_obs1 + N_obs2 + N_obs3
    Nsamples <- dim(posts[[1]])[1]
    
    L <- array(dim=c(N_all, Nsamples))
    
    L_bar <- rep(NA, N_all)
    var_LL <- rep(NA, N_all)
    
    ##########################
    # LIKELIHOOD CALCULATION #
    ##########################
    

    ##################
    # SUMMARIZE WAIC #
    ##################
    
    for(i in 1:N_all){
      L_bar[i] <- mean(exp(L[i, ]))
      var_LL[i] <- var(L[i, ])
    }  
    
    
    lppd <- sum(log(L_bar), na.rm=T)
    p_WAIC <- sum(var_LL, na.rm=T)
    WAIC <- -2 * (lppd - p_WAIC)
    return(list(lppd=lppd, p_WAIC=p_WAIC, WAIC=WAIC, L_bar=L_bar))
  })
}