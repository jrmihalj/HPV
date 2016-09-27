
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
simulate_data = FALSE 
use_complete_data = TRUE 
args = commandArgs(trailingOnly=TRUE)
k = as.numeric(args[1])


if(simulate_data){
  source('sim.R')
  stan_d <- list(n = n, 
                 d = d, 
                 y = y,
                 eta = 2,
                 patient = patient, n_patient = n_patients,
                 n_visit = n_visits, visit = visit,
                 dir_prior = c(1, 1))
}

if(!simulate_data){
  load("test_data_HIM_full.rda")
  d = stan_d$d
  n_patients = stan_d$n_patient
  n_complete_obs = stan_d$n
  tbv = stan_d$tbv
  
  if(!use_complete_data){
    # find missing visits
    missing <- which(is.na(rowSums(stan_d$y)))
    # find number and indices of complete observations 
    n_complete_obs <- stan_d$n - length(missing)
    stan_d$n_complete_obs <- n_complete_obs
    complete_obs <- which(stan_d$visit != -1)
    stan_d$complete_obs <- complete_obs
    # stan doesn't accept "NA" so set missing 0/1 observations to be .5 
    # (these will not be counted in the model but they need a placeholder)
    y = stan_d$y
    y[is.na(y)] <- .5
    stan_d$y = y
  }
  
}

inits_f <- function(){
  list(betas_phi = array(rnorm(d*d,0,1), dim=c(d,d)),
       betas_gam = array(rnorm(d*d,0,1), dim=c(d,d)),
       #alphas = rnorm(d, 0, 2),
       betas_tbv_phi = array(rnorm(d,0,1)),
       betas_tbv_gam = array(rnorm(d,0,1)),
       e_patient = array(rnorm(n_patients*d,0,1), dim=c(n_patients, d)),
       abs_ystar = array(abs(rnorm(n_complete_obs*d,0,2)), dim=c(n_complete_obs,d))
       )
}

params <- c('Rho_patient', 'Rho_visit',
            'sd_visit', 'sd_patient', 'var_mat',
            'betas_phi', 'betas_gam','betas_tbv_phi','betas_tbv_gam', 'alphas')

fit_model <- function(i){

  test <- stan('twolevel.stan',
                data = stan_d, chains = 1, iter = 10,
                init = inits_f,
                pars = params)
  
  start <- Sys.time()
  m_fit <- stan(fit = test,
                data = stan_d,
                init = inits_f,
                chains = 1, iter = 2000, warmup = 1000, thin=3,
                pars = params,
                control=list(max_treedepth=13))
  end <- Sys.time()
  time_taken = end - start
  print(time_taken)
  
  filename <- paste0("fit_chain_", i, "10_strains_tbv.rda")
  save(m_fit, file = filename)
}

test <- fit_model(k)


