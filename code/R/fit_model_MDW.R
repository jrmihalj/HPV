library(rstan)
library(parallel)
rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores())
simulate_data = FALSE 
use_complete_data = TRUE 


if(simulate_data){
  source('code/R/sim.R')
  stan_d <- list(n = n, 
                 n_strains = n_strains, 
                 y = y,
                 tbv = tbv,
                 eta = 2,
                 patient = patient, n_patient = n_patients,
                 n_visit_max = n_visit_max, visit = visit,
                 dir_prior = c(1, 1))
}

if(!simulate_data){
  load("data/test_data_HIM_200_pat_2_strains.rda") # does not exist?
  n_strains <- stan_d$n_strains
  n_patients <- stan_d$n_patient
  n <- stan_d$n
}

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



fit_model <- function(i){
  
 # test <- stan('twolevel.stan',
  #             data = stan_d, chains = 1, iter = 10,
   #            init = inits_f,
    #           pars = params)
  
  start <- Sys.time()
  m_fit <- stan('code/stan/twolevel.stan',
		#fit = test,
                data = stan_d,
                init = inits_f,
                chains = 1, iter = 1500, warmup = 500,
                pars = params,
                control=list(max_treedepth=13))
  end <- Sys.time()
  time_taken = end - start
  print(time_taken)
  
  filename <- paste0("output/fit_chain_", i, "_2_strains_tbv_200_pat.rda")
  save(m_fit, file = filename)
}

mclapply(c(1:3), FUN = fit_model, mc.cores = 3)

