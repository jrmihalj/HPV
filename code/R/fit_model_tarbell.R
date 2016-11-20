library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
simulate_data = FALSE 
use_complete_data = TRUE 

setwd("/scratch/jmihaljevic1/HPV/")

if(!simulate_data){
  load("data/full_HIM_data_10_strains_all_patients.rda")
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
            'betas_tbv_phi', 'betas_tbv_gam', 'alphas', 'log_lik')


test <- stan('code/stan/twolevel.stan',
             data = stan_d, chains = 1, iter = 10,
             init = inits_f,
             pars = params)

start <- Sys.time()
m_fit <- stan(fit = test,
              data = stan_d,
              init = inits_f,
              chains = 3, iter = 3000, warmup = 2000,
              pars = params,
              control=list(max_treedepth=13))
end <- Sys.time()
time_taken = end - start
print(time_taken)

object.size(m_fit)

R_hats = summary(m_fit)$summary[,"Rhat"]
filename = "output/R_hats_full.rds"
saveRDS(R_hats, file = filename)

posts = extract(m_fit)
filename <- "output/fit_full.rds"
saveRDS(posts, file = filename)


