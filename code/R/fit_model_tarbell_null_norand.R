library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("/scratch/jmihaljevic1/HPV/")

# Load the data:
load("data/full_HIM_data_10_strains_all_patients.rda")
n_strains <- stan_d$n_strains
n_patients <- stan_d$n_patient
n <- stan_d$n

# Initial values:
inits_f <- function(){
  list(betas_phi = array(rnorm(n_strains*n_strains,0,1), dim=c(n_strains,n_strains)),
       betas_gam = array(rnorm(n_strains*n_strains,0,1), dim=c(n_strains,n_strains)),
       betas_tbv_phi = rnorm(n_strains, 0, 1),
       betas_tbv_gam = rnorm(n_strains, 0, 1),
       alphas = rnorm(n_strains, -2, 1),
       abs_ystar = array(abs(rnorm(n*n_strains,0,2)), dim=c(n,n_strains))
       )
}

params <- c('betas_phi', 'phi_mean', 'phi_sd',
            'betas_gam', 'gam_mean', 'gam_sd',
            'betas_tbv_phi', 'tbv_phi_mean', 'tbv_phi_sd',
            'betas_tbv_gam', 'tbv_gam_mean', 'tbv_gam_sd',
            'alphas', 'alpha_mean', 'alpha_sd',
            'log_lik')

test <- stan('code/stan/twolevel_null_norand.stan',
             data = stan_d, chains = 1, iter = 10,
             init = inits_f,
             pars = params)

start <- Sys.time()
m_fit <- stan(fit = test,
              data = stan_d,
              init = inits_f,
              chains = 3, iter = 1500, warmup = 500,
              pars = params,
              control=list(max_treedepth=13))
end <- Sys.time()
time_taken = end - start
print(time_taken)

R_hats = summary(m_fit)$summary[,"Rhat"]
filename = "output/R_hats_null_norand.rds"
saveRDS(R_hats, file = filename)

posts = extract(m_fit, pars = params[1:15])
filename = "output/posts_full_null_norand.rds"
saveRDS(posts, file = filename)

# extract log_lik
log_lik = extract(m_fit, pars = "log_lik")$log_lik

# clear some memory:
rm(posts, R_hats, m_fit)

# Re-dimensionalize to a Sample x Observation matrix:
n_sample = 3000
dim(log_lik) = c(n_sample, n*n_strains)

# Now split the matrix into manageable pieces:
# Split by row (nrow = n_sample)
n_files = 10
max_row = n_sample / n_files
these_splits = split(1:n_sample, ceiling((1:n_sample)/max_row))

for(i in 1:length(these_splits)){
  
  lower = range(these_splits[[i]])[1]
  upper = range(these_splits[[i]])[2]
  
  saveRDS(log_lik[lower:upper,], 
          file = paste("output/loglik_full_null_norand_", i, ".rds", sep=""))
}

