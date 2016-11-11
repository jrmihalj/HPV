library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("~/Sylvia_HPV")

# Load the data:
load("data/full_HIM_data_10_strains_all_patients.rda")
n_strains <- stan_d$n_strains
n_patients <- stan_d$n_patient
n <- stan_d$n

# Initial values:
inits_f <- function(){
  list(betas_tbv_phi = rnorm(n_strains, 0, 1),
       betas_tbv = rnorm(n_strains, 0, 1),
       alphas = rnorm(n_strains, -2, 1),
       abs_ystar = array(abs(rnorm(n*n_strains,0,2)), dim=c(n,n_strains))
       )
}

params <- c('betas_tbv_phi', 'betas_tbv_gam', 'alphas', 'log_lik')

test <- stan('code/stan/twolevel_null_null.stan',
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

filename <- "output/fit_full_null_null.rda"
save(m_fit, file = filename)

