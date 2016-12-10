library(rstan)
library(RSQLite)
rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores())
simulate_data = FALSE 
use_complete_data = TRUE 

# Get chain id from batch script ##
args = commandArgs(trailingOnly=TRUE)
k = as.numeric(args[1]) 
dbResultsFilename <- paste0("results_chain_",k,".sqlite")

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
            'betas_phi', 'phi_mean', 'phi_sd',
            'betas_gam', 'gam_mean', 'gam_sd',
            'betas_tbv_phi', 'tbv_phi_mean', 'tbv_phi_sd',
            'betas_tbv_gam', 'tbv_gam_mean', 'tbv_gam_sd',
            'alphas', 'alpha_mean', 'alpha_sd',
            'log_lik')


test <- stan('code/stan/twolevel.stan',
             data = stan_d, chains = 1, iter = 10,
             init = inits_f,
             pars = params)

start <- Sys.time()
m_fit <- stan(fit = test,
              data = stan_d,
              init = inits_f,
              chains = 1, iter = 5000, warmup = 2000, thin = 7,
              pars = params,
              control=list(max_treedepth=13))
end <- Sys.time()
time_taken = end - start
print(time_taken)

# Save Rhat 
R_hat_df <- data.frame(parameter = rownames(as.data.frame((summary(m_fit)$summary[,"Rhat"]))),
                       R_hat = as.numeric((summary(m_fit)$summary[,"Rhat"])))
table_name <- "R_hat"
db <- dbConnect(SQLite(), dbResultsFilename)
dbWriteTable(db,table_name,R_hat_df, append=TRUE)
dbDisconnect(db)


# Get posterior parameter values 
posts =extract(m_fit, pars = params[1:14])

for( i in 1:length(posts)){
  post_df <- as.data.frame(posts[[i]])
  names(post_df) <- paste0(names(posts)[[i]], names(post_df))
  post_df$chain = k
  post_df$sample = c(1:nrow(post_df))
  table_name <- paste0(names(posts)[[i]])
  db <- dbConnect(SQLite(), dbResultsFilename)
  dbWriteTable(db,table_name,post_df, append=TRUE)
  dbDisconnect(db)
}

# extract log_lik
log_lik_df = as.data.frame(extract(m_fit, pars = "log_lik")$log_lik)
log_lik_df$chain = k
table_name <- "loglik"
db <- dbConnect(SQLite(), dbResultsFilename)
dbWriteTable(db,table_name,log_lik_df, append=TRUE)
dbDisconnect(db)

