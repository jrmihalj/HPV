library(RSQLite)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("/scratch/jmihaljevic1/HPV/")

# set chain ID to PID
k = Sys.getpid()
print(k)
dbResultsFilename <- paste0("./output/fit_null_null_test_chain_",k,".sqlite")

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

params <- c('betas_tbv_phi', 'tbv_phi_mean', 'tbv_phi_sd',
            'betas_tbv_gam', 'tbv_gam_mean', 'tbv_gam_sd',
            'alphas', 'alpha_mean', 'alpha_sd',
            'log_lik')

test <- stan('code/stan/twolevel_null_null.stan',
             data = stan_d, chains = 1, iter = 10,
             init = inits_f,
             pars = params)

start <- Sys.time()
m_fit <- stan(fit = test,
              data = stan_d,
              init = inits_f,
              chains = 1, iter = 1500, warmup = 500, #thin = 7,
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
posts =extract(m_fit, pars = params[1:9])

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
log_lik = extract(m_fit, pars = "log_lik")$log_lik

# Make it an observation x sample data frame
log_lik = t(log_lik)
log_lik_df = as.data.frame(log_lik)

n_col_max = 200
col_indices <- c(1:ncol(log_lik_df))
col_splits <- split(col_indices, ceiling(seq_along(col_indices)/n_col_max))

db <- dbConnect(SQLite(), dbResultsFilename)
for( i in c(1:length(col_splits))){
  table_name <- paste0("loglik_",i)
  ind <- unlist(col_splits[[i]]) 
  log_lik_df_subset <- log_lik_df[,ind]
  dbWriteTable(db,table_name,log_lik_df_subset, append=TRUE)
  print(paste("Done with split", i))
}
dbDisconnect(db)
