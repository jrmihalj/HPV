# fitting the simpler model
library(rstan)
library(ggplot2)
library(reshape)
library(RSQLite)
dbResultsFilename <- "results_hybrid.sqlite"
args = commandArgs(trailingOnly=TRUE)
k = as.numeric(args[1])

m <- 2                   # species
n_timesteps <- c(3,5,10,20,30)          # visits/repeat observations
n_site <- c(20,40,60,80,100,150,200,300,400,500) 
version_corr <- c(1,2,3)

sweep_df <- data.frame( n_time = rep(n_timesteps, each = length(n_site)), n_s = rep(n_site, each = length(version_corr)), version = rep(version_corr, times = length(n_timesteps)*length(n_site)))
n_timesteps <- sweep_df[k,]$n_time
n_site <- sweep_df[k,]$n_s
Rp_filename <- paste0("Rp_",sweep_df[k,]$version,".csv" )
Ro_filename <- paste0("Ro_",sweep_df[k,]$version,".csv" )


# locations/patients
n <- n_timesteps * n_site # number of observations


source('sim_data_metah.R')

stan_d <- list(n_unit = n_site, 
               n_time = n_timesteps,
               n = n, 
               m = m,
               unit = site,
               time = time,
               y = y,
               y_mat = y)
m_init <- NULL
m_init <- stan('metah.stan', data=stan_d, iter=10, chains=1)

m_fit <- stan(fit = m_init, data=stan_d, warmup=100, iter=2000, chains=2, cores=2, 
              pars=c('R_p', 'R_o', "beta_phi", "beta_gam"),
              control=list(adapt_delta=0.9))

#traceplot(m_fit, pars=c('lp__'))

post <- extract(m_fit)[1:4]
post_dfs <- lapply(post,as.data.frame)
for( i in 1:length(post_dfs)){
  post_dfs[[i]]$trial <- k
}
tableNames <- names(post)

db <- dbConnect(SQLite(), dbResultsFilename)
for( i in 1:length(tableNames)){
  df <- post_dfs[[i]]
  tableName <- tableNames[[i]]
  dbWriteTable(db,tableName,df, append=TRUE)
}
dbDisconnect(db)

# # evaluate recovery of R_p
# traceplot(m_fit, pars = 'cor_patient', inc_warmup = TRUE)
# Rp_df <- data.frame(Rp = c(Rp), y = ((m*2)^2):1)
# plot(m_fit, pars = 'cor_patient') + 
#   geom_point(data = Rp_df, aes(x = Rp, y = y), 
#              size = 4, col = 'blue') + 
#   ggtitle('Recovery of patient-level correlations')
# 
# # evaluate recovery of R_o
# traceplot(m_fit, pars = 'cor_time', inc_warmup = TRUE)
# Ro_df <- data.frame(Ro = c(Ro), y = ((m*2)^2):1)
# plot(m_fit, pars = 'cor_time') + 
#   geom_point(data = Ro_df, aes(x = Ro, y = y), 
#              size = 4, col = 'blue') + 
#   ggtitle('Recovery of observation-level correlations')
# 
# # evaluate recovery of z_tot (psi)
# mean_psi <- apply(post$psi, c(2,3), mean)
# 
# plot(mean_psi[,2] ~ z_tot[,2], xlab="data", ylab="model")
# abline(b=1, a=0)
