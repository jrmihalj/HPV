# fitting the simpler model
library(rstan)
library(ggplot2)
library(reshape)
library(RSQLite)
dbResultsFilename <- "results.sqlite"
args = commandArgs(trailingOnly=TRUE)
k = as.numeric(args[1])

m <- 2                   # species
n_timesteps <- c(3,5,10,20,30)          # visits/repeat observations
n_site <- c(20,40,60,80,100,150,200,300,400,500)  

sweep_df <- data.frame( n_time = rep(n_timesteps, each = length(n_site)), n_s = rep(n_site, times = length(n_timesteps)))
n_timesteps <- sweep_df[k,]$n_time
n_site <- sweep_df[k,]$n_s

# locations/patients

n <- n_timesteps * n_site # number of observations


source('sim_data_meta.R')

stan_d <- list(n_patients = n_site, 
               n_obs = n, 
               n_strains = m,
               patient = site,
               visit = time,
               z = y)

m_init <- NULL
m_init <- stan('longitudinal_model.stan', data=stan_d, iter=10, chains=1)
m_fit <- stan(fit = m_init, data=stan_d, warmup=400, iter=2000, chains=2, cores=2, 
              pars=c('cor_patient', 'cor_time', 'psi'))

#traceplot(m_fit, pars=c('lp__'))

post <- extract(m_fit)[1:3]
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
