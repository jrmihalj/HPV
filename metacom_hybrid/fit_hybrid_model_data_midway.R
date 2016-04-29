# fitting the simpler model
library(rstan)
library(ggplot2)
library(reshape)
library(RSQLite)
dbResultsFilename <- "results_hybrid_4_29_200_pat_LKJ10.sqlite"
#source('HIM_data_formatting.R')
load("test_data_HIM.rda")

m_init <- NULL
m_init <- stan('metah.stan', data=stan_d, iter=10, chains=1)

m_fit <- stan(fit = m_init, data=stan_d, warmup=100, iter=2000, chains=4, cores=4, 
              pars=c('R_p', 'R_o', "beta_phi", "beta_gam"),
              control=list(adapt_delta=0.9))

summary <- as.data.frame(summary(m_fit)$summary)


traceplot(m_fit, pars=c('lp__'))

post <- extract(m_fit)[1:4]
post_dfs <- lapply(post,as.data.frame)
tableNames <- names(post)

db <- dbConnect(SQLite(), dbResultsFilename)
for( i in 1:length(tableNames)){
  df <- post_dfs[[i]]
  tableName <- tableNames[[i]]
  dbWriteTable(db,tableName,df, overwrite=TRUE)
}
dbWriteTable(db,"summary",summary, overwrite=TRUE)
dbDisconnect(db)


