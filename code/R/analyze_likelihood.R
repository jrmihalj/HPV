# From outputted likelihood
library(rstan)
library(grid)
library(gridExtra)
library(RSQLite)
library(loo)

setwd("./output/")

files <- list.files(pattern = "norand_chain")
ll_all <- data.frame()
for( i in 1:length(files)){
  dbResultsFilename <- files[i]
  print(dbResultsFilename)
  db <- dbConnect(SQLite(), dbResultsFilename)
  loglik_tables <- dbListTables(db)[grep("loglik", dbListTables(db))]
  ll <- data.frame()
  for(j in c(1:length(loglik_tables))){
    print(loglik_tables[j])
    ll_sub <- dbReadTable(db,loglik_tables[j])
    if( j == 1){
      ll <- ll_sub
    }
    if(j > 1){
      ll <- cbind(ll,ll_sub)
    }
    rm(ll_sub)
  }
  if( i == 1){
    ll_all <- ll
  }
  if(i > 1){
    ll_all <- cbind(ll_all, ll)
  }
  rm(ll)
  dbDisconnect(db)
}

# WAIC and LOO-IC functions need an S x N matrix
# S = number posterior draws; N = number of data points/observations

# WAIC
waic_norand <- waic(t(as.matrix(ll_all)))
save(waic_norand, file = "waic_norand.rda")

# LOO-IC
looic_norand = loo(t(as.matrix(ll_all)))
save(looic_norand, file = "looic_norand.rda")

#------------------------------------------------------------------------
#------------------------------------------------------------------------

# Summarize WAIC

## Once all Lhoods have been calculated: ### 
waic_files <- list.files(pattern = "waic")
#waic_all <- data.frame()
for( i in 1:length(waic_files)){
  load(waic_files[i])
}

waic_tab =
  rbind(sapply(list(waic_null, waic_norand, waic_nofixed, waic_full), function(x) round(x$waic,2)),
      sapply(list(waic_null, waic_norand, waic_nofixed, waic_full), function(x) round(x$se_waic,2)))

colnames(waic_tab) = c("Null", "NoRand", "NoFixed", "Full")
rownames(waic_tab) = c("WAIC", "SE_WAIC")

t(waic_tab)

#------------------------------------------------------------------------
#------------------------------------------------------------------------

# Summarize LOO-IC

## Once all Lhoods have been calculated: ### 
loo_files <- list.files(pattern = "loo")
#loo_all <- data.frame()
for( i in 1:length(loo_files)){
  load(loo_files[i])
}

loo_tab =
  rbind(sapply(list(looic_null, looic_norand, looic_nofixed, looic_full), function(x) round(x$looic,2)),
        sapply(list(looic_null, looic_norand, looic_nofixed, looic_full), function(x) round(x$se_looic,2)))

colnames(loo_tab) = c("Null", "NoRand", "NoFixed", "Full")
rownames(loo_tab) = c("loo-ic", "SE_loo-ic")

t(loo_tab)
