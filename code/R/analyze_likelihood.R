# From outputted likelihood
library(rstan)
library(grid)
library(gridExtra)
library(RSQLite)
library(loo)

setwd("./output/")

files <- list.files(pattern = "fit_null_nofixed")
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

# WAIC function needs an S x N matrix
# S = number posterior draws; N = number of data points/observations
waic_nofixed <- waic(t(as.matrix(ll_all)))
save(waic_nofixed, file = "waic_null_nofixed.rda")

## Once all Lhoods have been calculated: ### 
waic_files <- list.files(pattern = "waic")
#waic_all <- data.frame()
for( i in 1:length(waic_files)){
  load(waic_files[i])
  #identifier <- strsplit(strsplit(waic_files[i],"[_]")[[1]][3],"[.]")[[1]][1]
  #df <- data.frame(waic = waic$waic, se = waic$se_waic, model = identifier)
  #waic_all <- rbind(waic_all, df)
}

waic_tab =
  rbind(sapply(list(waic_null, waic_norand, waic_nofixed), function(x) round(x$waic,2)),
      sapply(list(waic_null, waic_norand, waic_nofixed), function(x) round(x$se_waic,2)))

colnames(waic_tab) = c("Null", "NoRand", "NoFixed")
rownames(waic_tab) = c("WAIC", "SE_WAIC")

waic_tab
