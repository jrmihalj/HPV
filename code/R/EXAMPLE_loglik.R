setwd("./output/")

library(rstan)
library(RSQLite)
library(loo)

db_names = system("ls *null_chain*", intern = T)

results = NULL

for(i in 1:length(db_names)){
  
  db = dbConnect(SQLite(), db_names[i])
  results_temp = dbReadTable(db, "loglik")
  results_temp = as.matrix(results_temp)
  
  if(i == 1){
    results = results_temp
  }else{
    results = cbind(results, results_temp)
  }
  
  if(i == length(db_names)){
    rm(results_temp)
  }
  
}


waic_test = waic(results)
loo_test = loo(results)
str(loo_test)
