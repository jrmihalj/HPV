setwd("./output/")

library(rstan)
library(RSQLite)
library(loo)

db = dbConnect(SQLite(), "output/fit_null_null_chain_14681.sqlite")
results = dbReadTable(db, "loglik")
results = as.matrix(results)

waic_test = waic(results)
loo_test = loo(results)
str(loo_test)
