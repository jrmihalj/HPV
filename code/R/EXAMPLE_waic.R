setwd("./output/")

library(loo)

# Import posterior samples and Rhat values:
posts = readRDS("fit_full_null_null.rds")
Rhats = readRDS("R_hats_null_null.rds")
# Check Rhats, generally:
hist(Rhats)

# Need to put the log_lik into SxN matrix...
log_lik = posts$log_lik
these_dims = dim(log_lik)
dim(log_lik) = c(these_dims[1], these_dims[2]*these_dims[3])

# Calc WAIC
waic_out = waic(log_lik)
waic_out$waic
