library(RSQLite)
# library(rstan)
# rstan_options(auto_write = TRUE)
# options(mc.cores = parallel::detectCores())

setwd("/scratch/jmihaljevic1/HPV/")

# set chain ID to PID
n_row = 1000
n_col = 30525

log_lik = matrix(rnorm(n_row * n_col), ncol = n_col, nrow = n_row)
log_lik = t(log_lik)
log_lik_df = as.data.frame(log_lik)

setwd("/scratch/jmihaljevic1/HPV/output/")

db = RSQLite::dbConnect(SQLite(), dbname = "test_df.sqlite")
dbWriteTable(db,"log_lik",log_lik_df, append = T)

# n_col_max = 200
# col_indices <- c(1:ncol(log_lik_df))
# col_splits <- split(col_indices, ceiling(seq_along(col_indices)/n_col_max))
# 
# db = RSQLite::dbConnect(SQLite(), dbname = "test_df.sqlite")
# 
# for( i in c(1:length(col_splits))){
#   print(i)
#   table_name = paste0("loglik_",i)
#   print("a")
#   ind = unlist(col_splits[[i]]) 
#   print("b")
#   log_lik_df_subset <- log_lik_df[,ind]
#   print("c")
#   dbWriteTable(db,table_name,log_lik_df_subset, append=TRUE)
#   print("d")
# }
# 
# dbDisconnect(db)
