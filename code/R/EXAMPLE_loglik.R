n_sample = 3000
n = 30525
n_strains = 10

# Simulate the same sized array:
log_lik = array(rnorm(n_sample * n * n_strains), dim = c(n_sample, n, n_strains))

# Re-dimensionalize to n_sample x n_total_obs matrix
dim(log_lik) = c(n_sample, n*n_strains)

# Split the matrix into manageable pieces, and save as .rds
# Split by row (nrow = n_sample)
n_files = 10
max_row = n_sample / n_files
these_splits = split(1:n_sample, ceiling((1:n_sample)/max_row))

for(i in 1:length(these_splits)){
  
  lower = range(these_splits[[i]])[1]
  upper = range(these_splits[[i]])[2]
  
  saveRDS(log_lik[lower:upper,], 
          file = paste("output/loglik_full_null_null_", i, ".rds", sep=""))
}