data {
  int n; // number of sample units
  int m; // number of species
  int k; // number of columns in design matrix
  matrix[n, k] X; // design matrix
  matrix[n, m] X_intxn; // interaction design matrix
  int<lower = 0, upper = 1> y[n, m]; // presence/absence observations
}

parameters {
  cholesky_factor_corr[m] L_R;
  matrix[k, m] beta;
  matrix[n, m] e_raw;
  matrix[m, m] beta_intxn;
}

transformed parameters {
  matrix[n, m] mu;
  matrix[n, m] z;
  matrix[n, m] e;
  
  e <- e_raw * L_R;
  
  mu <- X * beta + X_intxn * beta_intxn;
  {
    // temporary scope for logit of cdf of e
    matrix[n, m] e_tmp;
    for (i in 1:n){
      for (j in 1:m){
        e_tmp[i, j] <- logit(Phi(e[i, j]));
      }
    }
    z <- mu + e_tmp;
  }
}

model {
  L_R ~ lkj_corr_cholesky(2);
  to_vector(e_raw) ~ normal(0, 1);
  to_vector(beta) ~ normal(0, 1);
  
  for (i in 1:m) {
    for (j in 1:m) {
      beta_intxn[i, j] ~ double_exponential(0, 1);
    }
  }
  
  for (i in 1:n){
    for (j in 1:m){
      y[i, j] ~ bernoulli(Phi(z[i, j]));
    }
  }
}

generated quantities {
  matrix[m, m] R;
  R <- multiply_lower_tri_self_transpose(L_R);
}
