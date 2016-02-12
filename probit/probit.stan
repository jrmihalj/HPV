data {
  int n; // number of observations
  int n_unit; // number of sample units
  int<lower=1, upper=n_unit> unit[n];
  int m; // number of species
  int k; // number of columns in design matrix
  matrix[n, k] X; // design matrix
  matrix[n, m] X_intxn; // interaction design matrix
  int<lower = 0, upper = 1> y[n, m]; // presence/absence observations
}

parameters {
  cholesky_factor_corr[m] L_R_o;
  cholesky_factor_corr[m] L_R_p;
  matrix[k, m] beta;
  matrix[n_unit, m] e_raw_p;
  matrix[n, m] e_raw_o;
  matrix[m, m] beta_intxn;
}

transformed parameters {
  matrix[n, m] mu;
  matrix[n, m] z;
  matrix[n, m] e_o;
  matrix[n_unit, m] e_p;

  e_o <- e_raw_o * L_R_o;
  e_p <- e_raw_p * L_R_p;
  
  mu <- X * beta + X_intxn * beta_intxn;
  {
    // temporary scope for logit of cdf of e
    matrix[n, m] e_tmp;
    matrix[n, m] e_tmp2;
    matrix[n, m] e_all;
    for (j in 1:m){
      for (i in 1:n){
        e_tmp[i, j] <- e_p[unit[i], j] + e_o[i, j];
      }
    }
    
    //Normalize
    e_tmp2 <- (e_tmp - mean(e_tmp))/sd(e_tmp);
    
    for (j in 1:m){
      for (i in 1:n){
        e_all[i, j] <- logit(Phi( e_tmp2[i, j] ));
      }
    }
    
    z <- mu + e_all;
  }
}

model {
  L_R_o ~ lkj_corr_cholesky(2);
  L_R_p ~ lkj_corr_cholesky(2);
  to_vector(e_raw_p) ~ normal(0, 1);
  to_vector(e_raw_o) ~ normal(0, 1);
  to_vector(beta) ~ normal(0, 1);
  
  for (j in 1:m) {
    for (i in 1:m) {
      beta_intxn[i, j] ~ double_exponential(0, 1);
    }
  }
  
  for (j in 1:m){
    for (i in 1:n){
      y[i, j] ~ bernoulli(Phi(z[i, j]));
    }
  }
}

generated quantities {
  matrix[m, m] R_p;
  matrix[m, m] R_o;

  R_o <- multiply_lower_tri_self_transpose(L_R_o);
  R_p <- multiply_lower_tri_self_transpose(L_R_p);

}
