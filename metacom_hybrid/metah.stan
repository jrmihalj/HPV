data {
  int n; // number of observations
  int n_unit; // number of sample units
  int n_time; // number of times sampled
  int<lower=1, upper=n_unit> unit[n];
  int<lower=1, upper=n_time> time[n];
  int m; // number of species
  matrix[n, m] y_mat;
  int<lower = 0, upper = 1> y[n, m]; // presence/absence observations
}

parameters {
  cholesky_factor_corr[m] L_R_o;
  cholesky_factor_corr[m] L_R_p;
  matrix[m, m] beta_phi;
  matrix[m, m] beta_gam;
  matrix[n_unit, m] e_raw_p;
  matrix[n, m] e_raw_o;
}

transformed parameters {
  matrix[n, m] mu_phi;
  matrix[n, m] mu_gam;
  vector[n] sum_phi;
  vector[n] sum_gam;
  matrix[n, m] z;
  matrix[n, m] e_o;
  matrix[n_unit, m] e_p;

  e_o <- e_raw_o * L_R_o;
  e_p <- e_raw_p * L_R_p;
  
  mu_phi <- y_mat * beta_phi;
  mu_gam <- y_mat * beta_gam;
  
  for(i in 1:n){
    sum_phi[i] <- sum(mu_phi[i,]);
    sum_gam[i] <- sum(mu_phi[i,]);
  }

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
    
    //Calculate the z-vals
    for (j in 1:m){
      for (i in 1:n){
        if (time[i] == 1) {
            z[i, j] <- e_all[i, j];
        } else {//Current visit depends on previous visit
            z[i, j] <- (y_mat[i - n_unit, j]) * sum_phi[i] + //persistence
            		       (1 - y_mat[i - n_unit, j]) * sum_gam[i] + //colonization
            		       e_all[i, j]; //correlated and nested ranef
        }
      }
    }
    
  }//End temporary scope
  
  
  
}

model {
  L_R_o ~ lkj_corr_cholesky(2);
  L_R_p ~ lkj_corr_cholesky(2);
  to_vector(e_raw_p) ~ normal(0, 1);
  to_vector(e_raw_o) ~ normal(0, 1);

  for(j in 1:m){
    for(i in 1:m){
      beta_phi[i,j] ~ double_exponential(0,1);
      beta_gam[i,j] ~ double_exponential(0,1);
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
