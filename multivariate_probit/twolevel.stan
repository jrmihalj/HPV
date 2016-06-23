data {
  int n;
  int d;
  matrix[n,d] y;
  real<lower = 1> eta;
  int n_patient;
  int<lower = 0, upper = n_patient> patient[n];
  vector[2] dir_prior;
  int n_visit;
  int<lower = 1, upper = n_visit> visit[n];
  
}

transformed data {
  matrix[n, d] sign;
  
  for (i in 1:n) {
    for (j in 1:d) {
      if (y[i, j] == 1) {
        sign[i, j] <- 1;
      } else {
        sign[i, j] <- -1;
      }
    }
  }
}

parameters {
  cholesky_factor_corr[d] L_Rho_patient;
  cholesky_factor_corr[d] L_Rho_visit;
  matrix<lower = 0, upper=9>[n, d] abs_ystar;
  matrix[n_patient, d] e_patientR;
  simplex[2] var_mat[d];
  matrix[d,d] betas;
  //real<lower = 1, upper=5> eta;
  //vector<lower=0>[d] y_star_means;
  //real<lower=0.0001> sd_beta;
}

transformed parameters {
  cholesky_factor_cov[d] L_Sigma_visit;
  matrix[n_patient, d] e_patient;
  vector<lower = 0>[d] sd_visit;
  vector<lower = 0>[d] sd_patient;
  matrix[n, d] fixedef_R;
  matrix[n, d] fixedef;

  for (i in 1:d) {
    sd_visit[i] <- sqrt(var_mat[i, 1]);
    sd_patient[i] <- sqrt(var_mat[i, 2]);
  }
  e_patient <- (diag_pre_multiply(sd_patient, L_Rho_patient) * e_patientR')';
  L_Sigma_visit <- diag_pre_multiply(sd_visit, L_Rho_visit);
  
  fixedef_R <- y * betas;
  
  for (i in 1:n) {
    for (j in 1:d) {
      if (visit[i] == 1){
        fixedef[i, j] <- 0;
      } else {
        fixedef[i, j] <- fixedef_R[i - n_patient, j];
      }
    }
  }
  
}

model {
  for (i in 1:d) {
    var_mat[i] ~ dirichlet(dir_prior);
  }
  
  //sd_beta ~ cauchy(0,3);
  //eta ~ normal(0,100);
  //y_star_means ~ normal(0,10);
  
  L_Rho_patient ~ lkj_corr_cholesky(eta);
  L_Rho_visit ~ lkj_corr_cholesky(eta);
  to_vector(e_patientR) ~ normal(0, 1);
  to_vector(betas) ~ normal(0, 5);
  //Might make sense to make this heirarchical
  //i.e. each strain has an average abs_star
  //to_vector(abs_ystar) ~ normal(0, 2);
  // for (i in 1:d){
  //   abs_ystar[,d] ~ normal(y_star_means[d], 2);
  // }
  
  {
    
    matrix[n, d] occur;
    matrix[n, d] mu_all;
    
    //For efficiency:
    mu_all <- fixedef + e_patient[patient, ];
  
    //For efficiency:
    for (i in 1:n) {
      occur[i, ] <- sign[i, ] .* abs_ystar[i, ];
    }
  
    for (i in 1:n) {
      occur[i, ] ~ multi_normal_cholesky(mu_all[i, ], L_Sigma_visit);
    }
    
  }
  
  
  
}

generated quantities {
  matrix[d, d] Rho_patient;
  matrix[d, d] Rho_visit;

  Rho_patient <- multiply_lower_tri_self_transpose(L_Rho_patient);
  Rho_visit <- multiply_lower_tri_self_transpose(L_Rho_visit);
}
