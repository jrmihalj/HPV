data {
  int n;
  int d;
  row_vector<lower = 0, upper = 1>[d] y[n];
  real<lower = 1> eta;
  int n_patient;
  int<lower = 1, upper = n_patient> patient[n];
  vector[2] dir_prior;
  int n_visit;
  int<lower = 1, upper = n_visit> visit[n];
  
}

transformed data {

  row_vector[d] sign[n];
  row_vector[d] zero_vec;
  
  for (i in 1:d) {
    zero_vec[i] <- 0;
  }

  for (j in 1:d) {
    for (i in 1:n) {
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
  row_vector<lower=0>[d] abs_ystar[n];
  row_vector[d] e_patient[n_patient];
  simplex[2] var_mat[d];
  matrix[d,d] betas;
  //real<lower = 1, upper=5> eta;
  //real<lower=0.0001> sd_beta;
}

transformed parameters {
  cholesky_factor_cov[d] L_Sigma_visit;
  cholesky_factor_cov[d] L_Sigma_patient;
  matrix[n_patient, d] e_patient2;
  vector<lower = 0>[d] sd_visit;
  vector<lower = 0>[d] sd_patient;
  matrix[n, d] fixedef_R;
  matrix[n, d] fixedef;
  matrix[n, d] mu_all;
  row_vector[d] occur[n];

    
  for (i in 1:d) {
    sd_visit[i] <- sqrt(var_mat[i, 1]);
    sd_patient[i] <- sqrt(var_mat[i, 2]);
  }
  
  //For efficiency, instead of this:
  //e_patient <- (diag_pre_multiply(sd_patient, L_Rho_patient) * e_patientR')';
  //Do below, and then have e_patient be distributed as a multi_normal_cholesky with zero mean
  L_Sigma_patient <- diag_pre_multiply(sd_patient, L_Rho_patient);
  L_Sigma_visit <- diag_pre_multiply(sd_visit, L_Rho_visit);
  
  for (i in 1:n) {
    fixedef_R[i,] <- y[i,] * betas;
  }
  
  //This is simply converting e_patient into a matrix for vectorization below
  for (j in 1:d) {
    for (i in 1:n_patient) {
      e_patient2[i, j] <- e_patient[i, j];
    }
  }
  
  for (j in 1:d) {
    for (i in 1:n) {
      if (visit[i] == 1){
        fixedef[i, j] <- 0;
      } else {
        fixedef[i, j] <- fixedef_R[i - n_patient, j];
      }
    }
  }
  
  //For efficiency:
  mu_all <- fixedef + e_patient2[patient, ];
  
  //For efficiency:
  for (i in 1:n) {
    occur[i, ] <- sign[i, ] .* abs_ystar[i, ];
  }

  
}

model {
  for (i in 1:d) {
    var_mat[i] ~ dirichlet(dir_prior);
  }
  
  //sd_beta ~ cauchy(0,3);
  //eta ~ normal(0,100);
  
  L_Rho_patient ~ lkj_corr_cholesky(eta);
  L_Rho_visit ~ lkj_corr_cholesky(eta);
  //to_vector(e_patientR) ~ normal(0, 1);
  to_vector(betas) ~ normal(0, 3);
  
  for (i in 1:n_patient){
    //random patient-level effects:
    e_patient[i, ] ~ multi_normal_cholesky(zero_vec, L_Sigma_patient);
  }

  for (i in 1:n) {
    occur[i, ] ~ multi_normal_cholesky(mu_all[i, ], L_Sigma_visit);
  }

  
}

generated quantities {
  matrix[d, d] Rho_patient;
  matrix[d, d] Rho_visit;

  Rho_patient <- multiply_lower_tri_self_transpose(L_Rho_patient);
  Rho_visit <- multiply_lower_tri_self_transpose(L_Rho_visit);
}
