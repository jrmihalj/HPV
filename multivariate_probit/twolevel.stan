data {
  int n;
  int d;
  int<lower = 0, upper = 1> y[n, d];
  real<lower = 1> eta;
  int n_patient;
  int<lower = 0, upper = n_patient> patient[n];
  vector[2] dir_prior;
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
  matrix<lower = 0>[n, d] abs_ystar;
  matrix[n_patient, d] e_patientR;
  simplex[2] var_mat[d];
}

transformed parameters {
  cholesky_factor_cov[d] L_Sigma_visit;
  matrix[n_patient, d] e_patient;
  vector<lower = 0>[d] sd_visit;
  vector<lower = 0>[d] sd_patient;

  for (i in 1:d) {
    sd_visit[i] <- sqrt(var_mat[i, 1]);
    sd_patient[i] <- sqrt(var_mat[i, 2]);
  }
  e_patient <- (diag_pre_multiply(sd_patient, L_Rho_patient) * e_patientR')';
  L_Sigma_visit <- diag_pre_multiply(sd_visit, L_Rho_visit);
}

model {
  for (i in 1:d) {
    var_mat[i] ~ dirichlet(dir_prior);
  }

  L_Rho_patient ~ lkj_corr_cholesky(eta);
  L_Rho_visit ~ lkj_corr_cholesky(eta);
  to_vector(e_patientR) ~ normal(0, 1);
  for (i in 1:n) {
    sign[i, ] .* abs_ystar[i, ] ~ multi_normal_cholesky(e_patient[patient[i], ],
                                                          L_Sigma_visit);
  }
}

generated quantities {
  matrix[d, d] Rho_patient;
  matrix[d, d] Rho_visit;

  Rho_patient <- multiply_lower_tri_self_transpose(L_Rho_patient);
  Rho_visit <- multiply_lower_tri_self_transpose(L_Rho_visit);
}
