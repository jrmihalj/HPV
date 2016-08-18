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

  for (j in 1:d) {
    zero_vec[j] = 0;
    for (i in 1:n) {
      if (y[i, j] == 1) {
        sign[i, j] = 1;
      } else {
        sign[i, j] = -1;
      }
    }
  }
}

parameters {
  cholesky_factor_corr[d] L_Rho_patient;
  cholesky_factor_corr[d] L_Rho_visit;
  row_vector<lower=0>[d] abs_ystar[n];
  matrix[n_patient,d] e_patient;
  simplex[2] var_mat[d];
  matrix[d, d] betas_phi;
  matrix[d, d] betas_gam_R;
  vector[d] alphas;
}

transformed parameters {
  cholesky_factor_cov[d] L_Sigma_visit;
  cholesky_factor_cov[d] L_Sigma_patient;
  vector<lower = 0>[d] sd_visit;
  vector<lower = 0>[d] sd_patient;
  matrix[n, d] fixedef_phi;
  matrix[n, d] fixedef_gam;
  matrix[n, d] fixedef;
  matrix[d,d] betas_gam;
  matrix[n, d] mu_all;

  for (i in 1:d) {
    sd_visit[i] = sqrt(var_mat[i, 1]);
    sd_patient[i] = sqrt(var_mat[i, 2]);
  }

  L_Sigma_patient = diag_pre_multiply(sd_patient, L_Rho_patient);
  L_Sigma_visit = diag_pre_multiply(sd_visit, L_Rho_visit);

  for(j in 1:d){
    for(i in 1:d){
      if(i == j){
        betas_gam[i,j] = 0;
      }else{
        betas_gam[i,j] = betas_gam_R[i, j];
      }
    }
  }

  for (i in 1:n) {
    fixedef_phi[i,] = y[i,] * betas_phi;
    fixedef_gam[i,] = y[i,] * betas_gam;
  }

  for (j in 1:d) {
    for (i in 1:n) {
      if (visit[i] == 1){
        fixedef[i, j] = alphas[j];
      } else {
        if( y[i - n_patient, j] == 1){
          fixedef[i, j] = alphas[j] + fixedef_phi[i - n_patient, j];
        } else {
          fixedef[i, j] = alphas[j] + fixedef_gam[i - n_patient, j];
        }
      }
    }
  }

  mu_all = fixedef + e_patient[patient, ];
}

model {
  for (i in 1:d) var_mat[i] ~ dirichlet(dir_prior);
  L_Rho_patient ~ lkj_corr_cholesky(eta);
  L_Rho_visit ~ lkj_corr_cholesky(eta);
  to_vector(betas_phi) ~ normal(0, 1);
  to_vector(betas_gam_R) ~ normal(0, 1);
  alphas ~ normal(-2, 1.5);

  for (i in 1:n_patient){
    e_patient[i, ] ~ multi_normal_cholesky(zero_vec, L_Sigma_patient);
  }
  for (i in 1:n) {
    sign[i, ] .* abs_ystar[i, ] ~ multi_normal_cholesky(mu_all[i, ], L_Sigma_visit);
  }
}

generated quantities {
  matrix[d, d] Rho_patient;
  matrix[d, d] Rho_visit;

  Rho_patient = multiply_lower_tri_self_transpose(L_Rho_patient);
  Rho_visit = multiply_lower_tri_self_transpose(L_Rho_visit);
}
