data {
  int n;
  int n_strains;
  row_vector<lower = 0, upper = 1>[n_strains] y[n];
  int n_patient;
  int n_visit_max;
  int<lower = 1, upper = n_patient> patient[n];
  int<lower = 1, upper = n_visit_max> visit[n];
  vector[n] tbv;
  vector[2] dir_prior;
  real<lower = 1> eta;
}

transformed data {
  row_vector[n_strains] sign[n];
  row_vector[n_strains] zero_vec;

  for (j in 1:n_strains) {
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
  cholesky_factor_corr[n_strains] L_Rho_patient;
  cholesky_factor_corr[n_strains] L_Rho_visit;
  row_vector<lower=0>[n_strains] abs_ystar[n];
  matrix[n_patient,n_strains] e_patient;
  simplex[2] var_mat[n_strains];
  vector[n_strains] betas_tbv_phi;
  vector[n_strains] betas_tbv_gam;
  vector[n_strains] alphas;
}

transformed parameters {
  cholesky_factor_cov[n_strains] L_Sigma_visit;
  cholesky_factor_cov[n_strains] L_Sigma_patient;
  vector<lower = 0>[n_strains] sd_visit;
  vector<lower = 0>[n_strains] sd_patient;
  matrix[n, n_strains] fixedef;
  matrix[n, n_strains] mu_all;

  for (i in 1:n_strains) {
    sd_visit[i] = sqrt(var_mat[i, 1]);
    sd_patient[i] = sqrt(var_mat[i, 2]);
  }

  L_Sigma_patient = diag_pre_multiply(sd_patient, L_Rho_patient);
  L_Sigma_visit = diag_pre_multiply(sd_visit, L_Rho_visit);


  for (j in 1:n_strains) {
    for (i in 1:n) {
      if (visit[i] == 1){
        fixedef[i, j] = alphas[j];
      } else {
        if( y[i - 1, j] == 1){
          fixedef[i, j] = alphas[j] + tbv[i] * betas_tbv_phi[j];
        } else {
          fixedef[i, j] = alphas[j] + tbv[i] * betas_tbv_gam[j];
        }
      }
    }
  }

  mu_all = fixedef + e_patient[patient, ];
}

model {
  for (i in 1:n_strains) {
    var_mat[i] ~ dirichlet(dir_prior);
    abs_ystar[i] ~ normal(0, 1);
  }
  
  L_Rho_patient ~ lkj_corr_cholesky(eta);
  L_Rho_visit ~ lkj_corr_cholesky(eta);
  betas_tbv_phi ~ normal(0, 1);
  betas_tbv_gam ~ normal(0, 1);
  alphas ~ normal(-2, 1.5);
  
  for (i in 1:n_patient){
    e_patient[i, ] ~ multi_normal_cholesky(zero_vec, L_Sigma_patient);
  }
  for (i in 1:n) {
    sign[i, ] .* abs_ystar[i, ] ~ multi_normal_cholesky(mu_all[i, ], L_Sigma_visit);
  }
}

generated quantities {
  matrix[n_strains, n_strains] Rho_patient;
  matrix[n_strains, n_strains] Rho_visit;

  Rho_patient = multiply_lower_tri_self_transpose(L_Rho_patient);
  Rho_visit = multiply_lower_tri_self_transpose(L_Rho_visit);
}