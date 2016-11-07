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
  row_vector<lower=0>[n_strains] abs_ystar[n];
  matrix[n_strains, n_strains] betas_phi;
  matrix[n_strains, n_strains] betas_gam_R;
  vector[n_strains] betas_tbv_phi;
  vector[n_strains] betas_tbv_gam;
  vector[n_strains] alphas;
}

transformed parameters {
  matrix[n, n_strains] fixedef_phi;
  matrix[n, n_strains] fixedef_gam;
  matrix[n, n_strains] fixedef;
  matrix[n_strains,n_strains] betas_gam;
  matrix[n, n_strains] mu_all;

  for(j in 1:n_strains){
    for(i in 1:n_strains){
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

  for (j in 1:n_strains) {
    for (i in 1:n) {
      if (visit[i] == 1){
        fixedef[i, j] = alphas[j];
      } else {
        if( y[i - 1, j] == 1){
          fixedef[i, j] = alphas[j] + fixedef_phi[i - 1, j] + tbv[i] * betas_tbv_phi[j];
        } else {
          fixedef[i, j] = alphas[j] + fixedef_gam[i - 1, j] + tbv[i] * betas_tbv_gam[j];
        }
      }
    }
  }

  mu_all = fixedef;
}

model {
  for (i in 1:n_strains) {
    abs_ystar[i] ~ normal(0, 1);
  }
  
  to_vector(betas_phi) ~ normal(0, 1);
  to_vector(betas_gam_R) ~ normal(0, 1);
  betas_tbv_phi ~ normal(0, 1);
  betas_tbv_gam ~ normal(0, 1);
  alphas ~ normal(-2, 1.5);
  
  for(j in 1:n_strains){
    for (i in 1:n) {
      sign[i, j] .* abs_ystar[i, j] ~ normal(mu_all[i, j], 1);
    }
  }
  
}
