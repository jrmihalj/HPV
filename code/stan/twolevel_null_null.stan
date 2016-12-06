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
  vector[n_strains] betas_tbv_phi;
  vector[n_strains] betas_tbv_gam;
  vector[n_strains] alphas;
  
  //Hyper-parameters
  real tbv_phi_mean;
  real tbv_gam_mean;
  real alpha_mean;
  
  real < lower = 0 > tbv_phi_sd; 
  real < lower = 0 > tbv_gam_sd;
  real < lower = 0 > alpha_sd; 
}

transformed parameters {
  matrix[n, n_strains] fixedef;
  matrix[n, n_strains] mu_all;


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

  mu_all = fixedef;
}

model {
  for (i in 1:n_strains) {
    abs_ystar[i] ~ normal(0, 1);
  }
  
  //Hyper-parameters
  tbv_phi_mean ~ normal(0, 5);
  tbv_gam_mean ~ normal(0, 5);
  alpha_mean ~ normal(0, 5);
  
  tbv_phi_sd ~ cauchy(0,2);
  tbv_gam_sd ~ cauchy(0,2);
  alpha_sd ~ cauchy(0,2);
  
  //Parameters
  betas_tbv_phi ~ normal(tbv_phi_mean, tbv_phi_sd);
  betas_tbv_gam ~ normal(tbv_gam_mean, tbv_gam_sd);
  alphas ~ normal(alpha_mean, alpha_sd);
  
  for(j in 1:n_strains){
    for (i in 1:n) {
      sign[i, j] .* abs_ystar[i, j] ~ normal(mu_all[i, j], 1);
    }
  }
  
}

generated quantities {
  row_vector[n_strains] log_lik[n];
  for(j in 1:n_strains){
  	for (i in 1:n){
    	log_lik[i,j] = normal_lpdf(sign[i, j] .* abs_ystar[i, j] | mu_all[i, j], 1);
	  }
  }
}
