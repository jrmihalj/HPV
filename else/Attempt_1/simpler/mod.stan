data {
	int<lower=0> n_strains;
	int<lower=0> n_patients;
	int<lower=0> n_vis; //total number of clinic visits across patients
	int<lower=0, upper=1> z[n_patients, n_vis, n_strains]; //occupancy observations
}

parameters {
  // random effects
	cholesky_factor_corr[n_strains * 2] L_patient; // * 2 b/c inc. phi & gamma
	matrix[n_strains * 2, n_patients] z_patient;
	vector<lower=0> [n_strains * 2] tau_patient; // across patient std deviations
	matrix<lower=0, upper=1>[n_patients, n_strains] psi0;
	real mu_psi0;
	real<lower=0> sd_psi0;
}

transformed parameters{
	real<lower=0, upper=1> phi[n_patients, n_vis, n_strains]; // colonization
	real<lower=0, upper=1> gam[n_patients, n_vis, n_strains]; // persistence
	matrix[n_patients, n_strains * 2] alpha_patient;
	real<lower=0, upper=1> psi[n_patients, n_vis, n_strains];
	
	//random effects
	alpha_patient <- (diag_pre_multiply(tau_patient, L_patient) * z_patient)';

	for (i in 1:n_patients) { 
    for (j in 1:n_vis){
      for (k in 1:n_strains){
        phi[i, j, k] <- inv_logit(alpha_patient[i, k]);
        gam[i, j, k] <- inv_logit(alpha_patient[i, k + n_strains]);
        if (j == 1) {
          psi[i, j, k] <- inv_logit(psi0[i, k]);
        } else {
          psi[i, j, k] <- z[i, j - 1, k] * phi[i, j - 1, k]
            + (1 - z[i, j - 1, k]) * gam[i, j - 1, k];
        }
      }
    }
  }
}

model {
  // prior on across patient random effects
  L_patient ~ lkj_corr_cholesky(5);
  tau_patient ~ cauchy(0, 2);
  to_vector(z_patient) ~ normal(0, 1);
  mu_psi0 ~ normal(0, 1);
  sd_psi0 ~ normal(0, 1);
  to_vector(psi0) ~ normal(mu_psi0, sd_psi0);
  
  for (i in 1:n_patients){
    for (j in 1:n_vis){
      for (k in 1:n_strains){
        z[i, j, k] ~ bernoulli(psi[i, j, k]);
      }
    }
  }
} // end model block

generated quantities {
  matrix[n_strains * 2, n_strains * 2] cor_patient;
  cor_patient <- multiply_lower_tri_self_transpose(L_patient);
}