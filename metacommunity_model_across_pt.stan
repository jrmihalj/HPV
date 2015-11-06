data {
	int<lower=0> n_strains;
	int<lower=0> n_obs;
	int<lower=0> n_patients;
	int<lower=0> n_visits_total; //total number of clinic visits across patients
	int<lower=0, upper=1> Y[n_obs]; //occupancy observations
	int<lower=1, upper=n_strains> strain[n_obs]; // strain for this observation
	int<lower=1, upper=n_patients> patient[n_obs]; // patient for this observation/strain 
	int<lower=1> visit_pat[n_obs]; // visit for this patient for this observation/strain
	int<lower=1> Visit[n_obs]; // observation (clinic visits across patients)
	vector[n_obs] X_time; // time varying covariates
	vector[n_obs] X_pat; // time invariant covariates
}
	
parameters {
  // intercepts
  real phi_mean ; // global mean of colonization probability
  real gam_mean; // global mean of vacancy probability 
  real psi_mean; // initial occupancy probability
  
  // params of fixed effects 
  //real beta_pat_mean; 
  //real beta_time_mean;
  //real<lower=0> beta_pat_sd;
  //real<lower=0> beta_time_sd;
  //vector[n_strains] beta_patR; // fixed across-patients covariate effect
  //vector[n_strains] beta_timeR; // fixed within patient covariate effect 

  // random effects
	cholesky_factor_corr[n_strains] L_patient; // 
	matrix[n_strains , n_patients] z_patient;
	vector<lower=0> [n_strains] tau_patient; // across patient std deviations

}

transformed parameters{
	vector[n_obs] phi; //time,patient,strain specific colonization probability
	vector[n_obs] gam; //time,patient,strain specific persistence probability
	matrix[n_patients, n_strains] alpha_patient;
	vector<lower=0, upper=1>[n_obs] psi;
	
	
	//vector[n_strains] beta_pat; // fixed across-patients covariate effect
  //vector[n_strains] beta_time; // fixed within patient covariate effect 
  
  //beta_pat <- beta_pat_mean + beta_pat_sd * beta_patR;
  //beta_time <- beta_time_mean + beta_time_sd * beta_timeR;
	
	//random effects
	alpha_patient <- (diag_pre_multiply(tau_patient, L_patient) * z_patient)';

	
	for(i in 1:n_obs){
		phi[i] <- inv_logit(phi_mean 
		                      + alpha_patient[patient[i], strain[i]] 
		                    //+ beta_pat[strain[i]]*X_pat[i]
		                   
		gam[i] <- inv_logit(gam_mean 
		                      + alpha_patient[patient[i], n_strains + strain[i]]
		                    //+ beta_pat[strain[i]]*X_pat[i] 
		                   
		                   
	  if(visit_pat[i] == 1){
  	  psi[i] <- inv_logit(psi_mean);
  	} else {
  	  psi[i] <- phi[i-n_strains] * psi[i-n_strains] 
  	              + gam[i-n_strains] * (1 - psi[i-n_strains]);  
  	  }
	}
}

model {
  phi_mean ~ normal(0, .5);
  gam_mean ~ normal(0, .5);
  psi_mean ~ normal(0, .5);
  
  // prior on fixed effects
  //beta_pat_mean ~ normal(0, 1.5);
  //beta_time_mean ~ normal(0, 1.5);
  //beta_pat_sd ~ normal(0, 2);
  //beta_time_sd ~ normal(0, 2);

  // coefficients (matt trick)
  //beta_patR ~ normal(0, 1); // normal(beta_pat_mean, beta_pat_sd)
  //beta_timeR ~ normal(0, 1);  // normal(beta_time_mean, beta_time_sd)

  // prior on across patient random effects
  L_patient ~ lkj_corr_cholesky(2);
  tau_patient ~ normal(0, 2);
  to_vector(z_patient) ~ normal(0, 1);


  Y ~ bernoulli(psi);
} // end model block

generated quantities {
  matrix[n_strains , n_strains ] cor_patient;
  // create the correlation matrix to draw random effects
  cor_patient <- multiply_lower_tri_self_transpose(L_patient);
}  