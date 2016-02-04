data 
{
]	int<lower=0> n_strains;
	int<lower=0> n_obs;
	int<lower=0> n_patients;
	int<lower=0> n_visits_total; //total number of clinic visits across patients
	int<lower=0, upper=1> Y[n_obs]; //occupancy observations
	int<lower=1, upper=n_strains> strain[n_obs]; // strain for this observation
	int<lower=1, upper=n_patients> patient[n_obs]; // patient for this observation/strain 
	int<lower=1> visit_pat[n_obs]; // visit for this patient for this observation/strain
	int<lower=1> Visit[n_obs]; // observation (clinic visits across patients)

}
	
parameters 
{
  // intercepts
  real psi_initial; // initial occupancy probability
  
  // params of fixed effects 
  vector[n_strains] betas[n_strians];

  // random  patient effects
	cholesky_factor_corr[n_strains] L_patient; // 
	matrix[n_strains , n_patients] z_patient;
	vector<lower=0> [n_strains] tau_patient; // across patient std deviations
	

	
}

transformed parameters
{
	
	// correlations 
	matrix[n_patients, n_strains] alpha_patient;
	
	// occupancy probability 
	vector<lower=0, upper=1>[n_obs] psi;
	
	//between-strain effects

	
	//random effects
	alpha_patient <- (diag_pre_multiply(tau_patient, L_patient) * z_patient)';
	
	for(i in 1:n_obs){
	
		if(strain[i] == 1){
			added_effect_phi <- alpha_patient_phi1_gamma2[patient[i], strain[i]];
			added_effect_gam <- alpha_patient_gamma1_phi2[patient[i], strain[i]];
		}
		
		if(strain[i] == 2){
			added_effect_phi <- alpha_patient_gamma1_phi2[patient[i], strain[i]];
			added_effect_gam <- alpha_patient_phi1_gamma2[patient[i], strain[i]];
		}
		
		phi[i] <- inv_logit(phi_mean 
		                      + alpha_patient_phi1_phi2[patient[i], strain[i]]
		                      + added_effect_phi);
		                      
		                    //+ beta_pat[strain[i]]*X_pat[i]
		                   
		gam[i] <- inv_logit(gam_mean 
		                      + alpha_patient_gamma1_gamma2[patient[i], strain[i]]
		                      + added_effect_gam);
		                       		                       
		                    //+ beta_pat[strain[i]]*X_pat[i] 
		                   
		                   
	  if(visit_pat[i] == 1){
  	  	psi[i] <- inv_logit(psi_mean);
  	  } // end if
  	  else {
  	 	 psi[i] <- phi[i-n_strains] * psi[i-n_strains] 
  	               + gam[i-n_strains] * (1 - psi[i-n_strains]);  
  	  }// end else
	}// end n_obs
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
  L_patient_phi1_phi2 ~ lkj_corr_cholesky(2);
  tau_patient_phi1_phi2 ~ normal(0, 2);
  to_vector(z_patient_phi1_phi2) ~ normal(0, 1);
  
  L_patient_gamma1_gamma2 ~ lkj_corr_cholesky(2);
  tau_patient_gamma1_gamma2 ~ normal(0, 2);
  to_vector(z_patient_gamma1_gamma2) ~ normal(0, 1);
  
  L_patient_phi1_gamma2 ~ lkj_corr_cholesky(2);
  //constrain taus
 // tau_patient_phi1_gamma2[1] <- tau_patient_phi1_phi2[1];
 // tau_patient_phi1_gamma2[2] <- tau_patient_gamma1_gamma2[2];
  to_vector(z_patient_phi1_gamma2) ~ normal(0, 1);
  
  L_patient_gamma1_phi2 ~ lkj_corr_cholesky(2);
  //tau_patient_gamma1_phi2[1] <- tau_patient_gamma1_gamma2[1];
 // tau_patient_gamma1_phi2[2] <- tau_patient_phi1_phi2[2];
  to_vector(z_patient_gamma1_phi2) ~ normal(0, 1);
  
  // likelihood of occupancy
  
  for( i in 1:nObs){
  	if( i == 1) {
  		Y[i] ~ bernoulli(psi);
  	}
  	else {
  	
  
  	}
  
  }
} // end model block

generated quantities {
  matrix[n_strains , n_strains ] cor_patient_phi1_phi2;
  matrix[n_strains , n_strains ] cor_patient_phi1_gamma2;
  matrix[n_strains , n_strains ] cor_patient_gamma1_phi2;
  matrix[n_strains , n_strains ] cor_patient_gamma1_gamma2;
  
  // create the correlation matrix to draw random effects
  cor_patient_phi1_phi2 <- multiply_lower_tri_self_transpose(L_patient_phi1_phi2);
  cor_patient_phi1_gamma2 <- multiply_lower_tri_self_transpose(L_patient_phi1_gamma2);
  cor_patient_gamma1_phi2 <- multiply_lower_tri_self_transpose(L_patient_gamma1_phi2);
  cor_patient_gamma1_gamma2 <- multiply_lower_tri_self_transpose(L_patient_gamma1_gamma2);
}  