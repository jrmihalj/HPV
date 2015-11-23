data 
{
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
	
parameters 
{
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
	cholesky_factor_corr[n_strains] L_patient_phi1_phi2; // 
	matrix[n_strains , n_patients] z_patient_phi1_phi2;
	vector<lower=0> [n_strains] tau_patient_phi1_phi2; // across patient std deviations
	
	cholesky_factor_corr[n_strains] L_patient_phi1_gamma2; // 
	matrix[n_strains , n_patients] z_patient_phi1_gamma2;
	//vector<lower=0> [n_strains] tau_patient_phi1_gamma2; // across patient std deviations
	
	cholesky_factor_corr[n_strains] L_patient_gamma1_phi2; // 
	matrix[n_strains , n_patients] z_patient_gamma1_phi2;
	//vector<lower=0> [n_strains] tau_patient_gamma1_phi2; // across patient std deviations
	
	cholesky_factor_corr[n_strains] L_patient_gamma1_gamma2; // 
	matrix[n_strains , n_patients] z_patient_gamma1_gamma2;
	vector<lower=0> [n_strains] tau_patient_gamma1_gamma2; // across patient std deviations
}

transformed parameters
{
	//persistence and colonization probabilities 
	vector[n_obs] phi; //time,patient,strain specific colonization probability
	vector[n_obs] gam; //time,patient,strain specific persistence probability
	
	vector<lower=0> [n_strains] tau_patient_phi1_gamma2;
  	vector<lower=0> [n_strains] tau_patient_gamma1_phi2;
	
	// correlations 
	matrix[n_patients, n_strains] alpha_patient_phi1_phi2;
	matrix[n_patients, n_strains] alpha_patient_gamma1_phi2;
	matrix[n_patients, n_strains] alpha_patient_phi1_gamma2;
	matrix[n_patients, n_strains] alpha_patient_gamma1_gamma2;
	
	// occupancy probability 
	vector<lower=0, upper=1>[n_obs] psi_col;
	vector<lower=0, upper=1>[n_obs] psi_per;
	
	
	//between-strain effects
	real added_effect_phi; // addition of between-strain phi/gam random effect 
	real added_effect_gam; // addition of between-strain phi/gam random effect 
	
	//fixed covariate effects 
	//vector[n_strains] beta_pat; // fixed across-patients covariate effect
  	//vector[n_strains] beta_time; // fixed within patient covariate effect 
  	//beta_pat <- beta_pat_mean + beta_pat_sd * beta_patR;
  	//beta_time <- beta_time_mean + beta_time_sd * beta_timeR;
	
	tau_patient_phi1_gamma2[1] <- tau_patient_phi1_phi2[1];
	tau_patient_phi1_gamma2[2] <- tau_patient_gamma1_gamma2[2];
	tau_patient_gamma1_phi2[1] <- tau_patient_gamma1_gamma2[1];
	tau_patient_gamma1_phi2[2] <- tau_patient_phi1_phi2[2];
	
	//random effects
	alpha_patient_phi1_phi2 <- (diag_pre_multiply(tau_patient_phi1_phi2, L_patient_phi1_phi2) * z_patient_phi1_phi2)';
	alpha_patient_gamma1_gamma2 <- (diag_pre_multiply(tau_patient_gamma1_gamma2, L_patient_gamma1_gamma2) * z_patient_gamma1_gamma2)';
	alpha_patient_phi1_gamma2 <- (diag_pre_multiply(tau_patient_phi1_gamma2, L_patient_phi1_gamma2) * z_patient_phi1_gamma2)';
	alpha_patient_gamma1_phi2  <- (diag_pre_multiply(tau_patient_gamma1_phi2, L_patient_gamma1_phi2) * z_patient_gamma1_phi2)';
	
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
  	  	psi_col[i] <- inv_logit(psi_mean); // p(colonization)
  	  	psi_per[i] <- inv_logit(psi_mean); // p(persistence)
  	  } // end if
  	  else {
  	 	 psi_per[i] <- phi[i-n_strains] * psi[i-n_strains] ;
  	 	 psi_col[i] <- gam[i-n_strains] * (1 - psi[i-n_strains]);  
  	  }// end else
	}// end n_obs
}

model {

  vector<lower=0, upper=1>[n_obs] psi;
  
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
  for(i in 1:n_obs){
  	if( i == 1){
  		psi[i] <- inv_logit(psi_mean);
  	}
  	if( i > 1){
  		if( Y[i - n_strains] == 1){
  			psi[i] <- psi_per[i];
  		}// end if 
  		if(Y[i - n_strains] == 0) {
  			psi[i] <- psi_col[i];
  		}// end if
  	}// end if 
  }// end for i in 1:n_obs
  
  Y ~ bernoulli(psi);
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