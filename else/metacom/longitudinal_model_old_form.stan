data {
	int<lower=0> n_strains;
	int<lower=0> n_obs;
	int<lower=0> n_patients;
	int<lower=0, upper=1> Y[n_obs]; //occupancy observations
	int<lower=1, upper=n_strains> strain[n_obs]; // strain for this observation
	int<lower=1, upper=n_patients> patient[n_obs]; // patient for this observation/strain 
	//int<lower=1> Visit[n_obs]; // observation (clinic visits across patients)
	//vector[n_obs] X_time; // time varying covariates
	//vector[n_obs] X_pat; // time invariant covariates
}
	
parameters {
  // intercepts
  real phi_mean ; // global mean of colonization probability
  real gam_mean; // global mean of vacancy probability 
  real psi_mean; // initial occupancy probability
  
  // params of fixed effects : ignore for now 
  /*real beta_pat_mean; 
  real beta_time_mean;
  real<lower=0> beta_pat_sd;
  real<lower=0> beta_time_sd;
  vector[n_strains] beta_patR; // fixed across-patients covariate effect
  vector[n_strains] beta_timeR; // fixed within patient covariate effect */

  // random effects
  matrix[n_patients, n_strains * 2] alpha_raw_patient; // * 2 b/c inc. phi & gamma
  matrix[n_visits_total, n_strains * 2] alpha_raw_time;
  cholesky_factor_corr[n_strains * 2] L_patient; 
  cholesky_factor_corr[n_strains * 2] L_time;

}

transformed parameters{
	vector[n_obs] phi; //time,patient,strain specific colonization probability
	vector[n_obs] gam; //time,patient,strain specific persistence probability
	vector[n_obs] z_phi;
	vector[n_obs] z_gam;
	matrix[n_patients, n_strains * 2] alpha_patient;
	matrix[n_obs, n_strains * 2] alpha_time;
	vector<lower=0, upper=1>[n_obs] psi;
	
  	//vector[n_strains] beta_pat; // fixed across-patients covariate effect
 	//vector[n_strains] beta_time; // fixed within patient covariate effect 
  
  	//beta_pat <- beta_pat_mean + beta_pat_sd * beta_patR;
  	//beta_time <- beta_time_mean + beta_time_sd * beta_timeR;
	
	//random effects
	alpha_patient <- alpha_raw_patient * L_patient;
	alpha_time <- alpha_raw_time * L_time; 
	
	{// temporary scope for logit of cdf of error terms
    real[n_obs] error_tmp_phi;
    real[n_obs] error_tmp2_phi;
    real[n_obs] error_all_phi;
    
    real[n_obs] error_tmp_gam;
    real[n_obs] error_tmp2_gam;
    real[n_obs] error_all_gam;
 
    for (i in 1:n_obs){
    	e_tmp_phi[i] <- alpha_patient[patient[i], strain[i]] + alpha_time[i, strain[i]];
    	e_tmp_gam[i] <- alpha_patient[patient[i], strain[i] + n_strains] + alpha_time[i, strain[i]  + n_strains];
    }

    //Normalize
    e_tmp2_phi <-  (e_tmp_phi - mean(e_tmp_phi))/sd(e_tmp_phi);
    e_tmp2_gam <-  (e_tmp_gam - mean(e_tmp_gam))/sd(e_tmp_gam);
    
   
    for (i in 1:n_obs){
    	e_all_phi[i] <- logit(Phi(e_tmp2_phi[i] ));
    	e_all_gam[i] <- logit(Phi(e_tmp2_gam[i] ));
    }
    
    
    z_phi <- phi_mean + e_all_phi; // ignore fixed effects for now 
	z_gam <- gam_mean + e_all_gam; // ignore fixed effects for now 
	}
	
	phi <- Phi(z_phi);
	gam <- Phi(z_gam);

	for( i in 1:n_obs){
		if(visit_pat[i] == 1){
  	  		psi[i] <- inv_logit(psi_mean);
  		} 
  		else {
  	  		psi[i] <- phi[i-n_strains] * psi[i-n_strains] + gam[i-n_strains] * (1 - psi[i-n_strains]);
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
  // prior on within patient random effects
  L_time ~ lkj_corr_cholesky(2);


  Y ~ bernoulli(psi);
} // end model block

generated quantities {
  matrix[n_strains * 2, n_strains * 2] cor_patient;
  matrix[n_strains * 2, n_strains * 2] cor_time;
  
  cor_patient <- multiply_lower_tri_self_transpose(L_patient);
  cor_time <- multiply_lower_tri_self_transpose(L_time);
}