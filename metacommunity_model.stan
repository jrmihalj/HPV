data
{
	int<lower=0> n_strains;
	int<lower=0> n_obs;
	int<lower=0> n_patients;
	int<lower=0> possible_visits; // store max possible visits for within-patient random effects
	
	int<lower=0> Y[n_obs]; //occupancy observations
	int<lower=0> strain[n_obs]; // strain for this observation
	int<lower=0> patient[n_obs]; // patient for this observation/strain 
	int<lower=0> visit_pat[n_obs]; // visit for this patient for this observation/strain
	
	real X_time[n_obs]; // time varying covariates
	real X_pat[n_obs]; // time invariant covariates

}
	
parameters
{
// intercepts
real phi_mean ; // global mean of colonization probability
real psi_mean; // initial occupancy probability

real beta_pat_mean;
real beta_time_mean;
real<lower=0> beta_pat_sd;
real<lower=0> beta_time_sd;

real beta_pat[n_strains]; // fixed across-patients covariate effect
real beta_time[ n_strains]; // fixed across


// random effects
	cholesky_factor_corr[n_strains] L_patient_phi;
	matrix[n_strains, n_patients] z_patient_phi;
	vector<lower=0> [n_strains] tau_patient_phi; // across patient std deviations
	
	cholesky_factor_corr[n_strains] L_patient_gam;
	matrix[n_strains, n_patients] z_patient_gam;
	vector<lower=0> [n_strains] tau_patient_gam; // across patient std deviations
	
	cholesky_factor_corr[n_strains] L_time_phi;
	matrix[n_strains, n_obs] z_time_phi;
	vector<lower=0> [n_strains] tau_time_phi; // within patient std deviations
	
	cholesky_factor_corr[n_strains] L_time_gam;
	matrix[n_strains, n_obs] z_time_gam;
	vector<lower=0> [n_strains] tau_time_gam; // within patient std deviations
	
	}
	
transformed parameters{
	vector[n_obs] phi; //time,patient,strain specific colonization probability
	vector[n_obs] gam; //time,patient,strain specific persistence probability
	
	//random effects
	alpha_patient_phi <- (diag_pre_multiply(tau_patient_phi, L_patient_phi) * z_patient_phi)';
	alpha_time_phi <- (diag_pre_multiply(tau_time_phi, L_time_phi) * z_time_phi)';
	
	alpha_patient_gam <- (diag_pre_multiply(tau_patient_gam, L_patient_gam) * z_patient_gam)';
	alpha_time_gam <- (diag_pre_multiply(tau_time_gam, L_time_gam) * z_time_gam)';
	
	for(i in 1:n_obs){
		phi[i] <- inv_logit(phi_mean + beta_pat[strain[i]]*X_pat[i] + alpha_patient_phi[patient[i],strain[i]] + beta_time[strain[i]]*X_time[i] + alpha_time_phi[i,strain[i]]);
		gam[i] <- inv_logit(gam_mean + beta_pat[strain[i]]*X_pat[i] + alpha_patient_gam[patient[i],strain[i]] + beta_time[strain[i]]*X_time[i] + alpha_time_gam[i,strain[i]]);
	}

model
{
  real psi[n_obs];
  
  // prior on fixed effects
  beta_pat_mean ~ normal(0,10);
  beta_time_mean ~ normal(0,10);
  beta_pat_sd ~ cauchy(0,5);
  beat_time_sd ~ cauchy(0,5);
  
  beta_patient ~ normal(beta_pat_mean, beta_pat_sd);
  beta_time ~ normal(beta_time_mean, beta_time_sd);  

  // prior on across patient random effects
  L_patient ~ lkj_corr_cholesky(2);
  tau_patient ~ cauchy(0, 3);
  to_vector(z_patient) ~ normal(0,1);
  
  // prior on within patient random effects
  L_time ~ lkj_corr_cholesky(2);
  tau_time ~ cauchy(0, 3);
  to_vector(z_time) ~ normal(0,1);
  
  // occupancy likelihood
  
  for( i in 1:n_obs){
  	// initial likelihood:
  	if(Visit_pat == 1){
  		psi[i] <- inv_logit(psi_mean);
  	}
  	
  	// subsequent likelihood:
  	if(Visit_pat >1){
		psi[i] <- gam[i-1]*psi[i-1] + phi[i-1]*(1 - psi[i-1]);
  	}
  
	Y[i] ~ bernoulli(psi[i]);
  }
  
  
  
	
} // end model block

