data
{
	int<lower=0> n_strains;
	int<lower=0> n_obs;
	int<lower=0> n_patients;
	
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

real beta_pat_mean;
real beta_time_mean;
real<lower=0> beta_pat_sd;
real<lower=0> beta_time_sd;

real beta_pat[n_strains]; // fixed across-patients covariate effect
real beta_time[ n_strains]; // fixed across


// random effects
	cholesky_factor_corr[n_strains] L_patient;
	matrix[n_strains, n_patients] z_patient;
	vector<lower=0> [n_strains] tau_patient; // across patient std deviations
	
	cholesky_factor_corr[n_strains] L_time;
	matrix[n_strains, n_obs] z_time;
	vector<lower=0> [n_strains] tau_time; // within patient std deviations
	
	}
	
transformed parameters{
	vector[n_obs] phi; //time,patient,strain specific colonization probability
	
	//random effects
	alpha_patient <- (diag_pre_multiply(tau_patient, L_patient) * z_patient)';
	alpha_time <- (diag_pre_multiply(tau_time, L_time) * time)';
	
	for(i in 1:n_obs){
		phi[i] <- phi_mean + beta_pat[strain[i]]*X_pat[i] + alpha_patient[patient[i],strain[i]] + beta_time[strain[i]]*X_time[i] + alpha_time[visit[i],strain[i]];
	}

model
{

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
  
  
	
} // end model block

