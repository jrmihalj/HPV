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
	
	vector[n_visits] X[n_patients,n_strains]; // 3D array holding profile at each visit with each strain : could do this in the stan model 
	
}
	
parameters 
{
  // intercepts
  real<lower = 0> psi_initial; // initial occupancy probability: equal for all strains
  
  // params of fixed effects 
  vector[n_strains] betas[n_strains];
  real<lower = 0> beta_mean; //params for betas (Laplace distributed)
  real<lower = 0> beta_scale;

  // random  patient effects
  cholesky_factor_corr[n_strains] L_patient; // 
  matrix[n_strains , n_patients] z_patient; //
  vector<lower=0> [n_strains] tau_patient; // across patient std deviations

}

transformed parameters
{
	
	// correlations 
	matrix[n_patients, n_strains] alpha_patient;
	vector[n_strains] cov_effects[n_patients, n_strains];	

	//random effects
	alpha_patient <- (diag_pre_multiply(tau_patient, L_patient) * z_patient)';
	
	for(i in 1:n_obs){
	//temporary params
	real lpsi;
	real covs;
	  if(visit_pat[i] == 1){
  	  	psi[i] <- inv_logit(psi_initial);
  	  } // end if
  	  
  	  else {
  	 	 covs <- X[visit_pat[i] - 1]; // n_pat x n_strain profile for previous visit 
		 cov_effect <- dot_product(betas[strain[i]], covs[patient[i]]);  
  	  	 lpsi <- logit(psi_initial) + cov_effect + alpha_patient[patient[i],strain[i]];
  	 	 psi[i] <- inv_logit(lpsi);
  	  }// end else
	}// end n_obs
}

model {
	
  // occupancy probability 
  vector<lower=0, upper=1>[n_obs] psi;
  
  //Priors: 
  //prior on psi 
  psi_initial ~ normal(0, .5);
  
  // prior on fixed covariate effects
  beta_mean ~ normal(0, 1.5);
  //beta_scale -???
  beta ~ double_exponential(beta_mean, beta_sigma);
  
  // prior on patient random effects
  L_patient ~ lkj_corr_cholesky(2);
  tau_patient ~ normal(0, 2);
  to_vector(z_patient) ~ normal(0, 1);

  // likelihood of occupancy
  for( i in 1:nObs){
  	if( i == 1) {
  		Y[i] ~ bernoulli(psi)
  } // end if 

} // end model block

generated quantities {
  matrix[n_strains , n_strains ] cor_patient;
  
  // create the correlation matrix to draw random effects
  cor_patient <- multiply_lower_tri_self_transpose(L_patient);
  }  