library(rstan)
library(MASS)
mcode <- '
data
{
  int n_strains;
  int n_obs;
  //int occupancy[n_obs, n_strains];
  vector[n_strains] patient_vals[n_obs];
	vector[n_strains] time_vals[n_obs];
}

transformed  data
{
  //params for prior on cov matrix: 0 mean, log normal prior on standard deviations
  vector [n_strains]  mu_ln_across;
  vector [n_strains]  sigma_ln_across;
  vector [n_strains]  mu_ln_within;
  vector [n_strains]  sigma_ln_within;
  real <lower=0> eta;
  for( i in 1:n_strains){
    mu_ln_across[i] <- 0.0;
    sigma_ln_across[i] <- 1.0;
    mu_ln_within[i] <- 0.0;
    sigma_ln_within[i] <- 1.0;
  }

  eta <- 1.0;
}

parameters
{

  vector[n_strains] b_0_phi;
  vector[n_strains] b_1t_phi;
  
  cov_matrix[n_strains] sig_across_phi;
  cov_matrix[n_strains] sig_within_phi;	

  vector[n_strains] beta_pat[n_obs];
  vector[n_strains] beta_time[n_obs];
}

model
{
  sig_across_phi ~ lkj_cov(mu_ln_across, sigma_ln_across, eta);
  sig_within_phi ~ lkj_cov(mu_ln_within, sigma_ln_within, eta);

  for( j in 1:n_strains){
    b_0_phi[j] ~ normal(0,1000); // relatively uninformative prior on mean 
    b_1t_phi[j] ~ normal(0,1000);
  }

  for( j in 1:n_obs){
    beta_pat[j] ~ multi_normal(b_0_phi, sig_across_phi);
    beta_time[j] ~ multi_normal(b_1t_phi, sig_within_phi);
  }
} // end model block

'


datalist = list(
  n_strains=n_strains,
  n_obs = n_obs,
  patient_vals = patient_vals,
  time_vals = time_vals
)

modfit <- stan(model_code = mcode, data = datalist,
            iter = 1000, chains = 1)

