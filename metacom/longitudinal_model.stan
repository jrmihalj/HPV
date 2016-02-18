data {
	int<lower=0> n_strains;
	int<lower=0> n_patients;
	int<lower=0> n_obs; //total number of clinic visits across patients
	int<lower=0, upper=1> z[n_patients, n_obs, n_strains]; //occupancy observations
}

parameters {
  // random effects
  matrix<lower=0, upper=1>[n_patients, n_strains] psi0;
  real mu_psi0;
  real<lower=0> sd_psi0;
	
  // random effects
  matrix[n_patients, n_strains * 2] alpha_raw_patient; // * 2 b/c inc. phi & gamma
  matrix[n_obs, n_strains * 2] alpha_raw_time;
  cholesky_factor_corr[n_strains * 2] L_patient; 
  cholesky_factor_corr[n_strains * 2] L_time;
}

transformed parameters{
	real<lower=0, upper=1> phi[n_patients, n_obs, n_strains]; // colonization
	real<lower=0, upper=1> gam[n_patients, n_obs, n_strains]; // persistence
	
	real<lower=0, upper=1> z_phi[n_patients, n_obs, n_strains]; 
	real<lower=0, upper=1> z_gam[n_patients, n_obs, n_strains]; 
	
	matrix[n_patients, n_strains * 2] alpha_patient;
	matrix[n_obs, n_strains * 2] alpha_time;
	real<lower=0, upper=1> psi[n_patients, n_obs, n_strains];
	
	//random effects
	//random effects
	alpha_patient <- alpha_raw_patient * L_patient;
	alpha_time <- alpha_raw_time * L_time; 
	
	
	{// temporary scope for logit of cdf of error terms
    real e_tmp_phi[n_patients, n_obs, n_strains];
    real e_tmp2_phi[n_patients, n_obs, n_strains];
    real e_all_phi[n_patients, n_obs, n_strains];
    
    real e_tmp_gam[n_patients, n_obs, n_strains];
    real e_tmp2_gam[n_patients, n_obs, n_strains];
    real e_all_gam[n_patients, n_obs, n_strains];
 
 	for(k in 1:n_patients){
 		for(j in 1:n_strains){
 			for (i in 1:n_obs){
    			e_tmp_phi[i,j,k] <- alpha_patient[k,j] + alpha_time[i,j];
    			e_tmp_gam[i,j,k] <- alpha_patient[k, j + n_strains] + alpha_time[i, j + n_strains];
    }

    //Normalize
    e_tmp2_phi <-  (e_tmp_phi - mean(e_tmp_phi))/sd(e_tmp_phi);
    e_tmp2_gam <-  (e_tmp_gam - mean(e_tmp_gam))/sd(e_tmp_gam);
    
   
   for(k in 1:n_patients){
 		for(j in 1:n_strains){
 			for (i in 1:n_obs){
    			e_all_phi[i,j,k] <- logit(Phi(e_tmp2_phi[i,j,k] ));
    			e_all_gam[i,j,k] <- logit(Phi(e_tmp2_gam[i,j,k] ));
    		}
    	}
    }
    
    
    z_phi <- e_all_phi; // ignore fixed effects for now 
	z_gam <- e_all_gam; // ignore fixed effects for now 
	}
	
	
	phi <- Phi(z_phi);
	gam <- Phi(z_gam);

	for (i in 1:n_patients) { 
		for (j in 1:n_obs){
      		for (k in 1:n_strains){
        		if (j == 1) {
          			psi[i, j, k] <- inv_logit(psi0[i, k]);
        		} else {
          			psi[i, j, k] <- z[i, j - 1, k] * phi[i, j - 1, k] + (1 - z[i, j - 1, k]) * gam[i, j - 1, k];
        		}
      		}
    	}
  	}
}

model {
  // prior on across patient random effects
  L_patient ~ lkj_corr_cholesky(5);
  L_time ~ lkj_corr_cholesky(5);
  mu_psi0 ~ normal(0, 1);
  sd_psi0 ~ normal(0, 1);
  to_vector(psi0) ~ normal(mu_psi0, sd_psi0);
  
  for (i in 1:n_patients){
    for (j in 1:n_obs){
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