data {
	int<lower=0> n_strains;
	int<lower=0> n_patients;
	int<lower=0> n_obs; //total number of clinic visits across patients
	int<lower=1> visit[n_obs];
	int<lower=1> patient[n_obs];
	int<lower=0, upper=1> z[n_obs, n_strains]; //occupancy observations
}

parameters {
  matrix[n_patients, n_strains] psi0; //initial occurrence
	
  // random effects
  matrix[n_patients, n_strains * 2] alpha_raw_patient; // * 2 b/c inc. phi & gamma
  matrix[n_obs, n_strains * 2] alpha_raw_time;
  cholesky_factor_corr[n_strains * 2] L_patient; 
  cholesky_factor_corr[n_strains * 2] L_time;
}

transformed parameters{
	matrix[n_obs, n_strains] phi; // colonization
	matrix[n_obs, n_strains] gam; // persistence
	
	matrix[n_obs, n_strains] z_phi; // colonization
	matrix[n_obs, n_strains] z_gam; // persistence
	
	matrix[n_patients, n_strains * 2] alpha_patient;
	matrix[n_obs, n_strains * 2] alpha_time;
	matrix[n_obs, n_strains] psi;
	
	//row_vector[n_strains] bern_phi[n_obs];
	matrix[n_obs, n_strains] bern_phi;

	
	//random effects
	//random effects
	alpha_patient <- alpha_raw_patient * L_patient;
	alpha_time <- alpha_raw_time * L_time; 
	
	
	{// temporary scope for logit of cdf of error terms
    matrix[n_obs, n_strains] e_tmp_phi;
    matrix[n_obs, n_strains] e_tmp2_phi;
    matrix[n_obs, n_strains] e_all_phi;
    
    matrix[n_obs, n_strains] e_tmp_gam;
    matrix[n_obs, n_strains] e_tmp2_gam;
    matrix[n_obs, n_strains] e_all_gam;
 
 	for(j in 1:n_strains){
 		for (i in 1:n_obs){
    		e_tmp_phi[i,j] <- alpha_patient[patient[i],j] + alpha_time[i,j];
    		e_tmp_gam[i,j] <- alpha_patient[patient[i], j + n_strains] + alpha_time[i, j + n_strains];
    	}
    }

    //Normalize
    e_tmp2_phi <-  (e_tmp_phi - mean(e_tmp_phi))/sd(e_tmp_phi);
    e_tmp2_gam <-  (e_tmp_gam - mean(e_tmp_gam))/sd(e_tmp_gam);
    
   
 	for(j in 1:n_strains){
 		for (i in 1:n_obs){
    		e_all_phi[i,j] <- logit(Phi(e_tmp2_phi[i,j] ));
    		e_all_gam[i,j] <- logit(Phi(e_tmp2_gam[i,j] ));
    	}
    }
    
  	z_phi <- 0 + e_all_phi; // ignore fixed effects for now 
	z_gam <- 0 + e_all_gam; // ignore fixed effects for now 
	}
	
	//This is redundant, but makes things a bit clearer below
	phi <- z_phi;
	gam <- z_gam;

 
	for (k in 1:n_strains){
      	for (j in 1:n_obs){
        	if (visit[j] == 1) {
          		psi[j, k] <- psi0[j, k];
        	} else {//Current visit depends on previous visit
          		psi[j, k] <- z[j - n_patients, k] * phi[j - n_patients, k] + 
          		                (1 - z[j - n_patients, k]) * gam[j - n_patients, k];
        	}
      	}
    }
  
  // Better to do this here for efficiency
  for (k in 1:n_strains){
    for (j in 1:n_obs){
      bern_phi[j, k] <- Phi(psi[j, k]);
    }
  }
  
}// end TP block

model {
  // prior on across patient random effects
  L_patient ~ lkj_corr_cholesky(2);
  L_time ~ lkj_corr_cholesky(2);
  to_vector(alpha_raw_patient) ~ normal(0, 1);
  to_vector(alpha_raw_time) ~ normal(0, 1);

  to_vector(psi0) ~ normal(0, 3);
  
  for(k in 1:n_strains){
  	for (j in 1:n_obs){
    	z[j, k] ~ bernoulli(bern_phi[j, k]);
    }
  }
  
} // end model block

generated quantities {
  matrix[n_strains * 2, n_strains * 2] cor_patient;
  matrix[n_strains * 2, n_strains * 2] cor_time;

  cor_patient <- multiply_lower_tri_self_transpose(L_patient);
  cor_time <- multiply_lower_tri_self_transpose(L_time);

}