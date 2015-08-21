data
{
	int n_patients;
    int n_visits; // max num visits 
    // int n_missing_values; we can add imputation for missing covariate responses when we start working with the real data 
    int n_strains;
    int n_covariates;
    
    real valid_observation [n_strains, n_individuals, n_visits]; // binary indicator of if this person was observed for this visit/strain
    real observations[n_strains, n_individuals, n_visits]; // infection status strian i, individual k, time t
	vector[n_covariates] constant_covariate_values[n_subjects]; // covs that don't vary across time
	vector[n_covariates] time_varying_covariate_values[n_subjects, n_visits]; // covs that vary across time 
	
parameters
{
	vector[n_covariates] beta_covariates;
	real<lower=0, upper=1> p_occupancy[n_strains, n_individuals, n_visits];
	real<lower=0, upper=1> p_persistence[n_strains, n_individuals, n_visits];
	real<lower=0, upper=1> p_colonization[n_strains, n_individuals, n_visits];
	
}

model
{
	for(i in 1:n_strains){
		for( j in 1:n_individuals{
			for( k in 1:n_visits){
				if(k == 1){
					observations[i,j,k] ~ bernoulli(p_occupancy[i,j,k]);
				} // end if 
				
				if(k > 1){
					p_occupancy[i,j,k] <- p_occoupancy[i,j,k-1]*p_persistence[i,j,k-1] + (1 - p_occupancy[i,j,k-1])*p_colonization[i,j,k-1];
					observations[i,j,k] ~ bernoulli(p_occupancy[i,j,k]);
				} // end if 
			} // end visits
		} // end individuals	
	} // end strains 
	
} // end model block

generated quantities{

}