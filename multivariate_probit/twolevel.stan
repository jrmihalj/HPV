data {
  int n;
  int d;
  //int<lower = 0, upper = 1> y[n, d];
  matrix[n, d] y;
  real<lower = 1> eta;
  int n_patient;
  int n_time; // number of visits 
  int<lower=1, upper=n_time> time[n]; // visit number 
  int<lower = 0, upper = n_patient> patient[n]; //patient number 
  vector[2] dir_prior;
}

transformed data {
  matrix[n, d] sign;

  for (i in 1:n) {
    for (j in 1:d) {
      if (y[i, j] == 1) {
        sign[i, j] <- 1;
      } else {
        sign[i, j] <- -1;
      }
    }
  }
}

parameters {
  cholesky_factor_corr[d] L_Rho_patient;
  cholesky_factor_corr[d] L_Rho_visit;
  matrix<lower = 0>[n, d] abs_ystar;
  matrix[n_patient, d] e_patientR;
  simplex[2] var_mat[d]; // patient level and obs-level random effect variances : simplex sum to one 
  matrix[d, d] beta_phi;
  matrix[d, d] beta_gam_raw;
   
}

transformed parameters {
  cholesky_factor_cov[d] L_Sigma_visit;
  matrix[n_patient, d] e_patient;
  vector<lower = 0>[d] sd_visit;
  vector<lower = 0>[d] sd_patient;
  matrix[n, d] mu_phi;
  matrix[n, d] mu_gam;
  matrix[d, d] beta_gam;
  matrix[n, d] fixed_effect;
  
  for(j in 1:d){
    for(i in 1:d){
      if(i==j){
        //This is because doesn't make sense for intra-specific effects on colonization.
        //By definition, species is not present in t-1
        beta_gam[i,j] <- 0;
      }else{
        beta_gam[i,j] <- beta_gam_raw[i,j];
      }
    }
  }
  
  mu_phi <- y * beta_phi;
  mu_gam <- y * beta_gam;
  
  for (j in 1:d){
      for (i in 1:n){
        if (time[i] == 1) {
            fixed_effect[i, j] <- 0;
        } else {//Current visit depends on previous visit
            fixed_effect[i, j] <- (y[i - n_patient, j]) * mu_phi[i-n_patient, j] + //persistence
            		       (1 - y[i - n_patient, j]) * mu_gam[i-n_patient, j];  //colonization
        }
      }
	}

  for (i in 1:d) {
    sd_visit[i] <- sqrt(var_mat[i, 1]);
    sd_patient[i] <- sqrt(var_mat[i, 2]);
  }
  e_patient <- (diag_pre_multiply(sd_patient, L_Rho_patient) * e_patientR')';
  L_Sigma_visit <- diag_pre_multiply(sd_visit, L_Rho_visit);
}

model {
  for (i in 1:d) {
    var_mat[i] ~ dirichlet(dir_prior); // dirichlet prior for the variances - why?
  }

  L_Rho_patient ~ lkj_corr_cholesky(eta); // specify tightness of corr prior 
  L_Rho_visit ~ lkj_corr_cholesky(eta);
  to_vector(e_patientR) ~ normal(0, 1);
  for (i in 1:n) {
    sign[i, ] .* abs_ystar[i, ] ~ multi_normal_cholesky((fixed_effect[i,] + e_patient[patient[i], ]),
                                                          L_Sigma_visit);
  }
}

generated quantities {
  matrix[d, d] Rho_patient;
  matrix[d, d] Rho_visit;

  Rho_patient <- multiply_lower_tri_self_transpose(L_Rho_patient);
  Rho_visit <- multiply_lower_tri_self_transpose(L_Rho_visit);
}
