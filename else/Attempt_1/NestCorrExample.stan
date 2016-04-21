data {
  int<lower=0> nsite; // number of sites
  int<lower=0> H; // number of host species
  int<lower=0> S; // number of symbiont species

  // host abundance data and indices
  int<lower=0> N_h; // number of host abundance observations
  int<lower=0> site_h[N_h]; // site indices for host abundance data
  int<lower=0> spec_h[N_h]; // species indices for host abundance data
  int<lower=0> y_h[N_h]; // host count data
  
  // symbiont abundance data and indices
  int<lower=0> N_s; // number of symbiont abundance observations
  int<lower=0> site_s[N_s]; // site where host capture
  int<lower=0> Hspec_s[N_s]; // host species
  int<lower=0> Sspec_s[N_s]; // symbiont species
  int<lower=0> hosts_sampled; // total number of hosts sampled w/ data
  int<lower=0, upper=hosts_sampled> host_id[N_s]; // individual host ids
  int<lower=0> y_s[N_s]; // parasite count data
}

parameters {
  // intercepts
  real<lower=0> sd_alpha_h;
  real<lower=0> sd_alpha_s;
  real mu_alpha_h;
  real mu_alpha_s;
  vector[H] alpha_h;
  vector[S] alpha_s;
  
  // random effects
  cholesky_factor_corr[H+S] L_site;
  matrix[H+S, nsite] z_site;
  vector<lower=0>[H+S] tau_site;  // among site standard deviations

  cholesky_factor_corr[S] L_host;
  matrix[S, H] z_host;
  vector<lower=0>[S] tau_host;  // among host species standard deviations

  cholesky_factor_corr[S] L_indiv;
  matrix[S, hosts_sampled] z_indiv;
  vector<lower=0>[S] tau_indiv;  // among individual standard deviations

}

transformed parameters{
  matrix[nsite, H+S] alpha_site;
  matrix[H, S] alpha_host;
  matrix[hosts_sampled, S] alpha_indiv;
  vector[N_h] loglambda_h;
  vector[N_s] loglambda_s;
  
  // random effects
  alpha_site <- (diag_pre_multiply(tau_site, L_site) * z_site)';
  alpha_host <- (diag_pre_multiply(tau_host, L_host) * z_host)';
  alpha_indiv <- (diag_pre_multiply(tau_indiv, L_indiv) * z_indiv)';
  
  for (i in 1:N_h){
    loglambda_h[i] <- alpha_h[spec_h[i]] + alpha_site[site_h[i], spec_h[i]];
  }
  
  for (i in 1:N_s){
    loglambda_s[i] <- alpha_s[Sspec_s[i]] + 
                        alpha_site[site_s[i], H + Sspec_s[i]] + 
                        alpha_host[Hspec_s[i], Sspec_s[i]] + 
                        alpha_indiv[host_id[i], Sspec_s[i]];
  }
  
}

model { 
  // intercepts
  mu_alpha_h ~ normal(0, 10);
  mu_alpha_s ~ normal(0, 10);
  sd_alpha_h ~ cauchy(0, 5);
  sd_alpha_s ~ cauchy(0, 5);
  alpha_h ~ normal(mu_alpha_h, sd_alpha_h);
  alpha_s ~ normal(mu_alpha_s, sd_alpha_s);
  
  // site level random effects
  L_site ~ lkj_corr_cholesky(2);
  tau_site ~ cauchy(0, 3);
  to_vector(z_site) ~ normal(0,1);
  
  // host-species random effects
  L_host ~ lkj_corr_cholesky(2);
  tau_host ~ cauchy(0, 3);
  to_vector(z_host) ~ normal(0, 1);
  
  // individual random effects
  L_indiv ~ lkj_corr_cholesky(2);
  tau_indiv ~ cauchy(0, 3);
  to_vector(z_indiv) ~ normal(0, 1);
  
  // host likelihood
  y_h ~ poisson_log(loglambda_h);
  
  // symbiont likelihood
  y_s ~ poisson_log(loglambda_s);
}

generated quantities {
  matrix[H+S,H+S] cor_site;
  matrix[S, S] cor_host;
  matrix[S, S] cor_indiv;
  
  cor_site <- multiply_lower_tri_self_transpose(L_site);
  cor_host <- multiply_lower_tri_self_transpose(L_host);
  cor_indiv <- multiply_lower_tri_self_transpose(L_indiv);
}