data {
  int n;
  int d;
  int<lower = 0, upper = 1> y[n, d];
  real<lower = 1> eta;
}

transformed data {
  matrix[n, d] sign;
  vector[d] mu;

  for (i in 1:n) {
    for (j in 1:d) {
      if (y[i, j] == 1) {
        sign[i, j] <- 1;
      } else {
        sign[i, j] <- -1;
      }
    }
  }

  for (j in 1:d) {
    mu[j] <- 0;
  }
}

parameters {
  corr_matrix[d] Rho_prior;
  cholesky_factor_corr[d] L_Rho;
  matrix<lower = 0>[n, d] abs_ystar;
}

model {
  L_Rho ~ lkj_corr_cholesky(eta);
  Rho_prior ~ lkj_corr(eta);

  // likelihood: latent variable probit model
  for (i in 1:n) {
    sign[i, ] .* abs_ystar[i, ] ~ multi_normal_cholesky(mu, L_Rho);
  }
}

generated quantities {
  matrix[d, d] Rho;

  Rho <- multiply_lower_tri_self_transpose(L_Rho);
}
