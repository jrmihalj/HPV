library(rstan)
library(MASS)
mcode <- '
  data{
    vector[2] x[100];
  }
  transformed data{
    vector [2]  mu_ln;
    vector [2]  sigma_ln;
    real <lower=0> eta;
    mu_ln[1] <- 0.0;
    mu_ln[2] <- 0.0;
    sigma_ln[1] <- 1.0;
    sigma_ln[2] <- 1.0;
    eta <- 1.0;
  }
  parameters{
    vector[2] mu;
    cov_matrix[2] Sigma;
  }
  model {
  Sigma ~ lkj_cov(eta);
    for (j in 1:2){
    mu[j]~normal(0,1000);
  }
  for(j in 1:100){
    x[j]~multi_normal(mu, Sigma);
  }
  }
'
Sigma <-matrix(c(10,3,3,2),2,2)
x<-mvrnorm(n=100, c(-7,3), Sigma)

datalist = list(
x=x
)

modfit <- stan(model_code = mcode, data = datalist,
            iter = 1000, chains = 4)
