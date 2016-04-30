# generating data from a multivariate probit model
n <- 300 # number of observations
d <- 2 # number of dimensions

corr <- runif(1, -.95, .95)
Rho <- matrix(c(1, corr, corr, 1), nrow = 2)
L_Rho <- chol(Rho)

y_star <- array(rnorm(n * d), dim = c(n, d)) %*% L_Rho
cor(y_star)

plot(y_star)

y <- ifelse(y_star > 0, 1, 0)

