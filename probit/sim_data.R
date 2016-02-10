# generative model from sebastian gonzalez paper
require(clusterGeneration)

n <- 400 # sampling units
m <- 5 # species

# among site corr matrix
R <- genPositiveDefMat(dim=m, #Number of columns/rows
                       covMethod = 'onion', 
                       rangeVar = c(1, 1), # Range of variances
                       eta=2)$Sigma

# among site ranefs
e <- matrix(nrow = n, ncol = m)
for (i in 1:n){
  e[i, ] <- rnorm(m) %*% t(chol(R))
}

k <- 1
X <- matrix(1, nrow = n, ncol = k)
beta <- matrix(rnorm(m * k), nrow = k)

mu <- X %*% beta

z <- mu + qlogis(pnorm(e))

y <- ifelse(z > 0, 1, 0)
