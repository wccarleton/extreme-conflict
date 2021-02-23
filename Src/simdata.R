Time <- 100
X <- rnorm(mean=0,n=Time,sd=1)
A <- 2
B_1 <- 0
B_2 <- 0

mu <- c()
mu[1] <- X[1] * (B_1 + (B_2 * (X[1] > A)))

l <- c()
l0 <- 1
rho = 0.9
sig <- 0.4
l[1] <- rnorm(n=1, mean = rho * l0, sd=sig)

for(j in 2:Time){
   l[j] <- rnorm(n=1, mean = rho * l[j-1],sd=sig)
   mu[j] <- X[j] * (B_1 + (B_2 * (X[j] > A)))
}

Y <- rpois(n=Time,lambda=exp(mu + l))
