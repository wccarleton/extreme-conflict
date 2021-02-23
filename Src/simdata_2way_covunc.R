Time <- 200
X <- rnorm(mean=0,n=Time,sd=1)
Ts <- do.call(rbind,lapply(1:200,function(x)rnorm(n=50,mean=x,sd=1)))
apply(Ts,2,function(j)approx(x=j,y=X,xout=1:200))

Xs <- apply(Ts,2,function(j){
        sample_X=approx(x=j,y=X,xout=1:200)
        return(sample_X$y)})

A1 <- 1.5
A2 <- -1.5
B_1 <- 0
B_2 <- 0
B_3 <- 0

mu <- c()
mu[1] <- X[1] * (B_1 + (B_2 * (X[1] > A1) + B_3 * (X[1] < A2)))

l <- c()
l0 <- 1
rho = 0.9
sig <- 0.4
l[1] <- rnorm(n=1, mean = rho * l0, sd=sig)

for(j in 2:Time){
   l[j] <- rnorm(n=1, mean = rho * l[j-1],sd=sig)
   mu[j] <- X[j] * (B_1 + (B_2 * (X[j] > A1) + B_3 * (X[j] < A2)))
}

Y <- rpois(n=Time,lambda=exp(mu + l))
