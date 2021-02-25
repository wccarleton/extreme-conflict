Time <- 200
X <- rnorm(mean=0,n=Time,sd=1)
ChronoUnc <- 0.5
K <- 5
Ts <- do.call(rbind,lapply(1:Time,function(x)rnorm(n=K,mean=x,sd=ChronoUnc)))
apply(Ts,2,function(j)approx(x=j,y=X,xout=1:Time))

Xs <- apply(Ts,2,function(j){
        sample_X=approx(x=j,y=X,xout=1:Time)
        return(sample_X$y)})

A1 <- 1.5
A2 <- -1.5
B_0 <- 1
B_1 <- 0
B_2 <- 1.5
B_3 <- 0

#background autocorrelation
l <- c()
l0 <- 0.5
rho = 0.9
sig <- 0.2
l[1] <- rnorm(n=1, mean = rho * l0, sd=sig)
for(j in 2:Time){
    l[j] <- rnorm(n=1, mean = rho * l[j-1],sd=sig)
}

mu <- c()
for(j in 1:Time){
   mu[j] <- B_0 + X[j] * (B_1 + (B_2 * (X[j] > A1) + B_3 * (X[j] < A2)))
}

Y <- rpois(n=Time,lambda=exp(mu + l))
