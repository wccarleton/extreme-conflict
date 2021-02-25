poisTSCode <- nimbleCode({
    ###top-level regression parameters
    B0 ~ dnorm(mean=0,sd=10)
    B1 ~ dnorm(mean=0,sd=10)
    B1ut ~ dnorm(mean=0,sd=5)
    B1lt ~ dnorm(mean=0,sd=5)
    sigB0 ~ dunif(1e-10,100)
    sigB1 ~ dunif(1e-10,100)
    sigB1ut ~ dunif(1e-10,100)
    sigB1lt ~ dunif(1e-10,100)
    ###upper and lower thresholds
    ut ~ dunif(-100,100)
    lt ~ dunif(-100,100)
    ###background autocorrelation
    sigma ~ dexp(0.1)
    lambda0 ~ dnorm(mean=0,sd=2)
    rho ~ dnorm(0,sd=2)
    lambda[1] ~ dnorm(rho * lambda0,sd=sigma)
    for (j in 2:Time){
        lambda[j] ~ dnorm(rho * lambda[j-1],sd=sigma)
    }
    ###low-level regression parameters
    for (k in 1:K){
        b0[k] ~ dnorm(mean=B0,sd=sigB0)
        b1[k] ~ dnorm(mean=B1,sd=sigB1)
        b1ut[k] ~ dnorm(mean=B1ut,sd=sigB1ut)
        b1lt[k] ~ dnorm(mean=B1lt,sd=sigB1lt)
        for (j in 1:Time){
            mu[j,k] <-  b0[k] + X[j,k] * (b1[k] + (b1ut[k] * (X[j,k] > ut) + b1lt[k] * (X[j,k] < lt)))
            Y[j,k] ~ dpois(exp(mu[j,k] + lambda[j]))
        }
    }
})

## Data
K <- ncol(Xs)
Y <- matrix(rep(Y,K),ncol=K)

poisTSData <- list(Y=Y[3:196,],
                  X=Xs[3:196,])

poisTSConsts <- list(Time=length(Y[3:196]),
                     K=K)

poisTSInits <- list(lambda0=0,
                     sigma=1,
                     rho=0,
                     ut=1,
                     lt=1,
                     B0=0,
                     B1=0,
                     B1ut=0,
                     B1lt=0,
                     sigB0=0.1,
                     sigB1=0.1,
                     sigB1ut=0.1,
                     sigB1lt=0.1,
                     b0=rep(0,K),
                     b1=rep(0,K),
                     b1ut=rep(0,K),
                     b1lt=rep(0,K))

poisTSModel <- nimbleModel(code=poisTSCode,
                        data=poisTSData,
                        inits=poisTSInits,
                        constants=poisTSConsts)

#compile nimble model to C++ code—much faster runtime
C_poisTSModel <- compileNimble(poisTSModel, showCompilerOutput = FALSE)

#configure the MCMC
poisTSModel_conf <- configureMCMC(poisTSModel,thin=1)
poisTSModel_conf$removeSampler(c(
                                "B0",
                                "B1",
                                "B1ut",
                                "B1lt",
                                "ut",
                                "lt",
                                "b0",
                                "b1",
                                "b1ut",
                                "b1lt"))
poisTSModel_conf$addSampler(target = c("B0", "B1", "B1ut","B1lt"), type = "AF_slice")

poisTSModel_conf$addSampler(target = c("ut","lt"), type = "RW_block")

for(k in 1:K){
   poisTSModel_conf$addSampler(target=c(paste("b0[",k,"]",sep=""),paste("b1[",k,"]",sep=""),paste("b1ut[",k,"]",sep=""),paste("b1lt[",k,"]",sep="")),type="AF_slice")
}

#select the variables that we want to monitor in the MCMC chain
poisTSModel_conf$monitors <- c("B0","B1","B1ut","B1lt","ut","lt","sigB0","sigB1","sigB1ut","sigB1lt")
poisTSModel_conf$addMonitors2(c("b0","b1","b1ut","b1lt","rho","sigma"))

#build MCMC
poisTSModelMCMC <- buildMCMC(poisTSModel_conf,enableWAIC=F)

#compile MCMC to C++—much faster
C_poisTSModelMCMC <- compileNimble(poisTSModelMCMC,project=poisTSModel)

#number of MCMC iterations
niter <- 400000

#set seed for replicability
set.seed(1)

#call the C++ compiled MCMC model
samples <- runMCMC(C_poisTSModelMCMC, niter=niter)

alarm()
