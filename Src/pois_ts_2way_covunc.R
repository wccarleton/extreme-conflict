poisTSCode <- nimbleCode({
   B_1 ~ dnorm(mean=0,sd=2)
   B_2 ~ dnorm(mean=0,sd=1)
   B_3 ~ dnorm(mean=0,sd=1)
   mu[1] <- X[1] * (B_1 + (B_2 * (X[1] > a1) + B_3 * (X[1] < a2)))
   a1 ~ dunif(-100,100)
   a2 ~ dunif(-100,100)
   sigma ~ dexp(0.1)
   lambda0 ~ dnorm(mean=0,sd=2)
   rho ~ dnorm(0,sd=2)
   lambda[1] ~ dnorm(rho * lambda0,sd=sigma)
   for (j in 2:Time){
      mu[j] <-  X[j] * (B_1 + (B_2 * (X[j] > a1) + B_3 * (X[j] < a2)))
      lambda[j] ~ dnorm(rho * lambda[j-1],sd=sigma)
   }
   for (j in 1:Time){
      Y[j] ~ dpois(exp(mu[j] + lambda[j]))
   }
})

##
K <- 1#ncol(X)

poisTSData <- list(Y=Y,
                  X=X)

poisTSConsts <- list(Time=Time,
                     K=K)

poisTSInits <- list(lambda0=0,
                     sigma=1,
                     rho=0,
                     a1=1,
                     a2=1,
                     B_1=rep(0,K),
                     B_2=rep(0,K),
                     B_3=rep(0,K))

poisTSModel <- nimbleModel(code=poisTSCode,
                        data=poisTSData,
                        inits=poisTSInits,
                        constants=poisTSConsts)

#compile nimble model to C++ code—much faster runtime
C_poisTSModel <- compileNimble(poisTSModel, showCompilerOutput = FALSE)

#configure the MCMC
poisTSModel_conf <- configureMCMC(poisTSModel,thin=1)
poisTSModel_conf$removeSampler(c("B_1","B_2","B_3","a1","a2"))
poisTSModel_conf$addSampler(target = c("B_1", "B_2", "B_3"), type = "AF_slice")
poisTSModel_conf$addSampler(target = c("a1","a2"), type = "RW_block")

#select the variables that we want to monitor in the MCMC chain
#poisTSModel_conf$addMonitors(c("mu"))

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
