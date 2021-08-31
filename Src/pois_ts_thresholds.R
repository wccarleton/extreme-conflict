poisTSCode <- nimbleCode({
    B1 ~ dnorm(mean=0,sd=100)
    B2 ~ dnorm(mean=0,sd=100)
    B3 ~ dnorm(mean=0,sd=100)
    mu[1] <- X[1] * (B1 + (B2 * (X[1] >= a1) + B3 * (X[1] <= a2)))
    a1 ~ dunif(0,5)
    d ~ dunif(1,15)
    a2 <- a1 - d
    sigma ~ dexp(0.75)
    lambda0 ~ dnorm(mean=0, sd=2)
    rho ~ dnorm(0, sd=2)
    lambda[1] ~ dnorm(rho * lambda0, sd = sigma)
    for (j in 2:Time){
      mu[j] <- X[j] * (B1 + (B2 * (X[j] >= a1) + B3 * (X[j] <= a2)))
      lambda[j] ~ dnorm(rho * lambda[j-1], sd = sigma)
    }
    for (j in 1:Time){
      Y[j] ~ dpois(exp(mu[j] + lambda[j]))
      y[j] ~ dpois(exp(mu[j] + lambda[j]))
    }
})

##
model <- "Glaser2009"

start = 1005
end = 1980

EuroClimCon <- read.csv("./Data/EuroClimCon.csv",head=T)
EuroClimCon <- subset(EuroClimCon, Year >= start & Year <= end)

Y <-EuroClimCon$Conflicts


if(model == "Glaser2009"){
    X <- scale(EuroClimCon$T_Glaser2009, scale = T)[,1]
}else if(model == "Luterbacher2016"){
    X <- scale(EuroClimCon$T_Luterbacher2016, scale = T)[,1]
}else if(model == "Buentgen2011"){
    X <- scale(EuroClimCon$T_Buntgen2011_JJA, scale = T)[,1]
}else if(model == "Buentgen2021"){
    X <- scale(EuroClimCon$T_Buntgen2021, scale = T)[,1]
}else if(model == "simulated"){
    X <- scale(EuroClimCon$T_Glaser2009, scale = T)[,1]

    Y <- simulate_ac_pois(n = length(X),
                    x = as.matrix(X),
                    b = 0,
                    l0 = 1,
                    r = 0.9,
                    s = 0.1)
}

Time <- length(Y)

niter <- 500000

poisTSData <- list(Y = Y,
                  X = X)

poisTSConsts <- list(Time = Time)

poisTSInits <- list(lambda0 = 0,
                     sigma = 1,
                     rho = 0,
                     a1 = 1,
                     d = 1,
                     B1 = 0,
                     B2 = 0,
                     B3 = 0)

poisTSModel <- nimbleModel(code = poisTSCode,
                        data = poisTSData,
                        inits = poisTSInits,
                        constants = poisTSConsts)

#compile nimble model to C++ code—much faster runtime
C_poisTSModel <- compileNimble(poisTSModel, showCompilerOutput = FALSE)

#configure the MCMC
poisTSModel_conf <- configureMCMC(poisTSModel,thin=1)
poisTSModel_conf$removeSampler(c("B1", "B2", "B3", "a1", "d"))
poisTSModel_conf$addSampler(target = c("B1", "B2", "B3"),
                            type = "AF_slice")
poisTSModel_conf$addSampler(target = c("a1", "d"),
                            type = "AF_slice")

#select the variables that we want to monitor in the MCMC chain
poisTSModel_conf$addMonitors2(c("y"))

#build MCMC
poisTSModelMCMC <- buildMCMC(poisTSModel_conf,enableWAIC=F)

#compile MCMC to C++—much faster
C_poisTSModelMCMC <- compileNimble(poisTSModelMCMC,project=poisTSModel)

#number of MCMC iterations
#niter <- 2000000

#set seed for replicability
#set.seed(1)

#call the C++ compiled MCMC model
samples <- runMCMC(C_poisTSModelMCMC, niter=niter)
save(samples,file=paste("./Results/mcmc_samples_",model,"_y.RData",sep=""))

alarm()
