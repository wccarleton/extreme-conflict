library(coda)
models <- c("Glaser2009",
            "Buentgen2011",
            "Buentgen2021",
            "Luterbacher2016")

geweke_z <- c()

for(j in models){
    p <- paste("./mcmc_samples_", j, ".RData", sep="")
    load(p)
    g <- geweke.diag(samples[-c(1:10000),])$z
    geweke_z <- rbind(geweke_z, g)
}

geweke_results <- data.frame(model = models)
geweke_results <- cbind(geweke_results, geweke_tests)

write.csv(geweke_results, file="./geweke_results.csv",row.names=F)
