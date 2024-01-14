args <- commandArgs(TRUE)
output <- args[1]

library(yaml)
config <- yaml.load_file("config/config.yaml")

#how many simulations you want to run?
sim <- config$sim_init
names <- config$name_parameter
nfactors = config$nfactors

results <- matrix(0, nrow=sim, ncol=nfactors)

#random sampling of initial parameters
for(index in 1:sim){
  MIR <- rbinom(1,1,0.5)
  n_bull <- sample(300:500,1)
  n_test <- floor((10000000-3000*n_bull)/(4000 + 1000*MIR))
  n_bull_sel <- sample(15:25, 1)
## sort the variable below as it appears in config$name_parameter 
  results[index,] <- c(n_test,n_bull,n_bull_sel,MIR)
}

#sample random seed
results <- cbind(sample(1:1e5,1),results)

#add the column's name
colnames(results) <- c("randomSeed",names)


write.csv(results, file = output, row.names=TRUE)


