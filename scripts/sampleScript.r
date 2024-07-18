#################################################################################
################################# Agreement #####################################
#################################################################################
# Authors
# Azadeh Hassanpour, azadeh.hassanpour@uni-goettingen.de
# Johannes Geibel, johannes.geibel@fli.de
# Torsten Pook, Torsten.pook@wur.nl
# Copyright © 2020 – 2024
# This program falls under a NonCommercial-NoDerivates-NoDistriubtion Public License.
# With use of these scripts, I confirm that I represent an academic institute and acknowledge that I shall use this script only for research. 
# I explicitly acknowledge the terms in the license agreement https://github.com/AHassanpour88/Evolutionary_Snakemake/blob/main/License.md.
# I understood that any commercial use needs a commercial license from the owner of the script. For more information about a commercial license please contact Torsten Pook.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#################################################################################
#################################################################################
#################################################################################

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
  n_bull <- sample(100:700,1)
  n_test <- floor((10000000-3000*n_bull)/(4000 + 1000*MIR))
  n_bull_sel <- sample(3:30, 1)
## sort the variable below as it appears in config$name_parameter 
  results[index,] <- c(n_test,n_bull,n_bull_sel,MIR)
}

#sample random seed
results <- cbind(sample(1:1e5,1),results)

#add the column's name
colnames(results) <- c("randomSeed",names)


write.csv(results, file = output, row.names=TRUE)


