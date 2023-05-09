args <- commandArgs(TRUE)
output <- args[1]
files <- args [-1]

library(yaml)
config <- yaml.load_file("config/config.yaml")

nfactors = config$nfactors
my_target = config$my_target
target = unlist(config$target)

load("EvoStatus/EvoStatus.RData")
# Calculate the number of rows and columns needed for the plot
nrows <- ceiling(sqrt(nfactors))
ncols <- ceiling(nfactors/nrows)

# Set the mfrow parameter based on the calculated values
par(mfrow=c(nrows, ncols))
pdf("Convergence.pdf")


if(my_target){
  for(index in 1:nfactors){
    plot(evolutionary_log$center[index,], type="l", ylim = c(min(target[index], min(evolutionary_log$center[index,]), na.rm = TRUE),
                                                             max(target[index], max(evolutionary_log$center[index,]), na.rm=TRUE)))
    
    # lines(evolutionary_log$est_pos[index,] , ylim=c(99000, 100500), col="green")
    if(!is.na(target[index])){abline(h=target[index], col="red", lwd=2)}
    
  }
  
}else{
  for(index in 1:3){
    plot(evolutionary_log$center[index,], type="l"
    #      , ylim = c(min(target[index], min(evolutionary_log$center[index,]), na.rm = TRUE),
    #                 max(target[index], max(evolutionary_log$center[index,]), na.rm=TRUE))
         )
    
    # lines(evolutionary_log$est_pos[index,] , ylim=c(99000, 100500), col="green")
    # if(!is.na(target[index])){abline(h=target[index], col="red", lwd=2)}
    
  }
}


dev.off()
