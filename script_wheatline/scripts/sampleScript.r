args <- commandArgs(TRUE)
output <- args[1]


library(yaml)
config <- yaml.load_file("config/config.yaml")

#how many simulations you want to run?
sim <- config$sim_init
names <- config$name_parameter
nfactors = config$nfactors

base_cost_ini <- unlist(config$base_cost_ini) 
cheapest_unit <- config$cheapest_unit
cost_cheapest <- config$cost_cheapest
cost_par <- config$cost_par 

my_functions <- c("cost_function")

if (any(!sapply(my_functions, exists))) {
  Path <- "./Functions"
  lapply(my_functions, function(f) source(file.path(Path, paste0(f, ".R"))))
}


results <- matrix(0, nrow=sim, ncol=nfactors)


#random sampling of initial parameters
for(index in 1:sim){
  
  valid <- FALSE
  attempt <- 1
  
  # check all parameters are reasonable if not valid at some points it set to false and rerun everything again
  while(valid==FALSE){
    
    nCrosses <- sample(50:250,1)
    nDH <- sample(50:nCrosses,1)
    nPYT <- sample(300:700,1)
    nAYT <- sample(20:80,1)
    nEYT <- sample(5:25,1)
    newParents_replace <-  sample(5:50,1)
    
    x <- c(nCrosses,nDH,nPYT,nAYT,nEYT,newParents_replace)
    
    
    base_cost <- calc_cost(base_cost_ini)     
    new_cost <- calc_cost(x)
    
    index2 <- 1
    while(abs(new_cost-base_cost)>5000 & index2 < 20){
      x[cost_par] <- round(x[cost_par] * base_cost / new_cost)
      new_cost <- calc_cost(x)
      index2 <- index2 + 1
      
    }
    
    
    if(new_cost!=base_cost){
      if(config$cheapest_unit > length(x)){
        stop("cheapest_unit of the config file has to be <= ",length(x))
      }
      
      x[cheapest_unit] <- x[cheapest_unit] - round((new_cost - base_cost) / cost_cheapest) -1
      
    }
    
    
    is_copy <- sum(colSums(t(results)==x)==nfactors)>0
    
    if(!is_copy){
      valid = TRUE
    }
    
    
    if(valid){
      if (x[1] > 1 && x[2] > 1 && x[3] > 1 && x[4] > 5 && x[5] > 5 && x[6] > 5 &&
          x[1] > x[2] && x[3] > x[4] && x[4] > x[5] && x[3] <= 900 && 
          x[6] <= 50 && x[1] <= 500){
        print("Conditions for contrains are satisfied.")
        valid = TRUE
        attempt <- 1
      } else {
        valid = FALSE
        attempt <- attempt + 1
        print("Conditions for contrains are not satisfied.")
      }
    }

    

    
  } 
  
  
  
  results[index,] <- x
  
}

sum(duplicated(results))
#sample random seed
results <- cbind(sample(1:1e5,1),results)

#add the column's name
colnames(results) <- c("randomSeed",names)


write.csv(results, file = output, row.names=TRUE)


