generate_new <- function(){

  x <- numeric(nfactors)
  
  x[1] <- sample(300:500,1)
  x[2] <- sample(15:25,1)
  x[3] <- sample(280:350,1)
  x[4] <- sample(70:100,1)
  x[5] <- sample(18:22,1)
  x[6] <- sample(35:40,1)

  x <- c(x[1],x[2],x[3],x[4],x[5],x[6])
  
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
  

  return(x)
}
