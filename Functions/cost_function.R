## define the cost of your breeding programs based on different parameters
calc_cost <- function(x){
 cost <-
   x[1] * (4000 + 1000*x[4]) +
   x[2] * 3000
  return(unname((cost)))
}

