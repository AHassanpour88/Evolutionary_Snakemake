## define the cost of your breeding programs based on different parameters
calc_cost <- function(x){
  
  cost <-
    # nCrosses 
    x[1] * 30 + # 100 * 30
    # Grow F1s 
    x[1] * 30 + # 100 * 30
    # nDH -> size nCrosses*nDH
    x[2] * x[1] * 30 + 
    # genotype
    x[2] * x[1] * 15 + 
    # nPYT
    x[3] * 5 * 20 + 
    # nAYT
    x[4] * 15 * 50 +
    # nAYT
    x[5] * 20 * 50 
  return(unname((cost)))
}



