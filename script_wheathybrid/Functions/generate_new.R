generate_new <- function(){

  x <- numeric(nfactors)
  
  x[1] <- rbinom(1,1,0.5)
  x[2] <- sample(seq(50, 200, by = 2), 1)
  x[3] <- sample(seq(50, 200, by = 2), 1)
  x[4] <- sample(8000:12000, 1)
  x[5] <- sample(8000:12000, 1)
  x[6] <- sample(400:600, 1)
  x[7] <- sample(400:600, 1)
  x[8] <- sample(200:400, 1)
  x[9] <- sample(200:400, 1)
  x[10] <- sample(15:50, 1)
  x[11] <- sample(15:50, 1)

  
  base_cost <- calc_cost(c(1, 100, 100, 10000, 10000, 500, 500, 300, 300, 30, 30))
  new_cost <- calc_cost(x)
  
  index2 <- 1
  while(abs(new_cost-base_cost)>5000 & index2 < 20){
    x[cost_par] <- round(x[cost_par] * base_cost / new_cost)
    new_cost <- calc_cost(x)
    index2 <- index2 + 1
  }
  

  if(new_cost!=base_cost){
    cheapest_unit_female <- config$cheapest_unit_female
    cheapest_unit_male <- config$cheapest_unit_male
    cost_cheapest_female <- config$cost_cheapest_female
    cost_cheapest_male <- config$cost_cheapest_male
    if(rbinom(1,1,0.5)==1){
      x[cheapest_unit_male] <- x[cheapest_unit_male] - round((new_cost - base_cost) / cost_cheapest_male) -1
    } else{
      x[cheapest_unit_female] <- x[cheapest_unit_female] - round((new_cost - base_cost) / cost_cheapest_female) -1
    }
  }  

  needed_recycle_f <- x[2] - round(x[10] * 0.2)
  needed_recycle_m <- x[3] - round(x[11] * 0.2)
  
  samples_m <- randomN(3, needed_recycle_m)
  samples_f <- randomN(3, needed_recycle_f)
  
  x[12] <- samples_f[1]
  x[13] <- samples_f[2]
  x[14] <- samples_f[3]
  x[15] <- samples_m[1]
  x[16] <- samples_m[2]
  x[17] <- samples_m[3]
  
  return(x)

}
