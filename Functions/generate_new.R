# generate_new <- function(){
#   x <- numeric(12)
#   x[1] <- rbinom(1,1,0.5)
#   x[2] <- rbinom(1,1,0.5)
#   x[3] <- sample(25000:40000,1)
#   x[4] <- sample(25000:40000,1)
#   x[5] <- sample(1000:3000, 1)
#   x[6] <- sample(1000:3000, 1)
#   x[7] <- sample(150:400, 1)
#   x[8] <- sample(150:400, 1)
#   x[9] <- sample(15:50, 1)
#   x[10] <- sample(15:50, 1)
#   x[11] <- sample(250:500, 1)
#   x[12] <- sample(250:500, 1)
#   
#   new_cost <- calc_cost(x)
#   
#   
#   # If costs are to high all parameters will be scaled down
#   #5000 as long as the cost is very different between this program and are base-line cost we need to make changes
#   #if we are really close we can stop
#   #this leads to not having the same exact cost,
#   index <- 1
#   while(abs(new_cost-base_cost)>5000 & index < 20){
#     x[cost_par] <- round(x[cost_par] * base_cost / new_cost)
#     new_cost <- calc_cost(x)
#     index <- index + 1
#   }
#   
#   
#   # For the fine tuning number will be adopted as they are least costly
#   # Cheapest unit in breeding is DH ones to compensate the rest for scaling
#   if(new_cost!=base_cost){
#     if(rbinom(1,1,0.5)==1){
#       x[cheapest_unit_male] <- x[cheapest_unit_male] - round((new_cost - base_cost) / cost_cheapest_male) -1
#     } else{
#       x[cheapest_unit_female] <- x[cheapest_unit_female] - round((new_cost - base_cost) / cost_cheapest_female) -1
#     }
#   }
#   
#   return(x)
# }


generate_new <- function(){

x <- numeric(3)
x[1] <- sample(300:500,1)
x[2] <- floor((10000000-3000)/4000)
x[3] <- sample(10:20, 1)

x <- c(x[2],x[1],x[3])

  return(x)
}
