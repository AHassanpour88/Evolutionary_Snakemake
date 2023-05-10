# 
# 
# args <- commandArgs(TRUE)
# output <- args[1]
# 
# library(yaml)
# 
# # Read the config file
# config <- yaml.load_file("config/config.yaml")
# cost_par <- config$cost_par
# 
# my_functions <- c("cost_function", "dnorm_approx", "generate_new","kernel_smoothing", "random_numbers")
# 
# if (any(!sapply(my_functions, exists))) {
#   Path <- "./Functions"
#   lapply(my_functions, function(f) source(file.path(Path, paste0(f, ".R"))))
# }
# 
# # How many simulations to run for the first round?
# # This usually should be big enough that we have a good coverage on search space!
# sim <- config$sim_init 
# 
# # How many factors we want to optimize?
# nfactors = config$nfactors
# results <- matrix(0, nrow=sim, ncol= nfactors)
# 
# #random sampling of initial parameters
# for(index in 1:sim){
#   
#   HT1_2 <- rbinom(1,1,0.5)
#   HT2_2 <- rbinom(1,1,0.5)
#   DH0_capacity_female <- sample(25000:35000,1)
#   DH0_capacity_male <- sample(25000:35000,1)
#   P0_female <- sample(1000:2000, 1) 
#   P0_male <- sample(1000:2000, 1) 
#   P1_female <- sample(150:400, 1) 
#   P1_male <- sample(150:400, 1) 
#   P2_female <- sample(15:50, 1) 
#   P2_male <- sample(15:50, 1) 
#   ncross_female <- sample(seq(250, 500, by = 2), 1) #for simulation in rounding process easier to be even number
#   ncross_male <- sample(seq(250, 500, by = 2), 1)
#   
#   x <- c(HT1_2,
#          HT2_2,
#          DH0_capacity_female,
#          DH0_capacity_male,
#          P0_female,
#          P0_male,
#          P1_female,
#          P1_male,
#          P2_female,
#          P2_male,
#          ncross_female,
#          ncross_male)
#   
#   
#   base_cost <- calc_cost(c(1,1,30000, 28000, 1400, 2000, 350, 350, 20,35, 300, 350))
#   new_cost <- calc_cost(x)
#   
#   # Check basic constraints # don't need them because the maximum of them not even reach the minimum of the other
#   # if(x[5]>x[3]){ #if P1 bigger than dh0 etc.
#   #   x[5] <- x[3]
#   # }
#   # If costs are to high all parameters will be scaled down
#   #5000 as long as the cost is very different between this program and are base-line cost we need to make changes
#   #if we are really close we can stop
#   #this leads to not having the same exact cost, 
#   index2 <- 1
#   while(abs(new_cost-base_cost)>5000 & index2 < 20){
#     x[cost_par] <- round(x[cost_par] * base_cost / new_cost)
#     new_cost <- calc_cost(x)
#     index2 <- index2 + 1
#   }
#   
#   
#   if(new_cost!=base_cost){
#     cheapest_unit_female <- config$cheapest_unit_female
#     cheapest_unit_male <- config$cheapest_unit_male
#     cost_cheapest_female <- config$cost_cheapest_female
#     cost_cheapest_male <- config$cost_cheapest_male
#     if(rbinom(1,1,0.5)==1){
#       x[cheapest_unit_male] <- x[cheapest_unit_male] - round((new_cost - base_cost) / cost_cheapest_male) -1
#     } else{
#       x[cheapest_unit_female] <- x[cheapest_unit_female] - round((new_cost - base_cost) / cost_cheapest_female) -1
#     }
#   }
#   ####################################################################################
#   ### make the contribution of recycling flexible 
#   # total number of recycle parents should be the number of total cross - 20% new materials
#   
#   num_linked_parameter <- config$num_linked_parameter
#   linked_param_male <- config$linked_param_male
#   linked_param_female <- config$linked_param_female
#   
#   
#   needed_recycle_f <- x[linked_param_female] - round(x[linked_param_female] * 0.2)
#   needed_recycle_m <- x[linked_param_male] - round(x[linked_param_male] * 0.2)
#   
#   samples_m <- randomN(num_linked_parameter, needed_recycle_m)
#   samples_f <- randomN(num_linked_parameter, needed_recycle_f)
#   
#   
#   new_par <- c(samples_f[1], samples_f[2],samples_f[3],
#                samples_m[1], samples_m[2],samples_m[3])
#   x <- c(x,new_par)
#   
#   results[index,] <- x
# }
# 
# ################################################################################ 
# 
# #sample random seed
# results <- cbind(sample(1:1e5,1),results)
# 
# #add the column's name
# colnames(results) <- c("randomSeed",
#                        "HT1_2",
#                        "HT2_2",
#                        "DH0_capacity_female",
#                        "DH0_capacity_male",
#                        "P0_female",
#                        "P0_male",
#                        "P1_female",
#                        "P1_male",
#                        "P2_female",
#                        "P2_male",
#                        "ncross_female",
#                        "ncross_male",
#                        "female_rec_p2",
#                        "female_rec_p1", 
#                        "female_rec_p0",
#                        "male_rec_p2",
#                        "male_rec_p1", 
#                        "male_rec_p0")
# 
# write.csv(results, file = output, row.names=TRUE)
# 
# 
# 

args <- commandArgs(TRUE)
output <- args[1]

#how many simulations you want to run?
sim = 250

results <- matrix(0, nrow=sim, ncol=3)

#random sampling of initial parameters
for(index in 1:sim){
  n_bull <- sample(300:500,1)
  n_test <- floor((10000000-3000*n_bull)/4000)
  n_bull_sel <- sample(10:20, 1)
  results[index,] <- c(n_test,n_bull,n_bull_sel)
}

#sample random seed
results <- cbind(sample(1:1e5,1),results)

#add the column's name
colnames(results) <- c("randomSeed","n_test","n_bull","n_bull_sel")


write.csv(results, file = output, row.names=TRUE)


