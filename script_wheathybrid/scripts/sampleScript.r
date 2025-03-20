args <- commandArgs(TRUE)
output <- args[1]

library(yaml)

# Read the config file
config <- yaml.load_file("config/config.yaml")
cost_par <- config$cost_par
# do we have linked parameter? which one?
linked_parameter <- config$linked_parameter
linked_parameter_adjustment <- config$linked_parameter_adjustment
is_hybrid_breeding <- config$Hybrid_Breeding

  if(is_hybrid_breeding){
    
    cheapest_unit_female <- config$cheapest_unit_female
    cheapest_unit_male <- config$cheapest_unit_male
    
    cost_cheapest_female <- config$cost_cheapest_female
    cost_cheapest_male <- config$cost_cheapest_male
    
    if(linked_parameter){
    num_linked_parameter_male <- config$num_linked_parameter_male
    num_linked_parameter_female <- config$num_linked_parameter_female
    
    linked_parameter_male <- config$linked_parameter_male
    linked_parameter_female <- config$linked_parameter_female

    param_link_male <- linked_parameter_male[-1]
    linked_param_male <- linked_parameter_male[1]
    
    param_link_female <- linked_parameter_female[-1]
    linked_param_female <- linked_parameter_female[1]
    
    }
  }

my_functions <- c("cost_function", "random_numbers")

if (any(!sapply(my_functions, exists))) {
  Path <- "./Functions"
  lapply(my_functions, function(f) source(file.path(Path, paste0(f, ".R"))))
}

# How many simulations to run for the first round?
# This usually should be big enough that we have a good coverage on search space!
sim <- config$sim_init 

# How many factors we want to optimize?
nfactors = config$nfactors
results <- matrix(0, nrow=sim, ncol= nfactors)

#random sampling of initial parameters
for(index in 1:sim){

  TC2.2_YT <- rbinom(1,1,0.5)
  ncross_female <- sample(seq(50, 200, by = 2), 1) 
  ncross_male <- sample(seq(50, 200, by = 2), 1)
  DH0_capacity_female <- sample(8000:12000,1)
  DH0_capacity_male <- sample(8000:12000,1)
  OBS1_female <- sample(400:600, 1) 
  OBS1_male <- sample(400:600, 1) 
  OBS2_female <- sample(200:400, 1) 
  OBS2_male <- sample(200:400, 1) 
  OBS3_female <- sample(15:50, 1) 
  OBS3_male <- sample(15:50, 1) 

  
  x <- c(TC2.2_YT,
        ncross_female,
        ncross_male,
        DH0_capacity_female,
        DH0_capacity_male,
        OBS1_female,
        OBS1_male,
        OBS2_female,
        OBS2_male,
        OBS3_female,
        OBS3_male
        )
  
  
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
  ####################################################################################
  ### make the contribution of recycling flexible 
  # total number of recycle parents should be the number of total cross - 20% new materials

  needed_recycle_f <- x[linked_param_female] - round(x[linked_param_female] * 0.2)
  needed_recycle_m <- x[linked_param_male] - round(x[linked_param_male] * 0.2)
  
  samples_m <- randomN(3, abs(needed_recycle_m))
  samples_f <- randomN(3, abs(needed_recycle_f))
  
  
  new_par <- c(samples_f[1], samples_f[2],samples_f[3],
               samples_m[1], samples_m[2],samples_m[3])
  x <- c(x,new_par)

  results[index,] <- x
}

################################################################################ 

#sample random seed
results <- cbind(sample(1:1e5,1),results)

#add the column's name
colnames(results) <- c("randomSeed",
                       "TC2.2_YT",
                       "ncross_female",
                       "ncross_male",
                       "DH0_capacity_female",
                       "DH0_capacity_male",
                       "P0_female",
                       "P0_male",
                       "P1_female",
                       "P1_male",
                       "P2_female",
                       "P2_male",
                       "female_rec_p2",
                       "female_rec_p1", 
                       "female_rec_p0",
                       "male_rec_p2",
                       "male_rec_p1", 
                       "male_rec_p0")

write.csv(results, file = output, row.names=TRUE)



