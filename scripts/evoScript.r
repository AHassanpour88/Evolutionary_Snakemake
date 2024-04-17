################## Do not change #################
args <- commandArgs(TRUE)
output <- args[1]
iteration <- as.integer(args[2])
files <- args [-1:-2]

sessionInfo()

library(yaml)
library(data.table)
library(stringr)

target_density = 30 ## in later iteration we are modifying some steps to achieve convergence nicer
################################################################################
################ Read the config file and Selection Rates ######################
################################################################################
config <- yaml.load_file("config/config.yaml") ### read the parameters of user
evoParameters <- fread(config$evoParams) ### read the selection intensity and number of parents and offspring that generates in each iteration

is_hybrid_breeding <- config$Hybrid_Breeding
is_line_breeding <- config$Line_Breeding

# Selection pressure can be adjusted by varying the candidate pool size in different iteration intervals, which can be defined in the 'config_iterations.csv' table. It is recommended to start with a larger candidate pool in the early generations and gradually reduce the number of candidates as the evolution progresses
sel_tmp_row <- which(evoParameters$iteration_start <= iteration & evoParameters$iteration_end >= iteration)

# we recommend to select most of the parents from top performing candidates
n_sel_high_target <- evoParameters$n_sel_high_target[sel_tmp_row]
# we recommend to select at least 10% of the parents from high performing candidates using kernel regression method
n_sel_kernel <- evoParameters$n_sel_kernel[sel_tmp_row]
# Use some high performing candidates from second last pool to don't loose diversity
n_sel_last2 <- evoParameters$n_sel_last2[sel_tmp_row]
# the number of total parents and also the first batch of offsprings will be exactly these parents
n_sel <- n_sel_high_target + n_sel_kernel + n_sel_last2
n_off1 <- n_sel

# how many offspring should be made with recombination approach
n_off_recombination <- evoParameters$n_off_recombination[sel_tmp_row]
# how many offspring should be made with mutation approach
n_off_mutation <- evoParameters$n_off_mutation[sel_tmp_row]
# how many offspring should be made with random approach
n_off_random <- evoParameters$n_off_random[sel_tmp_row]

# how many times simulation sould be replicated?
n_rep <- evoParameters$nrep[sel_tmp_row]

#how many parameters we have?
nfactors <- config$nfactors

#treshold for termination criteria
thresh_target_ave <- config$thresh_target_ave
thresh_target_sd <- config$thresh_target_sd

#how many binary parameters we have? Select the binary variables
number_binary_parameter <- config$number_binary_parameter
binary_parameter <- config$binary_parameter
is_binary <- binary <- config$binary
is_integer <- config$integer
cost_par <- config$cost_par

# do we have linked parameter? which one?
linked_parameter <- config$linked_parameter
linked_parameter_adjustment <- config$linked_parameter_adjustment

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

    if(is_line_breeding){
    
    cheapest_unit <- config$cheapest_unit
    cost_cheapest <- config$cost_cheapest
    
    if(linked_parameter){
    num_linked_parameter <- config$num_linked_parameter
    linked_parameter_line <- config$which_linked_parameter
    
    param_link <- linked_parameter_line[-1]
    linked_param <- linked_parameter_line[1]

    }
  }
# minimum amount of mutation that you allow for each parameter
min_range_mut <- unlist(config$min_range_mut)
#how close the solutions should be
distance_factor <- 0.4
diff = 2 #bandwidth of convergence
################################################################################
####################### EvoStatus: load or create ##############################
################################################################################
Status <- "./EvoStatus"
if(!dir.exists(Status)){
  dir.create(Status)
}
################################################################################
############################## General functions ##############################
################################################################################
# Kernel regression estimation
approx_f <- function(x, bw=1, results_smooth){

  distance <- numeric(nrow(results_smooth))

  for(index in 1:(length(param_cols))){ ### only the parameter include in cost - 6 recycle should come out
    distance <- distance + ((results_smooth[,param_cols[index]] - x[param_cols[index]])/bw[index])^2
  }

    weight <- dnorm_approx(sqrt(distance))



  target <- sum(results_smooth[,ncol(results_smooth)] * weight) / sum(weight)

  return (target)
}

# This is a faster version of dnorm()
dnorm_approx <- function(distance, bw = 2){
  exp(-distance^2/(2*(bw^2))) / bw
}

# Kernel regression density estimation
approx_density <- function(x, bw=1, results_smooth){

  distance <- numeric(nrow(results_smooth))

  for(index in 1:(length(param_cols))){
    distance <- distance + ((results_smooth[,param_cols[index]] - x[param_cols[index]])/bw[index])^2
  }

  weight <- dnorm_approx(sqrt(distance))

  target <-  sum(weight)

  return (target)
}
################################################################################
############################## Load the functions ##############################
################################################################################
# design by user
my_functions <- c("cost_function", "generate_new", "termination_criteria")

if (any(!sapply(my_functions, exists))) {
  Path <- "./Functions"
  lapply(my_functions, function(f) source(file.path(Path, paste0(f, ".R"))))
}
################################################################################
############# Read_in data from the last cycle && EvoStatus-file ###############
################################################################################
results_full <- list()
for(index in 1:length(files)){
  load(files[index])
  results_full[[index]] <- results
  #print(index)
}
results_full <- as.data.table(do.call(rbind,results_full))
cols_group <- colnames(results_full)[3:(nfactors+2)]
cols_av <- colnames(results_full)[(nfactors+3):ncol(results_full)]
results_new <- results_full[,c(.(randomSeed=randomSeed),lapply(.SD, mean)),
                            by = cols_group,.SDcols = cols_av]
setcolorder(results_new,colnames(results_full)[-1])
results_new <- as.matrix(results_new)
results_new <- results_new[,-c(1)] # remove randomseed

################################################################################
############################ Load evolutionary_log list ########################
################################################################################
if (iteration > 2) {
  load("EvoStatus/EvoStatus.RData")
  ##############################################################################
  ###################### Termination Criteria ##################################
  ##############################################################################
  # termination_criteria(evolutionary_log, iteration)
  # if you wish to create another termination criteria you can also use your own idea and just make sure you have at the end file="finish.txt" which is used by Snakemake to stop
  # if(FALSE){
  #   cat(paste0(iteration,"\n"),file="finish.txt")
  # }
} else {
  ################################################################################
  ###################### Initialize evolutionary_log list#########################
  ### Store all information related to the state of the evolutionary algorithm ###
  ################################################################################
  evolutionary_log = list()
  evolutionary_log$nfactors = config$nfactors

  if(binary_parameter){
    evolutionary_log$binary = config$binary
    number_binary_parameter = config$number_binary_parameter
    evolutionary_log$binary_freq = matrix(0, nrow = 2^sum(evolutionary_log$binary), ncol=0)
  }
  evolutionary_log$distance1 = list()
  evolutionary_log$distance2 = list()
  evolutionary_log$distance3 = list()
  evolutionary_log$distance4 = list()
  evolutionary_log$target_value = list()
  evolutionary_log$non_random = list()
  evolutionary_log$result_list <- list()
}

results_old = evolutionary_log$last_results
results_smooth <- rbind(evolutionary_log$all_results,results_new)
################################################################################
######################## Calculate new parameters to be used ###################
################################################################################
if(iteration == 2){

  redo_probability <- 0.2
  use_second_last <- FALSE

  mut_offspring <- rep(config$mut_offspring, nfactors)
  mut_parent <- rep(config$mut_parent, nfactors)

  #linked parameters that have reduced mutation rate
  if(linked_parameter){
    if(is_hybrid_breeding){
    mut_male <- param_link_male
    mut_offspring[mut_male] = mut_offspring[mut_male]/2
    mut_parent[mut_male] = mut_parent[mut_male]/2
    
    mut_female <- param_link_female
    mut_offspring[mut_female] = mut_offspring[mut_female]/2 
    mut_parent[mut_female] = mut_parent[mut_female]/2 
    }else{
      mut_line <- param_link
      mut_offspring[mut_line] = mut_offspring[mut_line]/2 
      mut_parent[mut_line] = mut_parent[mut_line]/2 
    }
  }

}  else if(iteration > 2){

  use_second_last <- TRUE

  #### only increase redo_probability if mating from mating between different binaries were not helpful ####

  if(binary_parameter){
    redo_probability = min(evolutionary_log$redo_probability[iteration - 2] + 0.1, 0.8) # never redo more than 80% but increase steadily
  }else{
    redo_probability <- config$redo_probability
  }

  mut_offspring = evolutionary_log$mut_offspring[,iteration - 2]
  mut_parent = evolutionary_log$mut_parent[,iteration - 2]

  #### when refinement is reached binary variable is fixed!
  if(binary_parameter){
    if(ncol(evolutionary_log$center)>4){
      to_test = which(evolutionary_log$binary)

      for(bven in to_test){
        fixed_center = sum(evolutionary_log$center[bven,(ncol(evolutionary_log$center)-4):ncol(evolutionary_log$center)]>0.8 |  evolutionary_log$center[bven,(ncol(evolutionary_log$center)-4):ncol(evolutionary_log$center)]<0.2)
        fixed_selected = sum(evolutionary_log$parent_center[bven,(ncol(evolutionary_log$parent_center)-4):ncol(evolutionary_log$parent_center)]>0.85 |  evolutionary_log$parent_center[bven,(ncol(evolutionary_log$parent_center)-4):ncol(evolutionary_log$parent_center)]<0.15)

        if(fixed_center==5 && fixed_selected){
          mut_offspring[bven] = 0.01
          mut_parent[bven] = 0.01
        }
      }
    }
  }
}

# Total number of offspring to generate
n_off <- n_off1 + n_off_recombination + n_off_mutation + n_off_random

if(iteration ==2){
  n_mix <- nrow(results_new)
} else{
  n_mix <- n_off1 + n_off_recombination + n_off_mutation
}

size_smooth <- nrow(results_smooth) # The smoothing function all simulations


param_cols <- (1:nfactors)
keep_parents <- TRUE

results_list <- list()
parents_lists <- list()
top10_list <- list()
topbinary_list <- list()
range_list <- list()
results_list[[1]] <- results_new

################################################################################
########################## Start of the EVO-pipeline############################
################################################################################

################################################################################
######### Sort results based on performance in the target function #############
# If binary parameters are available the results should resort based on binary #
################################################################################
cat("  Ranking the Best performing Iteration \n")

sorting <- sort(results_new[,ncol(results_new)], index.return=TRUE, decreasing = TRUE)
cand <- results_new[sorting$ix,]
# put the top 10% in the list, first candidates for the next round
top10_list[[length(top10_list)+1]] <- results_new[sorting$ix[1:(length(sorting$ix)/10)],]

# Optima in binary parameters impacts the optimum in the non-binary parameters
# Only include those simulations that were performed in settings with the ideal binary settings
# Optima means use most often in the new population
# If we have TP1.2 and TP2.2 the cost of the breeding program is larger, that's why the other cohorts are smaller in size
# if we calculate the optimum, the average performance of parents, is really in a different area of space
if(binary_parameter){
  binary_settings <- (results_new[,param_cols][,is_binary])
  # Which binary setting is available, all different binary settings which are available 2^2=4
  ubinary <- unique(binary_settings)
  #in case we have only one binary parameter change from a vector to matrix
  ubinary <- as.matrix(ubinary)
  # print(ubinary)

  # How frequently each binary setting appear, which one is the most frequent one (most promising one)
  nbinary <- numeric(nrow(ubinary))

  for(index4 in 1:nrow(ubinary)){
    temp1 <- t(binary_settings)==ubinary[index4,]
    nbinary[index4] <- sum(colSums(temp1)==sum(is_binary))
  }
  # include only the results with most common binary settings based on the average performance of only binary settings
  best_binary <- ubinary[which.max(nbinary),]
  results_binary <- results_new[colSums(t(results_new[,param_cols][,is_binary])==best_binary)==sum(is_binary),]

  sortingb <- sort(results_binary[,ncol(results_binary)], index.return=TRUE, decreasing = TRUE)
  topbinary_list[[length(topbinary_list)+1]] <- results_binary[sortingb$ix[1:(length(sortingb$ix)/10)],]
}

# calculate the sd of each parameter
range <- numeric(length(param_cols))
for(index in 1:length(param_cols)){
  #range[index] <- sd(results_new[1:(n_mix),param_cols[index]])
  range[index] <- max(0.04,sd(results_new[1:(n_mix),param_cols[index]]))
}

range_smooth <- numeric(length(param_cols))
for(index in 1:length(param_cols)){
  #range[index] <- sd(results_new[1:(n_mix),param_cols[index]])
  range_smooth[index] <- max(0.01,sd(results_new[1:(n_mix),param_cols[index]]))
}

range_list[[length(range_list)+1]] <- range
evolutionary_log$non_random[[iteration]] = 1:(n_mix)

################################################################################
##################### Select parents for the next cycle ########################
################################################################################
# The best will always be selected because there is nothing to compare against
parents_list1 <- cand[1,, drop= FALSE]
# Which of the candidates do I have to check
activ_indi <- 2
#how many were selected till far, this will be the first candidate that will be selected
selected <- 1

# The first one is already included, so I have to check the second one
which_selected <- 1

################################################################################
############ Select parents1: based on the highest performance #################
################################################################################
# How similar is this candidate to all parents that are already been selected
cat("  Start collecting Parent 1 \n")
while(selected < n_sel_high_target){
  cand[activ_indi,]

  distance <- numeric(nrow(parents_list1))
  for(index in 1:length(param_cols)){
    distance <- distance + ((parents_list1[,param_cols[index]] - cand[activ_indi, param_cols[index]])/range[index])^2
  }
  distance = sqrt(distance)

  # to make sure no similar setting will be selected
  # If one of the parents has a distance of 0.15 is high this will be added to the parent list
  # and number of candidates are adding 1 and the next parents will be checked!

  if(min(distance) > (length(param_cols) * 0.1* distance_factor)){ # relatedness is not as critical here but still filters out extremely similar lines

    # print(rbind(cand[activ_indi,], parents_list1[which.min(distance),]))
    # print(min(distance))
    parents_list1 <- rbind(parents_list1, cand[activ_indi,])
    selected <- selected + 1
    which_selected <- c(which_selected, activ_indi)
  }
  activ_indi <- activ_indi + 1
}
################################################################################
####################### Select parents3: second to last run ####################
################################################################################
cat("  Start collecting Parent 3 \n")

if(use_second_last){
  selected <- 0
  results_old <- results_old[sort(results_old[,ncol(results_old)], index.return=TRUE, decreasing = TRUE)$ix,]
  activ_indi <- 1
  parents_list3 <- NULL
  parents_sofar <- rbind(parents_list1[,param_cols])
  while(selected < n_sel_last2){
    if(sum(colSums(t(parents_sofar)==results_old[activ_indi,param_cols])==length(param_cols))==0){
      parents_list3 <- rbind(parents_list3, results_old[activ_indi,])
      selected <- selected + 1
    }
    activ_indi <- activ_indi +1
  }
} else{
  parents_list3 <- NULL
}
################################################################################
####################### Select parents2: based on smoothing ####################
################################################################################
cat("  Start collecting Parent 2 \n")

id_selected = sorting$ix[which_selected]
# We are removing the lines that we already selected
cand <- results_new[sorting$ix[-which_selected],] ## Not selected already selected individuals from the previous step
# apply smoothing on candidates, put the smooth value in the last columns, # expected performance of the candidates then will be added
# It should be smaller compared to before due to smoothing
for(index in 1:nrow(cand)){ # compute smoothing using kernel regression
  cand[index, ncol(cand)] <- approx_f(cand[index, param_cols], range_smooth/2, results_smooth)
}
sort_cand <- sort(cand[,ncol(cand)], index.return=TRUE, decreasing = TRUE)
cand <- cand[sort_cand$ix,]
if(binary_parameter){
  sortingb <- sort(results_binary[,ncol(results_binary)], index.return=TRUE, decreasing = TRUE)
}
# The same procedure as before
selected <- 1
if(n_sel_kernel>0){
  parents_list2 <- cand[1,, drop= FALSE]
  activ_indi <- 2
  while(selected < n_sel_kernel){
    cand[activ_indi,]
    consider_distance = rbind(parents_list1, parents_list3, parents_list2)
    consider_distance =  parents_list2
    distance <- numeric(nrow(consider_distance))
    for(index in 1:length(param_cols)){
      distance <- distance + ((consider_distance[,param_cols[index]] - cand[activ_indi, param_cols[index]])/range[index])^2
    }
    distance = sqrt(distance)
    # just a motivation being smoothing is already looking into the areas
    # So relatively close values will be more similar, to just have some coverage of the population we increased a bit
    #
    if(min(distance) > (length(param_cols) * 0.2 * distance_factor)){

      # print(rbind(cand[activ_indi,], parents_list1[which.min(distance),]))
      # print(min(distance))

      parents_list2 <- rbind(parents_list2, cand[activ_indi,])
      selected <- selected + 1
      id_selected =  c(id_selected, which.max(colSums(cand[activ_indi,]==t(results_new))))
    }
    activ_indi <- activ_indi + 1
  }
} else{
  parents_list2 <- NULL
}

################################################################################
######################### Collect the list of parents ##########################
################################################################################
parents_list <- rbind(parents_list1, parents_list2, parents_list3)
################################################################################
############### Prepare selected parents based on other criteria ###############
################################################################################
new_setting <- matrix(0, nrow=n_off, ncol=length(param_cols))

if(is.null(parents_list3)) {
  n_sel <- nrow(parents_list)
  new_setting[1:n_sel,] <- parents_list[,param_cols]
}else{
  new_setting[1:n_sel,] <- parents_list[,param_cols]
}

prob = rep(0.5, nfactors)

# If increasing/decreasing the parameter improves the objective function value (based on the kernel regression), the mutation probability for that parameter is increased by 0.02 (and vice versa)
# If the kernel regression suggests that high values for a parameter are helpful the change for the mutation to increase this value is increased to up to 75% (instead of 50%)
# If the kernel regression suggests that both lower and higher values are works then the overall mutation rate is reduced the parameter already seems to be in a good position
mut_offspring_temp = mut_offspring
mut_parent_temp = mut_parent
change1 = change2 = numeric(nfactors)
if(iteration > 2){
  to_check = max(1,(iteration-5)):(iteration-1)
  if(iteration>10){
    to_check[1] = iteration-9
    to_check[2] = iteration-6
  }
  # most important parameters
  nfreq = sum(!is_binary)
  for(index5 in to_check){
    x = evolutionary_log$center[,index5]
    x_value = approx_f(x, bw = range, results_smooth)
    for(index in which(!is_binary)){
      x_new = x
      x_new[index] = x_new[index] + max(range_smooth[index] , min_range_mut[index])
      x_new_value = approx_f( x_new, bw = range_smooth, results_smooth)
      x_new2 = x
      x_new2[index] = x_new2[index] - max(range_smooth[index] , min_range_mut[index])
      x_new2_value = approx_f( x_new2, bw = range_smooth, results_smooth)
      change1[index] = x_value - x_new_value
      change2[index] = x_value - x_new2_value
      if(x_new_value > x_value){
        prob[index] = prob[index] + 0.03
      }
      if(x_new2_value > x_value){
        prob[index] = prob[index] - 0.03
      }
      if(nfreq > 5){
        if((x_new_value < x_value) & (x_new2_value < x_value)){
          mut_offspring_temp[index] = mut_offspring_temp[index] * 0.75
          mut_parent_temp[index] = mut_parent_temp[index] * 0.75
        }
      }
    }
  }
  change = abs(change1 - change2)
  to_change = sort(change, index.return=TRUE, decreasing = TRUE)
  if(nfreq > 5){
    mut_offspring_temp[to_change$ix[1:round(nfreq/5)]] = mut_offspring_temp[to_change$ix[1:round(nfreq/5)]] * 2
    mut_parent_temp[to_change$ix[1:round(nfreq/5)]] = mut_parent_temp[to_change$ix[1:round(nfreq/5)]] * 2
    # less important parameters
    mut_offspring_temp[to_change$ix[round(nfreq/2+1):nfreq]] = mut_offspring_temp[to_change$ix[round(nfreq/2+1):nfreq]] / 2
    mut_parent_temp[to_change$ix[round(nfreq/2+1):nfreq]] = mut_parent_temp[to_change$ix[round(nfreq/2+1):nfreq]] / 2
  }
  # Any mutation probability greater than 0.5 is automatically set to 0.5
  mut_offspring_temp[mut_offspring_temp>0.5] = 0.5
  mut_parent_temp[mut_parent_temp>0.5] = 0.5
}
# Check mean/sd to avoid testing parameter settings that are extremely bad
# Currently this is only applied to n_off_random individuals
sd_results <- sd(results_new[1:(n_mix),ncol(results_new)]) ## column change?
mean_results <- mean(results_new[1:(n_mix),ncol(results_new)])
################################################################################
################### Generate offspring from selected parents ###################
################################################################################
cat("  Start generating offspring \n")
re <- NULL
for(index in (n_off1 + 1):n_off){
  set.seed(index)
  print(index)
  thres <- mean_results - 1*sd_results # just to avoid testing very very bad settings

  valid <- FALSE
  attempt <- 1

  # check all parameters are reasonable if not valid at some points it set to false and rerun everything again
  while(valid==FALSE){

    redo = TRUE # it is only for binary variables, if they are not available it would be set to FALSE

    while(redo){

      parents <- sample(1:n_sel, 2)
      if(index > (n_off1 + n_off_recombination)){
        parents <- rep(index%%n_sel+1,2)
        mut <- mut_parent_temp # 30%
      } else{
        mut <- mut_offspring_temp # 10
      }

      # Each parent used the same with the same frequency
      # Sampling two random parents for offspring generation
      parent1 <- parents_list[parents[1], param_cols]
      parent2 <- parents_list[parents[2], param_cols]

      if(binary_parameter==FALSE){
        redo = FALSE
      }

      if(binary_parameter){
        binary_check = prod(parent1[param_cols[is_binary]] == parent2[param_cols[is_binary]])
        if(binary_check==1){
          redo = FALSE
        } else{
          redo = as.logical(rbinom(1,1, prob = redo_probability))
        }
      }
    }
    ################################################################################
    ##################### Generate offspring with Recombination ####################
    ################################################################################
    share <- runif(1,0,1)
    off <- ((share*parent1 + (1-share)*parent2) )

    if(binary_parameter){
      for(index2 in param_cols[is_binary]){
        off[index2] <- sample(c(parent1[index2], parent2[index2]), 1)
      }


      for(index2 in param_cols[is_binary]){
        if(rbinom(1,1,mut)==1){
          off[index2] <- 1- off[index2]
        }
      }
    }
    ################################################################################
    ######################## Generate offspring with Mutation ######################
    ################################################################################
    # For continuous variable first calculate the range
    for(index2 in param_cols[!is_binary]){
      if(rbinom(1,1,mut)==1){
        ## Mutation can be by +two times of it the sd of a parameter and -two times  of it
        range_mut <- max((range_smooth[index2]*2), min_range_mut[index2])

        sign = (rbinom(1,1,prob[index2])-0.5)*2
        # max(round(range[index2]/2), 1)
        # range_mut <- max(round(range[index2]), 1)
        off[index2] <- off[index2] + sign * runif(1,0,range_mut)

        if(linked_parameter){
          if(is_hybrid_breeding){
            if(rbinom(1,1,0.5)){
              
              
              if(index2 %in% param_link_male){
                off[linked_param_male] = off[linked_param_male] + sign * sample(range_mut,1)
              }
              
              if(index2 %in% param_link_female){
                off[linked_param_female] = off[linked_param_female] + sign * sample(range_mut,1)
              }
              
            }
          }else{
            if(rbinom(1,1,0.5)){
              
              if(index2 %in% param_link){
                off[linked_param] = off[linked_param] + sign * sample(range_mut,1)
              }
              
              
            }
          }
          
        }
      }
    }
    
    off <- round(off, digits = 3)

        # Round only the columns identified by is_integer
    off[which(is_integer)] <- round(off[which(is_integer)])
    ################################################################################
    ##################### Rounding in light of the constraints ######################
    ################################################################################
    {
      base_cost_ini <- unlist(config$base_cost_ini)
      base_cost <- calc_cost(base_cost_ini) 
      new_cost <- calc_cost(off)
      
      index2 <- 1
      while(abs(new_cost-base_cost)>5000 & index2 < 20){
        off[cost_par] <- round(off[cost_par] * base_cost / new_cost)
        new_cost <- calc_cost(off)
        index2 <- index2 + 1
      }
   
      ################################################################################
      ######### Please adopt this part depending on your breeding program ############
      ########################## Adjust linked_parameters ############################
      ################################################################################

      if(linked_parameter_adjustment){
        # add the necessary adjustment
      }

      ################################################################################
      ################## Compensate the rest in the cheapest unit #####################
      ################################################################################
      if(is_hybrid_breeding){
        if(new_cost!=base_cost){
          if(config$cheapest_unit_female > length(off)){
            stop("cheapest_unit_female of the config file has to be <= ",length(off))
          }
          if(config$cheapest_unit_male > length(off)){
            stop("cheapest_unit_male of the config file has to be <= ",length(off))
          }
          
          if(rbinom(1,1,0.5)==1){
            off[cheapest_unit_male] <- off[cheapest_unit_male] - round((new_cost - base_cost) / cost_cheapest_male) -1
          } else{
            off[cheapest_unit_female] <- off[cheapest_unit_female] - round((new_cost - base_cost) / cost_cheapest_female) -1
          }
        }
      }
      if(is_line_breeding){
        if(new_cost!=base_cost){
          if(config$cheapest_unit > length(off)){
            stop("cheapest_unit of the config file has to be <= ",length(off))
          }
          
          off[cheapest_unit] <- off[cheapest_unit] - round((new_cost - base_cost) / cost_cheapest) -1
          
        }
      }
    }
    ################################################################################
    ################################################################################


    ################################################################################
    ####################### Managing diversity of offspring  #######################
    ################################################################################
    activ_indi <- 2
    selected <- 1

    if(index > (n_off1 + n_off_recombination) || index < (n_off1+1)){
      distance <- Inf
    } else{

      #
      distance <- numeric(nrow(new_setting))
      for(index2 in 1:length(param_cols)){
        distance <- distance + ((new_setting[,index2] - off[index2])/range[index2])^2
      }
      distance = sqrt(distance)

      re <- c(re, min(distance))
    }

    # Test some random settings
    # This should avoid long-term convergence into a local maxima instead of a global maxima
    if(index > (n_off1 + n_off_recombination + n_off_mutation)){

      est_value <- 0
      tt <- 1

      while(est_value < thres){
        off <- generate_new()


        if(size_smooth<Inf){
          est_value <- approx_f(off, bw = range_smooth, results_smooth)
        } else{
          est_value <- approx_f(off, bw = range_smooth, results_smooth)
        }


        if(tt%%100==0){
          thres <- thres - sd_results
        }

        tt <- tt + 1

        if(is.na(est_value)){
          est_value = thres # no points in the area - investigate
        }
      }
    }


    is_copy <- sum(colSums(t(new_setting)==off)==nfactors)>0
    # is_cor <- (sum(off[3:8] > 0)==6) & (off[3]>=off[5]) & (off[4]>=off[6]) & (off[5]>=off[7]) & (off[6]>=off[8])

    is_cor <- (sum(off[which(!is_binary)] > 0)==nfactors-number_binary_parameter)




    #    if(((sum(off>0)==3 & off[2]<=700 & off[2]>=100) & !is_copy)||attempt > 100 & min(distance)>0.1){
    if((!is_copy  && min(distance)>(0.05 * nfactors * distance_factor) && is_cor)||attempt > 100){
      valid <- TRUE
      attempt <- 1
    } else{
      attempt <- attempt + 1
    }

    ################################################################################
    ######### Please adopt this part depending on your breeding program ############
    ############################## keeping constraints #############################
    ################################################################################
    
# if (some constraints that should be hold) {
#   # print("All conditions for contrains are satisfied.")
#   valid <- TRUE
#   attempt <- 1
# } else {
#   valid <- FALSE
#   attempt <- attempt + 1
#   # print("Not all conditions for contrains are satisfied.")
#   
# }

    ################################################################################

    
    
  }

  new_setting[index,] <- off

}


new_setting <- cbind(sample(1:1e5,1),new_setting)

if(any(rowSums(new_setting[,-1])== 0)){
  warning("some new settings are zero --> check")
  new_setting <- new_setting[rowSums(new_setting[,-1]) != 0,]
}

#add the column's name
names <- config$name_parameter
colnames(new_setting) <- c("randomSeed", names)
rownames(new_setting) <- NULL
print(new_setting)
write.csv(new_setting, file = output)
################################################################################
#### Putting together information from the evolutionary pipeline performed #####
################################################################################

evolutionary_log$selection.rate = cbind(evolutionary_log$selection.rate ,c(n_sel_high_target, n_sel_kernel, n_sel_last2))
evolutionary_log$n.iteration = iteration

evolutionary_log$n_off = cbind(evolutionary_log$n_off ,c(n_off1, n_off_recombination, n_off_mutation, n_off_random))

evolutionary_log$redo_probability = c(evolutionary_log$redo_probability, redo_probability)
evolutionary_log$mut_offspring = cbind(evolutionary_log$mut_offspring, mut_offspring)
evolutionary_log$mut_parent = cbind(evolutionary_log$mut_parent, mut_parent)
evolutionary_log$use_second_last = c(evolutionary_log$use_second_last , use_second_last)

# evolutionary_log$n_repetitions = c(evolutionary_log$n_repetitions, n_repetitions) # number of times each setting is simulated (search / refinement phase)
evolutionary_log$n_rep = c(evolutionary_log$n_rep, n_rep)
evolutionary_log$last_results =  results_new # This is were the last cycle is stored
evolutionary_log$all_results = rbind(evolutionary_log$all_results ,
                                     results_new) # This is where all prior cycles are stored

# evolutionary_log$nrep = c(evolutionary_log$nrep, nrep)

if(evolutionary_log$n_off[4, iteration-1] > 0){
  new_setting = new_setting[1: sum(evolutionary_log$n_off[1:3, iteration-1]), ]
}

if(iteration > 2 && evolutionary_log$n_off[4, iteration-2] > 0){
  results_new = results_new[1: sum(evolutionary_log$n_off[1:3, iteration-2]),]
}
if(iteration > 2 && evolutionary_log$n_off[4, iteration-2] > 0){
  results_old = results_old[1: sum(evolutionary_log$n_off[1:3, iteration-3]),]
}

if(iteration==2){
  evolutionary_log$center = cbind(colMeans(results_new[,1:nfactors]))
}
evolutionary_log$center = cbind(evolutionary_log$center, colMeans(new_setting[,-1]))
evolutionary_log$parent_center = cbind(evolutionary_log$parent_center, colMeans(parents_list))

evolutionary_log$sd = cbind(evolutionary_log$sd, sqrt(diag(var( rbind(results_new[,1:nfactors], results_old[,1:nfactors], new_setting[,c(-1)])))))


distance1 = matrix(0, ncol = iteration, nrow=iteration) # l1 absolute distance
## More emphasis on small changes
distance2 = matrix(0, ncol = iteration, nrow=iteration) # l2 euclidean distance
# motivation of Euclidean distance is that when we have a large change in one direction and a small change in another direction the smaller change is not important
for(index in 1:iteration){
  for(index2 in 1:iteration){
    distance1[index,index2] = (sum((abs(evolutionary_log$center[,index] - evolutionary_log$center[,index2]) / evolutionary_log$sd[,iteration-1])))
    distance2[index,index2] = sqrt(sum((abs(evolutionary_log$center[,index] - evolutionary_log$center[,index2]) / evolutionary_log$sd[,iteration-1])^2))
  }
}

distance3 = matrix(0, ncol = iteration-1, nrow=iteration-1) # number of parameters that moved in the same direction

for(index in 1:(iteration-1)){
  for(index2 in 1:(iteration-1)){
    dir1 = (evolutionary_log$center[,index] - evolutionary_log$center[,index+1]) > 0
    dir2 = (evolutionary_log$center[,index2] - evolutionary_log$center[,index2+1]) > 0
    distance3[index,index2] = sum((dir1&dir2) | (!dir1&!dir2) )
  }

}

distance4 = matrix(0, ncol = iteration-1, nrow=iteration-1) # angel

for(index in 1:(iteration-1)){
  for(index2 in 1:(iteration-1)){
    dir1 = ((evolutionary_log$center[,index] - evolutionary_log$center[,index+1])) / evolutionary_log$sd[,iteration-1]
    dir2 = ((evolutionary_log$center[,index2] - evolutionary_log$center[,index2+1])) /  evolutionary_log$sd[,iteration-1]
    suppressWarnings({distance4[index,index2] = acos(sum(dir1 * dir2) / sqrt(sum(dir1^2)) / sqrt(sum(dir2^2)))})
  }

}

evolutionary_log$distance1[[iteration-1]] = distance1
evolutionary_log$distance2[[iteration-1]] = distance2
evolutionary_log$distance3[[iteration-1]] = distance3
evolutionary_log$distance4[[iteration-1]] = distance4

if(iteration>2){

  prop_noff1 = prop_noff2 = prop_noff3 = prop_noff4 = 0
  if(evolutionary_log$n_off[1,iteration-2] > 0){
    prop_noff1 = sum(id_selected <= evolutionary_log$n_off[1,iteration-2]) / evolutionary_log$n_off[1,iteration-2]
  }

  if(evolutionary_log$n_off[2,iteration-2] > 0){
    prop_noff2 = sum((id_selected <= sum(evolutionary_log$n_off[1:2,iteration-2])) & (id_selected > sum(evolutionary_log$n_off[1,iteration-2])) ) / evolutionary_log$n_off[2,iteration-2]
  }

  if(evolutionary_log$n_off[3,iteration-2] > 0){
    prop_noff3 = sum((id_selected <= sum(evolutionary_log$n_off[1:3,iteration-2])) & (id_selected > sum(evolutionary_log$n_off[1:2,iteration-2])) ) / evolutionary_log$n_off[3,iteration-2]
  }

  if(evolutionary_log$n_off[4,iteration-2] > 0){
    prop_noff4 = sum((id_selected <= sum(evolutionary_log$n_off[1:4,iteration-2])) & (id_selected > sum(evolutionary_log$n_off[1:3,iteration-2])) ) / evolutionary_log$n_off[4,iteration-2]
  }
  evolutionary_log$proportion_sel = cbind(evolutionary_log$proportion_sel, c(prop_noff1, prop_noff2, prop_noff3, prop_noff4))
}


evolutionary_log$target_value[[iteration-1]] = results_new[,ncol(results_new)]

evolutionary_log$target_value_avg = c(evolutionary_log$target_value_avg, mean(results_new[,ncol(results_new)]))
evolutionary_log$target_value_b10= c(evolutionary_log$target_value_b10, quantile(results_new[,ncol(results_new)], probs = 0.9))


if(binary_parameter){

  evolutionary_log$est_pos = cbind(evolutionary_log$est_pos, colMeans(topbinary_list[[length(topbinary_list)]]))
  evolutionary_log$est_approx = c(evolutionary_log$est_approx, approx_f(colMeans(topbinary_list[[length(topbinary_list)]]), range_smooth, results_smooth)) # this is possible in practise
  avg_approx = numeric(nrow(topbinary_list[[length(topbinary_list)]]))
  for(index in 1:length(avg_approx)){
    avg_approx[index] = approx_f(topbinary_list[[length(topbinary_list)]][index,], range_smooth/2, results_smooth)
  }
  evolutionary_log$est_approx_max = c(evolutionary_log$est_approx_max, max(avg_approx)) #
  evolutionary_log$est_approx_sd = c(evolutionary_log$est_approx_sd, sd(avg_approx)) #
  evolutionary_log$est_approx_avg = c(evolutionary_log$est_approx_avg, mean(avg_approx)) #

  evolutionary_log$est_value = c(evolutionary_log$est_value, eval(rbind(colMeans(topbinary_list[[length(topbinary_list)]])))) # this is possible in practise

  avg_approx = numeric(nrow(topbinary_list[[length(topbinary_list)]]))
  for(index in 1:length(avg_approx)){
    avg_approx[index] = eval(rbind(topbinary_list[[length(topbinary_list)]][index,]))
  }

  evolutionary_log$est_value_max = c(evolutionary_log$est_value_max, max(avg_approx)) #
  evolutionary_log$est_value_sd = c(evolutionary_log$est_value_sd, sd(avg_approx)) #
  evolutionary_log$est_value_avg = c(evolutionary_log$est_value_avg, mean(avg_approx)) #

  to_test = new_setting[,-1]
  to_test <- to_test[colSums(t(to_test[,param_cols][,is_binary])==best_binary)==sum(is_binary),]
}else{

  evolutionary_log$est_pos = cbind(evolutionary_log$est_pos, colMeans(top10_list[[length(top10_list)]]))
  evolutionary_log$est_approx = c(evolutionary_log$est_approx, approx_f(colMeans(top10_list[[length(top10_list)]]), range_smooth, results_smooth)) # this is possible in practise
  avg_approx = numeric(nrow(top10_list[[length(top10_list)]]))
  for(index in 1:length(avg_approx)){
    avg_approx[index] = approx_f(top10_list[[length(top10_list)]][index,], range_smooth/2, results_smooth)
  }
  evolutionary_log$est_approx_max = c(evolutionary_log$est_approx_max, max(avg_approx)) #
  evolutionary_log$est_approx_sd = c(evolutionary_log$est_approx_sd, sd(avg_approx)) #
  evolutionary_log$est_approx_avg = c(evolutionary_log$est_approx_avg, mean(avg_approx)) #

  evolutionary_log$est_value = c(evolutionary_log$est_value, eval(rbind(colMeans(top10_list[[length(top10_list)]])))) # this is possible in practise

  avg_approx = numeric(nrow(top10_list[[length(top10_list)]]))
  for(index in 1:length(avg_approx)){
    avg_approx[index] = eval(rbind(top10_list[[length(top10_list)]][index,]))
  }

  evolutionary_log$est_value_max = c(evolutionary_log$est_value_max, max(avg_approx)) #
  evolutionary_log$est_value_sd = c(evolutionary_log$est_value_sd, sd(avg_approx)) #
  evolutionary_log$est_value_avg = c(evolutionary_log$est_value_avg, mean(avg_approx)) #
}


evolutionary_log$prob = cbind(evolutionary_log$prob, prob)
evolutionary_log$mut_temp = cbind(evolutionary_log$mut_temp, mut_offspring_temp)
###########################################################################################
################# Parameter convergence using kernel regression ###########################
###########################################################################################
cat('Parameter convergence using kernel regression')
if(iteration >=2){


  for (it in (iteration - 1)) {

    print(it)

    files <- paste0(paste0("results/iteration_", it, "/"), dir(paste0("results/iteration_", it)))
    files = files[substr(files, start = nchar(files)-4, stop = nchar(files))=="RData"]
    
    result <- list()
    for (i in 1:length(files)) {
      load(files[i])
      result[[i]] <- cbind(
        as.integer(str_extract(files[i], "[[:digit:]]+")),
        as.integer(str_remove(str_extract(files[i], "results_[[:digit:]]+"), "results_")),
        results
      )
    }

    results <- as.data.table(do.call(rbind, result))
    results <- results[,c(.(rep=rep),lapply(.SD, mean)),
                                by = cols_group,.SDcols = cols_av]
    
    rep_del= nfactors+1
    results <- results[, -rep_del, with = FALSE]
    evolutionary_log$result_list[[it]] <- results
  }

  param_cols <- 1:nfactors
  

    i = iteration
    print(it)
    include = 1:it
    include2 = max(1,it-5): it
    results_smooth = NULL
    for(index in include){

      results_smooth = rbind(results_smooth, evolutionary_log$result_list[[index]])
      results_smooth = as.matrix(results_smooth)
    }

    # only consider non-randomly sampled settings ((particularly for range calculation!))
    results_smooth2 = NULL
    for(index in include2){
      results_smooth2= rbind(results_smooth2, evolutionary_log$result_list[[index]][evolutionary_log$non_random[[iteration]],])
      results_smooth2 = as.matrix(results_smooth2)
    }

    range <- numeric(length(param_cols))
    for(index in 1:length(param_cols)){
      range[index] <- max(0.00001,sd(results_smooth2[,param_cols[index]]))
    }


    dens = numeric(nrow(results_smooth2))
    for(index in 1:nrow(results_smooth2)){
      dens[index] = approx_density(results_smooth2[index,], range/2, results_smooth2)
    }

    # only areas that in the last few iterations are tested frequently are considered for the global optima
    # Avoid picking an outlier in an area with very few simulations
    consider = dens > quantile(dens, probs = 0.2)
    
    
    if(length(target_density)>0){
      diff = 2
      check = TRUE
      while(check){
        diff = diff + 0.2
        # just check edge cases to as this otherwise can take some time
        to_test = (which(dens > quantile(dens, probs = 0.2) & dens < quantile(dens, probs = 0.25)))
        dens_test = numeric(nrow(results_smooth2))
        for(index in to_test){
          dens_test[index] = approx_density(results_smooth2[index,], range/diff, results_smooth)
        }
        if(min(dens_test[to_test]) > target_density){
          
        } else{
          check = FALSE
          diff = diff - 0.2
        }
      }
    } else{
      diff = 2
    }
    
    
    
    # evaluate the candidates based on ALL simulations and not only the last few cycles
    target = numeric(nrow(results_smooth2))
    for(index in 1:nrow(results_smooth2)){
      target[index] = approx_f(results_smooth2[index,], range/diff, results_smooth)
    }
    
    target[!consider] = -Inf
    convergence =  c(results_smooth2[which.max(target),], target[which.max(target)])
    
    evolutionary_log$convergence = cbind(evolutionary_log$convergence, convergence)
}

evo_target = numeric(ncol(evolutionary_log$convergence))
for(index in 1:ncol(evolutionary_log$convergence)){
  evo_target[index] = approx_f(evolutionary_log$convergence[,index], range/diff, results_smooth)
}
evolutionary_log$convergence_target = evo_target


save(file = "EvoStatus/EvoStatus.RData", list = c("evolutionary_log"))
save(list = c("evolutionary_log"), file = paste0("EvoStatus/EvoStatus", iteration, ".RData"))

