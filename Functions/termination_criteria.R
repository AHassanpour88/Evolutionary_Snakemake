## this function can be changed to any termination criteria the user wants to implement
# the general termination criteria is when the number of final iteration is reached by the algorthm 

termination_criteria <- function(evolutionary_log, 
                                iteration, 
                                thresh_target_ave = 0.015, 
                                thresh_target_sd = 0.02, 
                                generation_back = 5) {
  
  num_generations <- iteration
  to_check <- (num_generations - generation_back):num_generations
  
  sd_condition_met <- FALSE
  avg_change_condition_met <- FALSE
  
  if (num_generations > 10) {
    
    last_sd_values <- evolutionary_log$est_approx_sd[to_check]
    last_avg_values <- evolutionary_log$est_approx_avg[to_check]
    
    sd_condition_met <- all(last_sd_values < thresh_target_sd)
    avg_change_condition_met <- abs(diff(last_avg_values)) < thresh_target_ave
    
  }
  if (sd_condition_met && all(avg_change_condition_met == TRUE)) {
    cat(paste("Both conditions are met at iteration", iteration, "\n"))
    cat(paste0(iteration, "\n"), file = "finish.txt")
  } else {
    cat(paste("Conditions are not met at iteration", iteration, "\n"))
  }
}

# Example usage:
# Assuming evolutionary_log, iteration are defined before calling the function
# terminate_evolution(evolutionary_log, iteration)

