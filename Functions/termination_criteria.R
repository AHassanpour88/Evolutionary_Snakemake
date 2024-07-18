#################################################################################
################################# Agreement #####################################
#################################################################################
# Authors
# Azadeh Hassanpour, azadeh.hassanpour@uni-goettingen.de
# Johannes Geibel, johannes.geibel@fli.de
# Torsten Pook, Torsten.pook@wur.nl
# Copyright © 2020 – 2024
# This program falls under a NonCommercial-NoDerivates-NoDistriubtion Public License.
# With use of these scripts, I confirm that I represent an academic institute and acknowledge that I shall use this script only for research. I explicitly acknowledge the terms in the license agreement https://github.com/AHassanpour88/Evolutionary_Snakemake/blob/main/License.md. I understood that any commercial use needs a commercial license from the owner of the script. For more information about a commercial license please contact Torsten Pook.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# 
#################################################################################
#################################################################################
#################################################################################



## this function can be changed to any termination criteria the user wants to implement
# the general termination criteria is when the number of final iteration is reached by the algorithm
# For the example provided we did not set any criteria and we visually observed if the target function is changing or not (this is more efficient because often you do not know what level of threshold is satisfactory!)

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

