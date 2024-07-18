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

generate_new <- function(){
  ### this can be as exact as the initial bound in sampleScript or new bounds
  ### As initial bounds exclude the area of optimum we add them here to check a broader range
  ### this will be activated only if n_off_random in iteration.csv in config folder will be > 0 
  x <- numeric(4)
  x[4] <- rbinom(1,1,0.5)
  x[1] <- sample(100:700,1)
  x[2] <- floor((10000000-3000*x[1])/(4000 + 1000*x[4]))
  x[3] <- sample(3:30, 1)

  x <- c(x[2],x[1],x[3],x[4])
  
  return(x)
}
