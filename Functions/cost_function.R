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

## define the cost of your breeding programs based on different parameters
calc_cost <- function(x){
 cost <-
   x[1] * (4000 + 1000*x[4]) +
   x[2] * 3000
  return(unname((cost)))
}

