## define the cost of your breeding programs based on different parameters
# Baseline
# x = c(TRUE, 100, 100, 10000, 10000, 500, 500, 300, 300, 30, 30)

calc_cost <- function(x){
  
  cost <-
    # nCrosses 
    x[2] * 30 + # 100 * 30
    x[3] * 30 + 
    # nDH -> size nDH female
    x[4] * 30 + 
    x[5] * 30 + 
    # nDH -> MAS HDRW female
    # 10/HDRW with cost of genotyping of 5 in 1 location
    x[4] * 10 * 5 * 1 +  
    x[5] * 10 * 5 * 1 + 
    # OBS1 -> GS 
    # 20/plot with cost of genotyping of 15 in 2 location 
    x[6] * 20 * 15 * 2 + 
    x[7] * 20 * 15 * 2 + 
    # TC1 production 
    x[6] * 20 + 
    x[7] * 20 + 
    # OBS2 
    # 20/plot in 2 location 
    x[8] * 20 * 2 + 
    x[9] * 20 * 2 + 
    # TC1 YT  
    # 50/plot in 2 location
    x[8] * 50 * 2 + 
    x[9] * 50 * 2 + 
    # TC2 production 
    # 50/plot 
    x[8] * 50 + 
    x[9] * 50 + 
    # OBS3 
    # 20/plot in 2 location 
    x[10] * 20 * 2 + 
    x[11] * 20 * 2 + 
    # TC2 YT  
    # 50/plot in 3 location
    x[10] * 50 * 3 + 
    x[11] * 50 * 3 + 
    # Hybrid prod 
    # half of the lines from OBS3 will be tested (baseline 15 out of 30)
    # percentage of share of lines go from OBS3.1 to OBS3.2
    15 * 100 * 0.5 + # decision of 0.5
    15 * 100 * 0.5 + # decision of 0.5
    # TC2 YT replication on female side
    # 50/plot in 3 location
    # half of the lines from OBS3 will be tested (baseline 20 out of 30)
    # 70% of share of lines go from OBS3.1 to OBS3.2
    ifelse(x[1] == 0, 0, (x[10] * 0.7 * 50 * 3 * x[1])) + 
    # Hybrid YT 
    # 50/plot in 10 location
    # half of the lines from OBS3 will be tested (baseline 15 out of 30)
    15 * 50 * 10 * 0.5 +
    15 * 50 * 10 * 0.5 
  return(unname((cost)))
}



