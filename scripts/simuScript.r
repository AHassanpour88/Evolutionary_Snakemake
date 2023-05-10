args <- commandArgs(TRUE)
print(args)
output <- args[1]
rep <- as.integer(args[2])
randomSeed <- as.integer(args[3]) * rep
n_test <- as.integer(args[4])
n_bull <- as.integer(args[5])
n_bull_sel <- as.integer(args[6])
results <- NULL

library(miraculix)
library(RandomFieldsUtils)
library(MoBPS)

RandomFieldsUtils::RFoptions(cores = 1)

print(RFoptions()$basic$cores)

n_simulation = 1
for(setting in 1:n_simulation){
  cat("Random seed: ", randomSeed, "\n")
  set.seed(randomSeed)
  
  
  h2 <- 0.3 # heritability
  n_cow <- 3000 
  
  
  
  {
    cat("MoBPS is starting")
    population <- creating.diploid(map=MoBPSmaps::map_cattle3,
                                   nindi= 100,
                                   n.additive = 1000,
                                   mean.target =100,
                                   var.target = 10, miraculix=TRUE)
    
    population <- breeding.diploid(population,
                                   heritability = h2,
                                   display.progress = FALSE)
    
    population <- breeding.diploid(population,
                                   breeding.size = c(n_bull, 0), ##hier
                                   name.cohort = "Testbullen_0",
                                   display.progress = FALSE)
    
    population <- breeding.diploid(population,
                                   breeding.size = c(n_bull, 0), ##hier
                                   name.cohort = "Testbullen_-1",
                                   display.progress = FALSE)
    
    population <- breeding.diploid(population,
                                   breeding.size = c(n_bull_sel, 0), ##hier
                                   name.cohort = "Zuchtbullen_0",
                                   display.progress = FALSE)
    
    population <- breeding.diploid(population,
                                   breeding.size = c(0,n_cow),
                                   name.cohort = "Kuehe_0",
                                   display.progress = FALSE)
    
    population <- breeding.diploid(population,
                                   breeding.size = c(0,n_cow),
                                   name.cohort = "Kuehe_-1",
                                   display.progress = FALSE)
    
    population <- breeding.diploid(population,
                                   breeding.size = c(0,n_cow),
                                   name.cohort = "Kuehe_-2",
                                   display.progress = FALSE)
  }
  
  # 5 Generations of Burn-in to initialize the breeding program
  for(index in 1:15){
    
    population <- breeding.diploid(population,
                                   breeding.size = c(n_bull, 0), #hier
                                   selection.m.cohorts = paste0("Zuchtbullen_", index-1),
                                   selection.f.cohorts = paste0("Kuehe_", index-1),
                                   max.offspring = c(Inf,1),
                                   name.cohort = paste0("Testbullen_", index),
                                   display.progress = FALSE)
    
    population <- breeding.diploid(population,
                                   breeding.size = c(0,n_test),
                                   selection.m.cohorts = paste0("Testbullen_", index-1),
                                   selection.f.cohorts = paste0("Kuehe_", index-1),
                                   # max.offspring = c(ceiling(n_test/n_bull),1),
                                   max.offspring = c(Inf,1),
                                   name.cohort = paste0("TestToechter_", index),
                                   display.progress = FALSE)
    
    if(index>1){ #phenoype the TestToechter if they are above one year old
      population <- breeding.diploid(population,
                                     phenotyping.cohorts = paste0("TestToechter_", index-1))
    }
    
    if(index>2){
      population <- breeding.diploid(population,
                                     offspring.bve.parents.cohorts = paste0("Testbullen_", index-2), #Cohorts to consider to derive phenotype from offspring phenotypes
                                     bve = TRUE,
                                     bve.cohorts = paste0("Testbullen_", index-2),
                                     input.phenotype = "off") #Select what to use in BVE, offspring phenotype ("off")
    }
    
    population <- breeding.diploid(population, breeding.size = c(n_bull_sel,0), #hier
                                   selection.size = c(n_bull_sel,0), #hier
                                   copy.individual.m = TRUE,
                                   selection.criteria = "bve",
                                   selection.m.cohorts =  paste0("Testbullen_", index-2),
                                   name.cohort = paste0("Zuchtbullen_",index),
                                   display.progress = FALSE)
    
    if(index>1){
      population <- breeding.diploid(population,
                                     bve=TRUE,
                                     relationship.matrix = "kinship",#Method to calculate relationship matrix for the breeding value estimation
                                     bve.cohorts = paste0("Kuehe_", c(index-1, index-2, index-3)))
    }
    
    
    population <- breeding.diploid(population, breeding.size = c(0, n_cow),
                                   selection.size = c(n_bull_sel, n_cow*2),
                                   selection.criteria = c("bve"),
                                   selection.m.cohorts = paste0("Zuchtbullen_", index-1),
                                   selection.f.cohorts = paste0("Kuehe_", c(index-1, index-2, index-3)),
                                   max.offspring = c(Inf,1),
                                   name.cohort=paste0("Kuehe_", index),
                                   display.progress = FALSE)
    
    population <- breeding.diploid(population,
                                   phenotyping.cohorts = paste0("Kuehe_", index-1))
    
    if(index==5){
      population <- bv.standardization(population,  #Function to get mean and genetic variance of a trait to a fixed value
                                       cohorts=paste0("Kuehe_", index),
                                       adapt.bve = TRUE, #Modify previous breeding value estimations by scaling
                                       adapt.pheno = TRUE, #Modify previous phenotypes by scaling
                                       mean.target = 100,
                                       var.target = 10)
    }
  }
  
  acc1 <- acc2 <- 0
  for(index3 in 4:13){
    acc1 <- acc1 + cor(get.bv(population, cohorts=paste0("Testbullen_",index3))[1,], get.pheno.off(population, cohorts=paste0("Testbullen_",index3))[1,], use="complete.obs")
    acc2 <- acc2 + cor(get.bv(population, cohorts=paste0("Testbullen_",index3))[1,], get.bve(population, cohorts=paste0("Testbullen_",index3))[1,], use="complete.obs")
  }
  bv <- mean(get.bv(population, cohorts=paste0("Kuehe_", 15)))
  kin <- rbind(kinship.emp.fast(population=population, cohorts=paste0("Kuehe_", c(15, 14, 13))))
  target_function <- bv - 50 * kin[1]
  results <- t(c(rep, randomSeed, n_test, n_bull, n_bull_sel, bv, kin[1], acc1, acc2, target_function))
  colnames(results) <- c("rep", "randomSeed", "n_test", "n_bull", "n_bull_sel", "bv", "kin", "acc1", "acc2", "target_function")
  print(results)
  save(list=c("results"), file = output)
  
}


# args <- commandArgs(TRUE)
# # args <- c(1,1,22122,0,1,30609,35223,2553,1930,247,317,16,39,421,348,144,93,114,35,12,217)
# print(args)
# output <- args[1]
# rep <- as.integer(args[2])
# randomSeed <- as.integer(args[3]) * rep
# HT1_2 <- as.numeric(args[4])
# HT2_2 <- as.numeric(args[5])
# DH0_capacity_female <- as.numeric(args[6])
# DH0_capacity_male <- as.numeric(args[7])
# P0_female <- as.numeric(args[8])
# P0_male <- as.numeric(args[9])
# P1_female <- as.numeric(args[10])
# P1_male <- as.numeric(args[11])
# P2_female <- as.numeric(args[12])
# P2_male <- as.numeric(args[13])
# ncross_female <- as.numeric(args[14])
# ncross_male <- as.numeric(args[15])
# female_rec_p2 <- as.numeric(args[16])
# female_rec_p1 <- as.numeric(args[17])
# female_rec_p0 <- as.numeric(args[18])
# male_rec_p2 <- as.numeric(args[19])
# male_rec_p1 <-  as.numeric(args[20])
# male_rec_p0 <- as.numeric(args[21])
# results <- NULL
# 
# 
# inbreeding_target = runif(1,0,1)
# target_function = 100000 + 500 * HT1_2 - 100* HT2_2 -
#   abs(male_rec_p2 - 100)  - abs(male_rec_p1 - 30)  - abs(male_rec_p0 - 110)  -
#   abs(female_rec_p2 - 250) - abs(female_rec_p1 - 35)  - abs(female_rec_p0 - 35)  -
#   abs(DH0_capacity_female - 34000) * 0.05 - abs(DH0_capacity_male - 27000) * 0.05 -
#   abs(ncross_female - 400) - abs(ncross_male - 300) + rnorm(1, 0, 250)
# 
# results <- t(c(rep,
#                randomSeed,
#                HT1_2,
#                HT2_2,
#                DH0_capacity_female,
#                DH0_capacity_male,
#                P0_female,
#                P0_male,
#                P1_female,
#                P1_male,
#                P2_female,
#                P2_male,
#                ncross_female,
#                ncross_male,
#                female_rec_p2,
#                female_rec_p1,
#                female_rec_p0,
#                male_rec_p2,
#                male_rec_p1,
#                male_rec_p0,
#                inbreeding_target,
#                target_function))
# 
# colnames(results) <- c("rep",
#                        "randomSeed",
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
#                        "male_rec_p0",
#                        "inbreeding_target",
#                        "target_function")
# 
# save(list=c("results"), file = output)
# 