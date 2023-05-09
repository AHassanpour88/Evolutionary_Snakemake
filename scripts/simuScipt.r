args <- commandArgs(TRUE)
# args <- c(1,1,22122,0,1,30609,35223,2553,1930,247,317,16,39,421,348,144,93,114,35,12,217)
print(args)
output <- args[1]
rep <- as.integer(args[2])
randomSeed <- as.integer(args[3]) * rep
HT1_2 <- as.numeric(args[4])
HT2_2 <- as.numeric(args[5])
DH0_capacity_female <- as.numeric(args[6])
DH0_capacity_male <- as.numeric(args[7])
P0_female <- as.numeric(args[8])
P0_male <- as.numeric(args[9])
P1_female <- as.numeric(args[10])
P1_male <- as.numeric(args[11])
P2_female <- as.numeric(args[12])
P2_male <- as.numeric(args[13])
ncross_female <- as.numeric(args[14])
ncross_male <- as.numeric(args[15])
female_rec_p2 <- as.numeric(args[16])
female_rec_p1 <- as.numeric(args[17])
female_rec_p0 <- as.numeric(args[18])
male_rec_p2 <- as.numeric(args[19])
male_rec_p1 <-  as.numeric(args[20])
male_rec_p0 <- as.numeric(args[21])
results <- NULL


inbreeding_target = runif(1,0,1)
target_function = 100000 + 500 * HT1_2 - 100* HT2_2 -
  abs(male_rec_p2 - 100)  - abs(male_rec_p1 - 30)  - abs(male_rec_p0 - 110)  -
  abs(female_rec_p2 - 250) - abs(female_rec_p1 - 35)  - abs(female_rec_p0 - 35)  -
  abs(DH0_capacity_female - 34000) * 0.05 - abs(DH0_capacity_male - 27000) * 0.05 -
  abs(ncross_female - 400) - abs(ncross_male - 300) + rnorm(1, 0, 250)

results <- t(c(rep,
               randomSeed,
               HT1_2,
               HT2_2,
               DH0_capacity_female,
               DH0_capacity_male,
               P0_female,
               P0_male,
               P1_female,
               P1_male,
               P2_female,
               P2_male,
               ncross_female,
               ncross_male,
               female_rec_p2,
               female_rec_p1,
               female_rec_p0,
               male_rec_p2,
               male_rec_p1,
               male_rec_p0,
               inbreeding_target,
               target_function))

colnames(results) <- c("rep",
                       "randomSeed",
                       "HT1_2",
                       "HT2_2",
                       "DH0_capacity_female",
                       "DH0_capacity_male",
                       "P0_female",
                       "P0_male",
                       "P1_female",
                       "P1_male",
                       "P2_female",
                       "P2_male",
                       "ncross_female",
                       "ncross_male",
                       "female_rec_p2",
                       "female_rec_p1",
                       "female_rec_p0",
                       "male_rec_p2",
                       "male_rec_p1",
                       "male_rec_p0",
                       "inbreeding_target",
                       "target_function")

save(list=c("results"), file = output)




