args <- commandArgs(TRUE)
output <- args[1]
rep <- as.integer(args[2])
randomSeed <- as.integer(args[3]) * rep
TC2.2_YT <- as.numeric(args[4])
ncross_female <- as.numeric(args[5]) 
ncross_male <- as.numeric(args[6]) 
DH0_capacity_female <- as.numeric(args[7]) 
DH0_capacity_male <- as.numeric(args[8])
OBS1_female <- as.numeric(args[9]) 
OBS1_male <- as.numeric(args[10])
OBS2_female <- as.numeric(args[11])
OBS2_male <- as.numeric(args[12]) 
OBS3_female <- as.numeric(args[13]) 
OBS3_male <- as.numeric(args[14]) 
female_rec_OBS3 <- as.numeric(args[15]) 
female_rec_OBS2 <- as.numeric(args[16]) 
female_rec_OBS1 <- as.numeric(args[17]) 
male_rec_OBS3 <- as.numeric(args[18])
male_rec_OBS2 <-  as.numeric(args[19]) 
male_rec_OBS1 <- as.numeric(args[20]) 
results <- NULL

cat("rep: ", rep, "\n")
cat("Random seed: ", randomSeed, "\n")
cat("TC2.2_YT: ", TC2.2_YT, "\n")
cat("ncross_female: ", ncross_female, "\n")
cat("ncross_male: ", ncross_male, "\n")
cat("DH0_capacity_female: ", DH0_capacity_female, "\n")
cat("DH0_capacity_male: ", DH0_capacity_male, "\n")
cat("OBS1_female: ", OBS1_female, "\n")
cat("OBS1_male: ", OBS1_male, "\n")
cat("OBS2_female: ", OBS2_female, "\n")
cat("OBS2_male: ", OBS2_male, "\n")
cat("OBS3_female: ", OBS3_female, "\n")
cat("OBS3_male: ", OBS3_male, "\n")
cat("female_rec_OBS3: ", female_rec_OBS3, "\n")
cat("female_rec_OBS2: ", female_rec_OBS2, "\n")
cat("female_rec_OBS1: ", female_rec_OBS1, "\n")
cat("male_rec_OBS3: ", male_rec_OBS3, "\n")
cat("male_rec_OBS2: ", male_rec_OBS2, "\n")
cat("male_rec_OBS1: ", male_rec_OBS1, "\n")

# args =c(1,1 ,34060 ,1 ,0, 149, 72, 11203, 8863, 574 ,409, 375, 233, 12 ,16 ,7, 88, 24 ,15 ,41 ,2)



## additional constraint from meeting 06.11.2023
OBS3_female <- ifelse(OBS3_female < 15, 15, OBS3_female)
OBS3_male <- ifelse(OBS3_male < 15, 15, OBS3_male)


 generation_start <- 21
 generation_end <- 40
 objective_function_year <- 40

 library(yaml)
 library(data.table)
 library(MoBPS)
 library(RandomFieldsUtils)
 library(miraculix)
 library(purrr)
 
 config <- yaml.load_file("config/config.yaml")
 
RandomFieldsUtils::RFoptions(cores = 1)
RandomFieldsUtils::RFoptions(warn_address=FALSE)
print(RFoptions()$basic$cores)

n_simulation = 1
for(setting in 1:n_simulation){

  set.seed(randomSeed)
  
  test_mode <- FALSE #if TRUE size of the big cohorts will be scaled!
  use.recalculate.manual = TRUE
  
  
  if(test_mode){ # test mode was not tested only with large population
    scale <- 0.5 # Scale the size of some cohorts a factor (mostly the large initial DHs)
    chr.nr <- NULL # single chromosome
    log <- TRUE
    nsnp <- 1000
    n.gen <- 1
    options(nwarnings = 10000)
  } else{
    scale <- 1
    chr.nr <- 21
    log <- FALSE
    nsnp <- 10000
    n.gen <- 10  ## Build up LD
  }
  
  bv_male <- list()
  bv_female <- list()
  bv_male_sd <- list()
  bv_female_sd <- list()
  
  cohorts_list_female1 <- list(
    "best_cross_f_",
    "DH0_female_",
    "DH1_female_",
    "OBS1_female_",
    "OBS2_female_"
  )
  cohorts_list_female2 <- list(
    "OBS1_female_rec_",
    "OBS2_female_rec_",
    "OBS3_female_rec_"
  )
  
  if (TC2.2_YT == 1) {
    cohorts_list_female <- c(cohorts_list_female1, "OBS3.1_female_", "OBS3.2_female_", cohorts_list_female2)
  } else {
    cohorts_list_female <- c(cohorts_list_female1, "OBS3.1_female_", cohorts_list_female2)
  }
  
  cohorts_list_male <- list(
    "best_cross_m_",
    "DH0_male_",
    "DH1_male_",
    "OBS1_male_",
    "OBS2_male_",
    "OBS3_male_",
    "OBS1_male_rec_",
    "OBS2_male_rec_",
    "OBS3_male_rec_"
  )
  
  # Founder generation /setup
  ### Creation of traits
  ### trait 1 = Grain Yield, Trait 2 = Protein Content, Trait 3 = Ressesive disease
  burn_in_generation <- 20
  n_traits <- 6 # three trait, 2 different expressions
  n_location <- 10 # number of total locations in breeding program
  n_traits_n_location <- n_traits * n_location # combination of traits and locations
  
  
  for(bv in 1:n_traits_n_location){
    bv_male[[bv]]= matrix(0, nrow=40, ncol=length(cohorts_list_male))
    bv_female[[bv]]= matrix(0, nrow=40, ncol=length(cohorts_list_female))
    bv_male_sd[[bv]]= matrix(0, nrow=40, ncol=length(cohorts_list_male))
    bv_female_sd[[bv]]= matrix(0, nrow=40, ncol=length(cohorts_list_female))
  }
  
  # # Load the burn-in phase to save 20 cycle 
  load("scripts/population_burnin.RData")
  
  {
    # avoid recalculation of the initial new diversity material
    # Later on we want to introduce materials to the breeding program with at least same performance than the initial materials, not worst
    
    pop1 = add.diversity(population,
                         pool.cohorts = paste0("best_cross_m_", -5),
                         export.pop1 = TRUE)
    
    pop2 = add.diversity(population,
                         pool.cohorts = paste0("best_cross_f_", -5),
                         breeding.size = ncross_female,
                         export.pop1 = TRUE)
    
    
    pop3 = add.diversity(population = population,
                         pool.cohorts = c("founder_female", "founder_male"),
                         export.pop1 = TRUE)
  }
  ############################# NEXT ROUND #######################################
  for(i in (generation_start):(generation_end)){
    i_prior1 <- i-1
    i_prior2 <- i-2
    i_prior3 <- i-3
    
    
    ### number of generated DH should be depend on the F1
    F1_male <- ncross_male
    F1_female <- ncross_female
    
    ## same number of DH1 is producing as DH0
    DH1_capacity_male <- DH0_capacity_male
    DH1_capacity_female <- DH0_capacity_female
    
    # Year1
    # Creating F1 progenies
    population <- breeding.diploid(population,
                                   breeding.size=c(F1_male,0),
                                   selection.size=c(ncross_male,0),
                                   selection.m.cohorts = paste0(c("best_cross_m_"), i),
                                   name.cohort = paste0(c("F1_male_dh_"), i),
                                   bv.ignore.traits = 1:(n_traits*n_location),
                                   display.progress = FALSE)
    
    population <- breeding.diploid(population,
                                   breeding.size=c(0,F1_female),
                                   selection.size=c(0,ncross_female),
                                   selection.f.cohorts = paste0(c("best_cross_f_"), i),
                                   name.cohort = paste0(c("F1_female_dh_"), i),
                                   bv.ignore.traits = 1:(n_traits*n_location),
                                   display.progress = FALSE)
    
    if(i >= (-4)){
      # Year2 - Year3
      # DH0 production from F1
      population <- breeding.diploid(population,
                                     breeding.size=round(c(DH1_capacity_male,0) * scale^4),
                                     selection.m.cohorts = paste0(c("F1_male_dh_"), i_prior1),
                                     name.cohort = paste0(c("DH0_male_"), i),
                                     share.genotyped = 1,
                                     dh.mating = TRUE,
                                     same.sex.activ = TRUE,
                                     same.sex.sex = 0,
                                     display.progress = FALSE,
                                     use.recalculate.manual = TRUE)
      
      population <- breeding.diploid(population,
                                     breeding.size=c(0,DH1_capacity_female) * scale^4,
                                     selection.f.cohorts = paste0(c("F1_female_dh_"), i_prior1),
                                     name.cohort = paste0(c("DH0_female_"), i),
                                     share.genotyped = 1,
                                     dh.mating = TRUE,
                                     same.sex.activ = TRUE,
                                     same.sex.sex = 1,
                                     display.progress = FALSE,
                                     use.recalculate.manual = TRUE)
      
    }
    
    if( i >= (-3)){
      
      # Year4
      ### Selection of DH1 is done at random!
      # DH1 production
      #on default all information of individuals is copied. here limiting that to only copy ID (14) and pedigree (46) information of the optional once.
      # gain.stats calculates how much improvement you are making relative to you parents with in 1.11 is automatically done on default.
      population <- breeding.diploid(population,
                                     breeding.size=c(DH1_capacity_male,0)*scale^4,
                                     selection.size=c(DH1_capacity_male,0)*scale^4,
                                     selection.m.cohorts = paste0(c("DH0_male_"), i_prior1),
                                     name.cohort = paste0(c("DH1_male_"), i),
                                     copy.individual.m = TRUE,
                                     copy.individual.use = 14, # copy ID (14) information
                                     copy.individual.use2 = 46, # copypedigree (46) information
                                     gain.stats = FALSE, ## does not work otherwise when parental lines are delete.gen
                                     selection.skip = TRUE,
                                     display.progress = FALSE)
      
      removes_gen <- get.database(population, cohorts = paste0(c("DH0_male_"), i_prior1))[,1]
      population <- breeding.diploid(population, delete.gen = removes_gen)
      
      population <- breeding.diploid(population,
                                     breeding.size=round(c(0,DH1_capacity_female)*scale^4),
                                     selection.size=round(c(0,DH1_capacity_female)*scale^4),
                                     selection.f.cohorts = paste0(c("DH0_female_"), i_prior1),
                                     name.cohort = paste0(c("DH1_female_"), i),
                                     copy.individual.f = TRUE,
                                     gain.stats = FALSE, ## does not work otherwise when parental lines are delete.gen
                                     display.progress = FALSE,
                                     copy.individual.use = 14,
                                     selection.skip = TRUE,
                                     copy.individual.use2 = 46)
      
      ## some clean up as we don't use these informations
      removes_gen <- get.database(population, cohorts = paste0(c("DH0_female_"), i_prior1))[,1]
      population <- breeding.diploid(population, delete.gen = removes_gen)
      
      # Simulate Phenotype for DH1 with heritability 0.4 for Trait 3 per se
      trait.weights <- c(0,0,1,0,0,0)
      index_perse_loc1 <- get.index(population, locations = c(1), location.weights = c(1,rep(0,9)), trait.weights=trait.weights)
      
      ## skip the pred GCA heritability of GY and Pc traits, check if MoBPS calculate the heritability in line with what you want
      ## prediction accuracy is heritability squared
      ## Correlation between genetic values and BVE:
      # 0.69548 0.7145493
      # h2 = 0.5 for pred GCA it is totally fine
      population <- breeding.diploid(population,
                                     heritability = 0.4,
                                     phenotyping.cohorts = paste0(c("DH1_male_","DH1_female_"), i),
                                     variance.correction = "generation.mean",
                                     n.observation = index_perse_loc1)
      
      # DH1 MAS
      # Doing MAS outside of MoBPS manually to avoid unnecessary actions / population list is overwritten multipe times etc.
      mas_geno1 = get.geno(population, cohorts = paste0("DH1_male_", i))
      mas_geno2 = get.geno(population, cohorts = paste0("DH1_female_", i))
      
      y_hat1 = matrix(0, nrow = n_traits*n_location, ncol = ncol(mas_geno1))
      y_hat2 = matrix(0, nrow = n_traits*n_location, ncol = ncol(mas_geno2))
      
      
      for(bven in which(index_perse_loc1 != 0)){
        n_qtl = nrow(population$info$real.bv.add[[bven]])
        mas_marker = sample(1:n_qtl,ceiling(n_qtl*0.15))
        mas.effects = population$info$real.bv.add[[bven]][mas_marker,5]
        
        active.snp = population$info$real.bv.add[[bven]][mas_marker,]
        mas_pos = c(0, population$info$cumsnp)[active.snp[,2]]+active.snp[,1]
        
        y_hat1[bven,] <- mas.effects %*% mas_geno1[mas_pos,]
        y_hat2[bven,] <- mas.effects %*% mas_geno2[mas_pos,]
        
        
      }
      
      population = insert.bve(population, cbind(colnames(mas_geno1), t(y_hat1)))
      population = insert.bve(population, cbind(colnames(mas_geno2), t(y_hat2)))
      
      #only in per se traits
      mas_part_m = colSums(get.bve(population, cohorts = paste0("DH1_male_", i))* index_perse_loc1)
      mas_part_f = colSums(get.bve(population, cohorts = paste0("DH1_female_", i))* index_perse_loc1)
      
      # Agronomic part
      agro_part_m = rnorm(length(mas_part_m))
      agro_part_f = rnorm(length(mas_part_f))
      
      # For all traits
      combined_index_m = mas_part_m / sd(mas_part_m) + agro_part_m / sd(agro_part_m)
      population <- insert.bve(population, bves = cbind(names(combined_index_m), matrix(rep(combined_index_m, times = n_traits_n_location), ncol = n_traits_n_location)))
      
      combined_index_f = mas_part_f / sd(mas_part_f) + agro_part_f / sd(agro_part_f)
      population <- insert.bve(population, bves = cbind(names(combined_index_f), matrix(rep(combined_index_f, times = n_traits_n_location), ncol = n_traits_n_location)))
      
    }
    
    if( i >= (-2)){
      
      # Year5
      # Creating OBS1
      population <- breeding.diploid(population,
                                     breeding.size=c(OBS1_male,0),
                                     selection.size=c(OBS1_male,0),
                                     selection.m.cohorts = paste0("DH1_male_", i_prior1),
                                     name.cohort = paste0(c("OBS1_male_"), i),
                                     copy.individual.m = TRUE,
                                     selection.criteria = "bve",
                                     max.offspring = 1,
                                     multiple.bve.weights.m = index_perse_loc1,
                                     multiple.bve.scale.m = "bve_sd",
                                     gain.stats = FALSE, ## does not work otherwise when parental lines are delete.gen
                                     display.progress = FALSE)
      
      removes_gen <- get.database(population, cohorts = paste0("DH1_male_", i_prior1))[,1]
      population <- breeding.diploid(population, delete.gen = removes_gen)
      
      population <- breeding.diploid(population,
                                     breeding.size=c(0,OBS1_female),
                                     selection.size=c(0,OBS1_female),
                                     selection.f.cohorts = paste0("DH1_female_", i_prior1),
                                     name.cohort = paste0("OBS1_female_", i),
                                     copy.individual.f = TRUE,
                                     selection.criteria = "bve",
                                     max.offspring = 1,
                                     multiple.bve.weights.f = index_perse_loc1,
                                     multiple.bve.scale.f =  "bve_sd",
                                     gain.stats = FALSE, ## does not work otherwise when parental lines are delete.gen
                                     display.progress = FALSE)
      
      removes_gen <- get.database(population, cohorts = paste0("DH1_female_", i_prior1))[,1]
      population <- breeding.diploid(population, delete.gen = removes_gen)
      
      # Genomic values of OBS1s are not automatically calculated because they are copies of DH1s
      # Calculate them here - otherwise genomic variance will be 0
      population <- recalculate.manual(population, cohorts= paste0(c("OBS1_male_","OBS1_female_"), i))
      
      # per se yield trial OBS1
      # variance.correction = different cohorts having different averaging performance
      # to make sure they are all on same level
      # when calculating variance, the values are high due to that,
      # If we estimate genetic variance without correction variance is higher
      
      population <- breeding.diploid(population,
                                     new.bv.observation.cohorts = paste0(c("OBS1_male_","OBS1_female_"), i),
                                     heritability = 0.4,
                                     variance.correction = "generation.mean",
                                     n.observation = index_perse_loc1)
      
      #Evaluation/selection parents (line level), GS
      trait.weights <- c(0,0,1,0,0,0)
      index_perse_loc2 <- get.index(population, locations = c(1,2), location.weights = c(1,1,rep(0,8)), trait.weights=trait.weights)
      all_traits <- 1:n_traits_n_location
      index_of_ones <- which(index_perse_loc2 != 0)
      
      perse_ignore_loc2 <- setdiff(all_traits, index_of_ones)
      
      population <- breeding.diploid(population,
                                     bve=TRUE,
                                     bve.cohorts = paste0(c("OBS1_male_", "OBS1_female_"), i),
                                     bve.insert.cohorts = paste0(c("OBS1_male_", "OBS1_female_"), i),
                                     bve.ignore.traits = perse_ignore_loc2)
      
      bve_parent_m = colSums(get.bve(population, cohorts =  paste0(c("OBS1_male_"), i)) * index_perse_loc2)
      bve_parent_f = colSums(get.bve(population, cohorts =  paste0(c("OBS1_female_"), i))* index_perse_loc2)
      
      #currently for all traits
      combined_index_m = bve_parent_m / sd(bve_parent_m)
      population <- insert.bve(population, bves = cbind(names(combined_index_m), matrix(rep(combined_index_m, times = n_traits_n_location), ncol = n_traits_n_location)))
      
      combined_index_f = bve_parent_f / sd(bve_parent_f)
      population <- insert.bve(population, bves = cbind(names(combined_index_f), matrix(rep(combined_index_f, times = n_traits_n_location), ncol = n_traits_n_location)))
      
      ### TC1 is more like symbolic act to create seeds but has no influence in the design and selection here
      ### the 300 lines that are tested next year for TC1 YT is actually OBS2 lines
      
    }
    
    #### Creation of OBS2 from OBS1 - selection is based on bve GCA
    if( i >= (-1)){
      population <- breeding.diploid(population,
                                     breeding.size=c(OBS2_male,0),
                                     selection.size=c(OBS2_male,0),
                                     selection.m.cohorts = paste0(c("OBS1_male_"), i_prior1),
                                     name.cohort = paste0(c("OBS2_male_"), i),
                                     copy.individual.m = TRUE,
                                     selection.criteria = "bve",
                                     multiple.bve.weights.m = index_perse_loc2,
                                     phenotyping.child ="zero",
                                     display.progress = FALSE)
      
      population <- breeding.diploid(population,
                                     breeding.size=c(0,OBS2_female),
                                     selection.size=c(0,OBS2_female),
                                     selection.f.cohorts = paste0(c("OBS1_female_"), i_prior1),
                                     name.cohort = paste0(c("OBS2_female_"), i),
                                     copy.individual.f = TRUE,
                                     selection.criteria = "bve",
                                     multiple.bve.weights.f = index_perse_loc2,
                                     phenotyping.child ="zero",
                                     display.progress = FALSE)
      
      population <- breeding.diploid(population,
                                     new.bv.observation.cohorts = paste0(c("OBS2_male_","OBS2_female_"), i),
                                     heritability = 0.5,
                                     variance.correction = "generation.mean",
                                     n.observation = index_perse_loc2)
      
      #Evaluation/selection parents (line level) - GCA of the parents
      population <- breeding.diploid(population,
                                     bve=TRUE,
                                     bve.cohorts = paste0(c("OBS2_male_", "OBS2_female_"), i),
                                     bve.insert.cohorts = paste0(c("OBS2_male_", "OBS2_female_"), i),
                                     bve.ignore.traits = perse_ignore_loc2)
      
      
      population <- breeding.diploid(population,
                                     breeding.size= c(2*OBS2_male,0),
                                     selection.size = c(OBS2_male,2),
                                     breeding.all.combination = T,
                                     selection.m.cohorts = paste0(c("OBS2_male_"), i),
                                     selection.f.cohorts = paste0(c("Tester_female_"), max(i,1)),
                                     selection.criteria = "bve",
                                     phenotyping.child ="zero",
                                     name.cohort = paste0(c("TC1_YT_male_"), i),
                                     display.progress = FALSE,
                                     use.recalculate.manual = use.recalculate.manual)
      
      # TC1.1_YT
      population <- breeding.diploid(population,
                                     breeding.size = c(0,OBS2_female*2),
                                     selection.size = c(2,OBS2_female),
                                     breeding.all.combination = T,
                                     selection.f.cohorts = paste0(c("OBS2_female_"), i),
                                     selection.m.cohorts = paste0(c("Tester_male_"), max(1,i)),
                                     selection.criteria = "bve",
                                     phenotyping.child ="zero",
                                     name.cohort = paste0(c("TC1_YT_female_"), i),
                                     display.progress = FALSE,
                                     use.recalculate.manual = use.recalculate.manual)
      
      ## in TC YT trial both 2 traits are present
      ## #Grain yiled has a weight of 3
      trait.weights_tc <- c(0,0,0,3,1,1)
      index_tc_loc2 <- get.index(population, locations = c(1,2), location.weights = c(1,1,rep(0,8)), trait.weights=trait.weights_tc)
      index_of_ones_tc <- which(index_tc_loc2 !=0)
      tc_ignore_loc2 <- setdiff(all_traits, index_of_ones_tc)
      
      population <- breeding.diploid(population,
                                     phenotyping.cohorts = paste0(c("TC1_YT_male_","TC1_YT_female_"), i),
                                     heritability = 0.6,
                                     variance.correction = "generation.mean",
                                     n.observation = index_tc_loc2)
      
      population = breeding.diploid(population,
                                    offspring.bve.parents.cohorts = paste0("OBS2_male_", i),
                                    offspring.bve.offspring.cohorts = paste0("TC1_YT_male_", i)) # import of offspring phenotypes
      
      
      suppressWarnings({
        population = breeding.diploid(population,
                                      bve=TRUE,
                                      bve.cohorts = paste0("OBS2_male_", i),
                                      input.phenotype = "off",
                                      bve.ignore.traits = tc_ignore_loc2)
      })
      
      population = breeding.diploid(population,
                                    offspring.bve.parents.cohorts = paste0("OBS2_female_", i),
                                    offspring.bve.offspring.cohorts = paste0("TC1_YT_female_", i))
      
      suppressWarnings({
        population = breeding.diploid(population,
                                      bve=TRUE,
                                      bve.cohorts = paste0("OBS2_female_", i),
                                      input.phenotype = "off",
                                      bve.ignore.traits = tc_ignore_loc2)
      })
    }
    
    
    if( i >= 0){
      #############################################################################################
      # Select the best lines for OBS3
      if(OBS3_male == 0) {
        OBS3_male <-1
      }
      
      best_lines_m <- breeding.diploid(population,
                                       selection.size = c(OBS3_male,0),
                                       selection.m.cohorts = paste0(c("OBS2_male_"), i_prior1),
                                       selection.criteria = "bve",
                                       multiple.bve.weights.m = index_tc_loc2,
                                       export.selected	= TRUE)[[1]]
      
      fixed_breeding_m = cbind(best_lines_m[,1:3,drop=FALSE],
                               best_lines_m[,1:3,drop=FALSE])
      
      best_lines_f <- breeding.diploid(population,
                                       selection.size = c(0,OBS3_female),
                                       selection.f.cohorts = paste0(c("OBS2_female_"), i_prior1),
                                       selection.criteria = "bve",
                                       multiple.bve.weights.f = index_tc_loc2,
                                       export.selected	= TRUE)[[2]]
      
      
      fixed_breeding_f = cbind(best_lines_f[,1:3,drop=FALSE],
                               best_lines_f[,1:3,drop=FALSE])
      
      # Creation of OBS3
      population <- breeding.diploid(population,
                                     breeding.size=c(OBS3_male,0),
                                     fixed.breeding = fixed_breeding_m,
                                     copy.individual.m = TRUE,
                                     name.cohort = paste0(c("OBS3_male_"), i),
                                     phenotyping.child = "zero",
                                     display.progress = FALSE)
      
      population <- breeding.diploid(population,
                                     breeding.size=c(0,OBS3_female),
                                     fixed.breeding = fixed_breeding_f,
                                     copy.individual.f = TRUE,
                                     name.cohort = paste0(c("OBS3.1_female_"), i),
                                     phenotyping.child = "zero",
                                     display.progress = FALSE)
      
      ## traits Protein content and disease
      trait.weights <- c(0,1,1,0,0,0)
      index_perse_loc3 <- get.index(population, locations = c(1,2,3), location.weights = c(1,1,1), trait.weights=trait.weights)
      
      population <- breeding.diploid(population,
                                     new.bv.observation.cohorts = paste0(c("OBS3_male_","OBS3.1_female_"), i),
                                     heritability = 0.6,
                                     variance.correction = "generation.mean",
                                     n.observation = index_perse_loc3)
      
      population <- breeding.diploid(population,
                                     breeding.size=c(OBS3_male*3,0),
                                     selection.size=c(OBS3_male,3),
                                     breeding.all.combination = T,
                                     selection.m.cohorts = paste0(c("OBS3_male_"), i),
                                     selection.f.cohorts = paste0(c("OBS3.1_female_"), max(0,i_prior1)),
                                     phenotyping.child ="zero",
                                     name.cohort = paste0(c("TC2_YT_male"), i),
                                     display.progress = FALSE,
                                     use.recalculate.manual = use.recalculate.manual)
      
      population <- breeding.diploid(population,
                                     breeding.size=c(0,OBS3_female*3),
                                     selection.size=c(3,OBS3_female),
                                     breeding.all.combination = T,
                                     selection.f.cohorts = paste0(c("OBS3.1_female_"), i),
                                     selection.m.cohorts = paste0(c("OBS3_male_"), i),
                                     phenotyping.child ="zero",
                                     name.cohort = paste0(c("TC2.1_YT_female_"), i),
                                     display.progress = FALSE,
                                     use.recalculate.manual = use.recalculate.manual)
      
      trait.weights_tc <- c(0,0,0,3,1,1)
      index_tc_loc3 <- get.index(population, locations = c(1,2,3), location.weights = c(1,1,1, rep(0,7)), trait.weights=trait.weights_tc)
      index_of_ones_tc <- which(index_tc_loc3 !=0)
      tc_ignore_loc3 <- setdiff(all_traits, index_of_ones_tc)
      
      population <- breeding.diploid(population,
                                     phenotyping.cohorts = paste0(c("TC2_YT_male","TC2.1_YT_female_"), i),
                                     heritability = 0.6,
                                     variance.correction = "generation.mean",
                                     n.observation = index_tc_loc3)
      
      population = breeding.diploid(population,
                                    offspring.bve.parents.cohorts = paste0("OBS3_male_", i),
                                    offspring.bve.offspring.cohorts = paste0(c("TC2_YT_male"), i))
      
      
      
      ## cross traits in location 1 and 2
      suppressWarnings({
        population = breeding.diploid(population,
                                      bve=TRUE,
                                      bve.cohorts =  c(paste0(c("OBS3_male_"), i),
                                                       paste0(c("OBS2_male_"), i_prior1)),
                                      input.phenotype = "off",
                                      bve.ignore.traits = tc_ignore_loc3)
      })
      
      population = breeding.diploid(population,
                                    offspring.bve.parents.cohorts = paste0("OBS3.1_female_", i),
                                    offspring.bve.offspring.cohorts = paste0(c("TC2.1_YT_female_"), i))
      
      suppressWarnings({
        population = breeding.diploid(population,
                                      bve=TRUE,
                                      bve.cohorts =  c(paste0(c("OBS3.1_female_"), i),
                                                       if(exist.cohort(population, paste0("OBS2.2_female_",i_prior1))){paste0("OBS2.2_female_", i_prior1)}
                                                       else{NULL},
                                                       paste0(c("OBS2_female_"), i_prior2)),
                                      input.phenotype = "off",
                                      bve.ignore.traits = tc_ignore_loc3)
        
      })
    }
    
    
    if(i >= (1)){
      
      if(TC2.2_YT == 1) {
        # share of % of what we use in OBS3.1_female
        share_OBS3_1 <- 70
        
        OBS3.2 = round((get.database(population, cohorts = paste0("OBS3.1_female_", i_prior1))[,4]) * share_OBS3_1 /100)
        
        if(OBS3_female == 0) {
          OBS3_female <-1
        }
        
        best_lines_f <- breeding.diploid(population,
                                         selection.size = c(0,OBS3.2),
                                         selection.f.cohorts = paste0(c("OBS3.1_female_"), i_prior1),
                                         selection.criteria = "bve",
                                         multiple.bve.weights.f = index_tc_loc3,
                                         export.selected	= TRUE)[[2]]
        
        fixed_breeding_f = cbind(best_lines_f[,1:3,drop=FALSE],
                                 best_lines_f[,1:3,drop=FALSE])
        
        population <- breeding.diploid(population,
                                       breeding.size=c(0,OBS3.2),
                                       fixed.breeding = fixed_breeding_f,
                                       copy.individual.f = TRUE,
                                       name.cohort = paste0(c("OBS3.2_female_"), i),
                                       phenotyping.child = "zero",
                                       display.progress = FALSE)
        
        # per se yield trial P2
        population <- breeding.diploid(population,
                                       new.bv.observation.cohorts = paste0(c("OBS3.2_female_"), i),
                                       heritability = 0.5,
                                       variance.correction = "generation.mean",
                                       n.observation = index_perse_loc3)
        
        population <- breeding.diploid(population,
                                       breeding.size=c(0,nrow(fixed_breeding_f)*3),
                                       selection.size=c(3,nrow(fixed_breeding_f)),
                                       breeding.all.combination = TRUE,
                                       selection.f.cohorts = paste0(c("OBS3.2_female_"), i),
                                       selection.m.cohorts = paste0(c("OBS3_male_"), i),
                                       phenotyping.child ="zero",
                                       name.cohort = paste0(c("TC2.2_YT_female_"), i),
                                       display.progress = FALSE)
        
        population <- breeding.diploid(population,
                                       phenotyping.cohorts = paste0(c("TC2.2_YT_female_"), i),
                                       heritability = 0.6,
                                       variance.correction = "generation.mean",
                                       n.observation = index_tc_loc3)
        
        population = breeding.diploid(population,
                                      offspring.bve.parents.cohorts = paste0("OBS3.2_female_", i),
                                      offspring.bve.offspring.cohorts = paste0(c("TC2.2_YT_female_"), i))
        
        suppressWarnings({
          population = breeding.diploid(population,
                                        bve=TRUE,
                                        bve.cohorts = c(paste0(c("OBS3.2_female_"), i),
                                                        paste0(c("OBS2.2_female_"), i_prior1),
                                                        if(exist.cohort(population, paste0("OBS2.2_female_",i_prior2))){paste0(c("OBS2.2_female_"), i_prior2)}
                                                        else{NULL},
                                                        paste0(c("OBS2_female_"), i_prior3)),
                                        input.phenotype = "off",
                                        bve.ignore.traits = tc_ignore_loc3)
        })
        
      }
      
      ### It is not required to go further back as this is just product development and not pool improvement
      i_prior1 <- max(i-1,1)
      i_prior2 <- max(i-2,1)
      i_prior3 <- max(i-3,1)
      #############################################################################################
      
    }
    
    if( i <= 0){
      # this number should be calculated based on the breeding scheme without and then check how they change
      # The motivation behind was that Generation-6 should be worse than Generation 0. The number put in there was the gain we observed in one breeding cycle.
      # Should be totally fine to use same setting even for different weights. This would then mimic a change in Index weights at time 0. Its mostly about have some gain to avoid that Generation 1 to 6 basically the same genetic values .
      target.value = rep(c(1.38,1.9,3.16,1.38,2.20,4.48), each = n_location) * (i + 6)/7 + rowMeans(get.bv(population, cohorts = paste0("best_cross_m_", -5)))
      population = add.diversity(population, target.value = target.value,
                                 pool.cohorts = paste0("best_cross_m_", -5),
                                 breeding.size = ncross_male,
                                 reduction.multiplier = 1, selection.rate = if(i<=(-5)){0.7} else {0.5},
                                 sex.quota = 0,
                                 name.cohort = paste0("best_cross_m_", i + 1),
                                 use.recalculate.manual = use.recalculate.manual,
                                 pop1 = pop1)
      
      target.value = rep(c(1.38,1.9,3.16,1.38,2.20,4.48), each = n_location) * (i + 6)/7 + rowMeans(get.bv(population, cohorts = paste0("best_cross_f_", -5)))
      population = add.diversity(population, target.value = target.value,
                                 pool.cohorts = paste0("best_cross_f_", -5),
                                 breeding.size = ncross_female,
                                 reduction.multiplier = 1, selection.rate = if(i<=(-5)){0.7} else {0.5},
                                 sex.quota = 1,
                                 name.cohort = paste0("best_cross_f_", i + 1),
                                 use.recalculate.manual = use.recalculate.manual,
                                 pop1 = pop2)
      
    } else{
      #######################################################################################
      ###########################   lines for the next cycle  ###############################
      #########################   competitor, prebreeding 20%  ##############################
      #######################################################################################
      Recycle_M <- ceiling(ncross_male * 0.2)*2
      Recycle_F <- ceiling(ncross_female * 0.2)*2
      
      population = add.diversity(population = population,
                                 pool.cohorts = c("founder_male", "founder_female"),
                                 target.cohorts = paste0(c("best_cross_m_","best_cross_f_"), i),
                                 reduction.multiplier = 1,
                                 selection.rate = 0.5,
                                 breeding.size = c(Recycle_M, Recycle_F),
                                 name.cohort = paste0(c("New_material_M_", "New_material_F_"), i),
                                 use.recalculate.manual = use.recalculate.manual,
                                 pop1 = pop3)
      
      #############################################################################
      number_line_f <- c(female_rec_OBS3,female_rec_OBS2,female_rec_OBS1)
      how_many_percent_f <- ncross_female - (number_line_f[1] + number_line_f[2] + number_line_f[3] + (Recycle_F/2))
      zzz <- which.max(number_line_f)
      number_line_f[zzz] <- number_line_f[zzz] + how_many_percent_f
      
      female_OBS3 <- number_line_f[1]*2
      female_OBS2 <- number_line_f[2]*2
      female_OBS1 <- number_line_f[3]*2
      
      
      number_line_m <- c(male_rec_OBS3,male_rec_OBS2,male_rec_OBS1)
      how_many_percent_m <- ncross_male - (number_line_m[1:1] + number_line_m[2:2] + number_line_m[3:3] + (Recycle_M/2))
      zz <- which.max(number_line_m)
      number_line_m[zz] <- number_line_m[zz] + how_many_percent_m
      
      male_OBS3 <- number_line_m[1]*2
      male_OBS2 <- number_line_m[2]*2
      male_OBS1 <- number_line_m[3]*2
      ###########################################################################>
      population <- breeding.diploid(population,
                                     selection.m.cohorts = paste0(c("OBS3_male_"), i),
                                     copy.individual.m = T,
                                     selection.criteria = "bve",
                                     name.cohort = paste0("OBS3_male_rec_", i),
                                     display.progress = FALSE)
      
      new1m <- get.pedigree(population,
                            cohorts = paste0(c("OBS3_male_rec_"), i),
                            raw=TRUE)[,1:3]
      
      export2 <- breeding.diploid(population,
                                  selection.size = c(min(male_OBS2,50),0),
                                  selection.m.cohorts = paste0(c("OBS2_male_"), i),
                                  selection.criteria = "bve",
                                  display.progress = FALSE,
                                  export.selected = TRUE)
      
      population <- breeding.diploid(population,
                                     fixed.breeding = cbind(export2[[1]][,1:3],
                                                            export2[[1]][,1:3],0),
                                     copy.individual = TRUE,
                                     name.cohort = paste0("OBS2_male_rec_", i),
                                     display.progress = FALSE)
      
      new2m <- get.pedigree(population,
                            cohorts = paste0(c("OBS2_male_rec_"), i),
                            raw=TRUE)[,1:3]
      
      export3 <- breeding.diploid(population,
                                  selection.size = c(min(male_OBS1,50),0),
                                  selection.m.cohorts = paste0(c("OBS1_male_"), i),
                                  selection.criteria = "bve",
                                  display.progress = FALSE,
                                  export.selected = TRUE)
      
      population <- breeding.diploid(population,
                                     fixed.breeding = cbind(export3[[1]][,1:3],
                                                            export3[[1]][,1:3],0),
                                     copy.individual = TRUE,
                                     name.cohort = paste0("OBS1_male_rec_", i),
                                     display.progress = FALSE)
      
      new3m <- get.pedigree(population,
                            cohorts = paste0(c("OBS1_male_rec_"), i),
                            raw=TRUE)[,1:3]
      
      new4m <- breeding.diploid(population,
                                selection.size = c((Recycle_M),0), # 20%
                                selection.m.cohorts = paste0(c("New_material_M_"), i),
                                export.selected = TRUE,
                                display.progress = FALSE)[[1]][,1:3]
      
      
      sel_p2 <- sample(nrow(new1m), male_OBS3, replace=TRUE)
      p2_lines_cross <- new1m[sel_p2,]
      
      ebvs = export2[[1]][,4]
      ebvs = (ebvs - min(ebvs)) / (max(ebvs) - min(ebvs)) * 4 + 1 # Line with the best EBV will be selected 5 times as frequently
      
      sel_p1 <- sample(nrow(new2m), male_OBS2, replace=TRUE, prob = ebvs)
      p1_lines_cross <- new2m[sel_p1,]
      
      ebvs = export3[[1]][,4]
      ebvs = (ebvs - min(ebvs)) / (max(ebvs) - min(ebvs)) * 4 + 1 # Line with the best EBV will be selected 5 times as frequently
      
      sel_p0 <- sample(nrow(new3m), male_OBS1, replace=TRUE, prob = ebvs)
      p0_lines_cross <- new3m[sel_p0,]
      
      all_rec_m <- rbind(p2_lines_cross,p1_lines_cross,p0_lines_cross)
      
      # shuffle the rows
      perm <- sample(nrow(all_rec_m))
      all_rec_m_shuffled <- all_rec_m[perm, ]
      
      all_lines_m <- rbind(new4m,all_rec_m_shuffled)
      
      fixed_breeding_cross_m <- cbind(all_lines_m[1:(nrow(all_lines_m)/2),],
                                      all_lines_m[((nrow(all_lines_m)/2)+1):nrow(all_lines_m),])
      
      population <- breeding.diploid(population,
                                     breeding.size = c(ncross_male,0),
                                     fixed.breeding = fixed_breeding_cross_m,
                                     name.cohort = paste0("best_cross_m_", i + 1),
                                     display.progress = FALSE,
                                     use.recalculate.manual = use.recalculate.manual)
      
      
      population <- breeding.diploid(population,
                                     selection.f.cohorts = paste0(c(if(exist.cohort(population, paste0("OBS3.2_female_",i))){paste0(c("OBS3.2_female_"))}else{paste0("OBS3.1_female_")}), i),
                                     copy.individual.f = T,
                                     selection.criteria = "bve",
                                     name.cohort = paste0("OBS3_female_rec_", i),
                                     display.progress = FALSE)
      
      
      new1f <- get.pedigree(population, cohorts = paste0(c("OBS3_female_rec_"), i), raw=TRUE)[,1:3]
      
      
      export2 <- breeding.diploid(population,
                                  selection.size = c(0,min(female_OBS2, 50)),
                                  selection.f.cohorts = paste0(c(if(exist.cohort(population, paste0("OBS2.2_female_",i))){paste0(c("OBS2.2_female_"))}else{paste0("OBS2_female_")}), i),
                                  selection.criteria = "bve",
                                  display.progress = FALSE,
                                  export.selected = TRUE)
      
      population <- breeding.diploid(population,
                                     fixed.breeding = cbind(export2[[2]][,1:3],
                                                            export2[[2]][,1:3],1),
                                     copy.individual = TRUE,
                                     name.cohort = paste0("OBS2_female_rec_", i),
                                     display.progress = FALSE)
      
      
      new2f <- get.pedigree(population,
                            cohorts = paste0(c("OBS2_female_rec_"), i),
                            raw=TRUE)[,1:3]
      
      export3 <- breeding.diploid(population,
                                  selection.size = c(0,min(female_OBS1, 50)),
                                  selection.f.cohorts = paste0(c("OBS1_female_"), i),
                                  selection.criteria = "bve",
                                  display.progress = FALSE,
                                  export.selected = TRUE)
      
      population <- breeding.diploid(population,
                                     fixed.breeding = cbind(export3[[2]][,1:3],
                                                            export3[[2]][,1:3],1),
                                     copy.individual = TRUE,
                                     name.cohort = paste0("OBS1_female_rec_", i),
                                     display.progress = FALSE)
      
      
      
      new3f <- get.pedigree(population,
                            cohorts = paste0(c("OBS1_female_rec_"), i),
                            raw=TRUE)[,1:3]
      
      new4f <- breeding.diploid(population,
                                selection.size = c(0,(Recycle_F)), # 20%
                                selection.f.cohorts = paste0(c("New_material_F_"), i),
                                export.selected = TRUE,
                                display.progress = FALSE)[[2]][,1:3]
      
      
      sel_p2 <- sample(nrow(new1f), female_OBS3, replace=TRUE)
      p2_lines_cross <- new1f[sel_p2,]
      
      ebvs = export2[[2]][,4]
      ebvs = (ebvs - min(ebvs)) / (max(ebvs) - min(ebvs)) * 4 + 1 # Line with the best EBV will be selected 5 times as frequently
      
      sel_p1 <- sample(nrow(new2f), female_OBS2, replace=TRUE, prob = ebvs)
      p1_lines_cross <- new2f[sel_p1,]
      
      ebvs = export3[[2]][,4]
      ebvs = (ebvs - min(ebvs)) / (max(ebvs) - min(ebvs)) * 4 + 1 # Line with the best EBV will be selected 5 times as frequently
      
      sel_p0 <- sample(nrow(new3f), female_OBS1, replace=TRUE, prob = ebvs)
      p0_lines_cross <- new3f[sel_p0,]
      
      all_rec_f <- rbind(p2_lines_cross,p1_lines_cross,p0_lines_cross)
      
      # shuffle the rows
      perm <- sample(nrow(all_rec_f))
      all_rec_f_shuffled <- all_rec_f[perm, ]
      
      all_lines_f <- rbind(new4f,all_rec_f_shuffled)
      
      fixed_breeding_cross_f <- cbind(all_lines_f[1:(nrow(all_lines_f)/2),],
                                      all_lines_f[((nrow(all_lines_f)/2)+1):nrow(all_lines_f),])

      
      population <- breeding.diploid(population,
                                     breeding.size = c(0,ncross_female),
                                     fixed.breeding = fixed_breeding_cross_f,
                                     name.cohort = paste0("best_cross_f_", i + 1),
                                     display.progress = FALSE,
                                     use.recalculate.manual = use.recalculate.manual)
      #################################################################################
      ####################################################################################
      # Testers for the next round
      # Testers can be changed to another that they be different
      # high chance that at the end this testers are the same!
      population <- breeding.diploid(population,
                                     selection.size = c(2,0),
                                     selection.m.cohorts = paste0(c("OBS3_male_"), i),
                                     copy.individual.m = T,
                                     selection.criteria = "bve",
                                     name.cohort = paste0("Tester_male_", i + 1),
                                     display.progress = FALSE)
      
      population <- breeding.diploid(population,
                                     selection.size = c(0,2),
                                     selection.f.cohorts = if(exist.cohort(population, paste0("OBS3.2_female_",i))){paste0(c("OBS3.2_female_"), i)}
                                     else{paste0(c("OBS3.1_female_"), i)},
                                     copy.individual.f = T,
                                     selection.criteria = "bve",
                                     name.cohort = paste0("Tester_female_", i + 1),
                                     display.progress = FALSE)
    }
    considered_index <- i


    # Calculate Results for only the current generation
    # objects to write results into have to be generated before!
    if(i>=1){
      for (index in considered_index) {
        for (p in 1:length(cohorts_list_female)) {
          cohorts <- paste0(cohorts_list_female[[p]], index)
          for (bv in 1:n_traits_n_location) {
            bvs = get.bv(population, cohorts = cohorts)[bv,]
            bvs[bvs==0] = NA
            suppressWarnings(bv_female[[bv]][index, p] <- mean(bvs, na.rm=TRUE))
            suppressWarnings(bv_female_sd[[bv]][index, p] <- sqrt(var(bvs, na.rm = TRUE)))
          }
        }
      }
    }
    
    if(i>=1){
      for (index in considered_index) {
        for (p in 1:length(cohorts_list_male)) {
          cohorts <- paste0(cohorts_list_male[[p]], index)
          for (bv in 1:n_traits_n_location) {
            bvs = get.bv(population, cohorts = cohorts)[bv,]
            bvs[bvs==0] = NA
            suppressWarnings(bv_male[[bv]][index, p] <- mean(bvs, na.rm=TRUE))
            suppressWarnings(bv_male_sd[[bv]][index, p] <- sqrt(var(bvs, na.rm = TRUE)))
          }
        }
      }
    }
  
  
  ############## Create the objective function of interest ###################
  # at the moment we have an average of all traits
  base_rate <- 0.05
  
  if(i==generation_end){
    #sum of all traits
    genetic_gain_male_all <- reduce(bv_male[1:n_traits_n_location], `+`)/n_traits_n_location
    genetic_gain_female_all <- reduce(bv_female[1:n_traits_n_location], `+`)/n_traits_n_location
    
    # population is in column 1
    genetic_gain_male_all <- genetic_gain_male_all[generation_start:generation_end,1]
    genetic_gain_female_all <- genetic_gain_female_all[generation_start:generation_end,1]
    
    gain_value_male <- genetic_gain_male_all / (1 + base_rate)^(generation_start:generation_end)
    final_value_male <- sum(gain_value_male)
    
    gain_value_female <- genetic_gain_female_all / (1 + base_rate)^(generation_start:generation_end)
    final_value_female <- sum(gain_value_female)
    
    
    cat("Gain population male overall", i, ":", (genetic_gain_male_all), "\n")
    cat("Gain population female overall", i, ":", (genetic_gain_female_all), "\n")
    
    gain_target_population <- (final_value_male+final_value_female)/2
  }
  
  
  
  
  # This can be changed to whatever cohorts we want to check the diversity for it and control the rate of inbreeding  
  if(i==generation_end){
    
    # geno_male <- get.geno(population, cohorts = paste0("best_cross_m_", objective_function_year ))
    # minor_allele_freq_male <- rowMeans(geno_male)/2
    # minor_allele_freq_male[minor_allele_freq_male>0.5] <- 1 - minor_allele_freq_male[minor_allele_freq_male>0.5]
    # 
    # 
    # geno_female <- get.geno(population, cohorts = paste0("best_cross_f_", objective_function_year))
    # minor_allele_freq_female <- rowMeans(geno_female)/2
    # minor_allele_freq_female[minor_allele_freq_female>0.5] <- 1 - minor_allele_freq_female[minor_allele_freq_female>0.5]
    # 
    # 
    kin_stock_f <- kinship.emp.fast(population=population, cohorts=paste0("best_cross_f_", i))
    kin_stock_m <- kinship.emp.fast(population=population, cohorts=paste0("best_cross_m_", i))
    kin_stock_f[is.na(kin_stock_f)] <- 0
    kin_stock_m[is.na(kin_stock_m)] <- 0
    
    ## average inbreeding coefficient Fx
    kin_stock_f_coe <- round(2*kin_stock_f[2] - 1, 3) * sum(1/(1+base_rate)^generation_start:generation_end)
    kin_stock_m_coe <- round(2*kin_stock_m[2] - 1, 3) * sum(1/(1+base_rate)^generation_start:generation_end)
    
    cat("Inbreeding population male year", i, ":", (kin_stock_m_coe), "\n")
    cat("Inbreeding population female year", i, ":", (kin_stock_f_coe), "\n")
    
    
    
    inbreeding_target <- (kin_stock_m_coe + kin_stock_f_coe) /2
    target_function <- gain_target_population  
  }
  
  }
  
  
  results <- t(c(rep, 
                 randomSeed, 
                 TC2.2_YT,
                 ncross_female,
                 ncross_male,
                 DH0_capacity_female,
                 DH0_capacity_male,
                 OBS1_female,
                 OBS1_male,
                 OBS2_female,
                 OBS2_male,
                 OBS3_female,
                 OBS3_male,
                 female_rec_OBS3,
                 female_rec_OBS2,
                 female_rec_OBS1,
                 male_rec_OBS3,
                 male_rec_OBS2,
                 male_rec_OBS1,
                 inbreeding_target, 
                 target_function))
  
  colnames(results) <- c("rep",
                         "randomSeed",
                         "TC2.2_YT",
                         "ncross_female",
                         "ncross_male",
                         "DH0_capacity_female",
                         "DH0_capacity_male",
                         "OBS1_female",
                         "OBS1_male",
                         "OBS2_female",
                         "OBS2_male",
                         "OBS3_female",
                         "OBS3_male",
                         "female_rec_OBS3",
                         "female_rec_OBS2", 
                         "female_rec_OBS1",
                         "male_rec_OBS3",
                         "male_rec_OBS2", 
                         "male_rec_OBS1",
                         "inbreeding_target",
                         "target_function")
  print(results)
  save(list=c("results"), file = output)

# Save the specified objects to a file named based on the variables `rep`, `randomSeed`, `TC1.2_YT`, and `TC2.2_YT`
#save(file = paste0("Pop_bve_", rep, "_", randomSeed, "_", TC1.2_YT, "_", TC2.2_YT, "_bvs.RData"), 
#     list = c("bv_male", "bv_female", "bv_male_sd", "bv_female_sd"))
}


