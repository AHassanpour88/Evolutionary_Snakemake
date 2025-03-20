rm(list = ls())
# install.packages(pkgs = "AlphaSimR")
library(package = "AlphaSimR")


nCrosses = 100 # Number of crosses per year
nDH      = 100 # DH lines produced per cross
nPYT     = 500 # Entries per preliminary yield trial
nAYT     = 50  # Entries per advanced yield trial
nEYT     = 10  # Entries per elite yield trial (edited) 

# ---- Number of simulation replications ----
nReps = 1

for(REP in 1:nReps){
  cat("Working on REP:", REP,"\n")
  set.seed(REP)
  # Global Parameters
  # ---- Number of simulation replications and breeding cycles ----
  nReps   = 1    # Number of simulation replicates
  nBurnin = 20   # Number of years in burnin phase
  nFuture = 20   # Number of years in future phase
  nCycles = nBurnin + nFuture
  startTP = 19   # Year to start training population
  
  # ---- Genome simulation ----
  nChr = 10
  nQtl = 1000 # Number of QTL per chromosome
  nSnp = 400  # Number of SNP per chromosome
  
  # ---- Initial parents mean and variance ----
  initMeanG  = 1
  initVarG   = 1
  initVarEnv = 1e-6 # Virtually zero for consistency with 2-Part paper
  initVarGE  = 2
  varE       = 4 # Yield trial error variance, bushels per acre
  # Relates to error variance for an entry mean
  
  # ---- Breeding program details ----
  nParents = 50  # Number of parents to start a breeding cycle
  famMax   = 10  # The maximum number of DH lines per cross to enter PYT

  # Effective replication of yield trials
  repHDRW  = 4/9 # h2 = 0.1
  repPYT   = 1   # h2 = 0.2
  repAYT   = 4   # h2 = 0.5
  repEYT   = 8   # h2 = 0.7
  

  # ---- Create initial parents ----
# Create founders

# Generate initial haplotypes
founderPop = runMacs(nChr = nChr,
                     nInd     = nParents,
                     segSites = nQtl + nSnp,
                     inbred   = TRUE,
                     species  = "WHEAT")
SP = SimParam$new(founderPop)

# Add SNP chip
SP$restrSegSites(nQtl, nSnp)
if (nSnp > 0) {
  SP$addSnpChip(nSnp)
}

# Add traits: trait represents yield
SP$addTraitAG(nQtlPerChr = nQtl,
              mean       = initMeanG,
              var        = initVarG,
              varEnv     = initVarEnv,
              varGxE     = initVarGE)

# Collect pedigree
SP$setTrackPed(TRUE)

# Create founder parents
Parents = newPop(founderPop)

# Add phenotype reflecting evaluation in EYT
Parents = setPheno(Parents, varE = varE, reps = repEYT)
rm(founderPop)


  # ---- Fill breeding pipeline with unique individuals from initial parents ----
# Fill breeding pipeline

# Set initial yield trials with unique individuals
for(cohort in 1:7){
  cat("  FillPipeline stage:",cohort,"of 7\n")
  if(cohort<7){
    # Stage 1
    F1 = randCross(Parents, nCrosses)
  }
  if(cohort<6){
    # Stage 2
    DH = makeDH(F1, nDH)
  }
  if(cohort<5){
    # Stage 3
    HDRW = setPheno(DH, varE = varE, reps = repHDRW)
  }
  if(cohort<4){
    # Stage 4
    PYT = selectWithinFam(HDRW, famMax)
    PYT = selectInd(PYT, nPYT)
    PYT = setPheno(PYT, varE = varE, reps = repPYT)
  }
  if(cohort<3){
    # Stage 5
    AYT = selectInd(PYT, nAYT)
    AYT = setPheno(AYT, varE = varE, reps = repAYT)
  }
  if(cohort<2){
    # Stage 6
    EYT = selectInd(AYT, nEYT)
    EYT = setPheno(EYT, varE = varE, reps = repEYT)
  }
  if(cohort<1){
    # Stage 7
  }
}


  # ---- Create a data frame to track key parameters----
  output = data.frame(year     = 1:nCycles,
                      meanG  = numeric(nCycles),
                      varG   = numeric(nCycles),
                      accSel = numeric(nCycles))

  # ---- Burn-in phase: Phenotypic selection program ----
  # NOTE: some scripts below are for phenoypic selection, but saved in genomic
  # selection folder!
  cat("--> Working on Burn-in with Phenotypic selection program \n")
  for(year in 1:nBurnin) {
    cat(" Working on burn-in year:",year,"\n")

# Update parents

# Replace 10 oldest inbred parents with 10 new parents from EYT stage
Parents = c(Parents[11:nParents], EYT)

# Advance year

# Advance breeding program by 1 year
# Works backwards through pipeline to avoid copying data

# Stage 7
# Release variety

# Stage 6
EYT = selectInd(AYT, nEYT)
EYT = setPheno(EYT, varE = varE, reps = repEYT)

# Stage 5
AYT = selectInd(PYT, nAYT)
AYT = setPheno(AYT, varE = varE, reps = repAYT)

# Stage 4
output$accSel[year] = cor(HDRW@gv, HDRW@pheno)
PYT = selectWithinFam(HDRW, famMax)
PYT = selectInd(PYT, nPYT)
PYT = setPheno(PYT, varE = varE, reps = repPYT)

# Stage 3
HDRW = setPheno(DH, varE = varE, reps = repHDRW)

# Stage 2
DH = makeDH(F1, nDH)

# Stage 1
F1 = randCross(Parents, nCrosses)

# Store training population

if (year == startTP){
  cat("  Start collecting training population \n")
  TrainPop = c(PYT, EYT, AYT)
}

if (year > startTP & year < nBurnin+1){
  cat("  Collecting training population \n")
  TrainPop = c(TrainPop,
               PYT, EYT, AYT)
}

if (year > nBurnin){
  cat("  Maintaining training population \n")
  TrainPop = c(TrainPop[-c(1:c(PYT, EYT, AYT)@nInd)],
               PYT, EYT, AYT)
}

    # Report results
    output$meanG[year] = meanG(DH)
    output$varG[year]  = varG(DH)
  }
  

  Scenario  <- output$scenario <- "Phenotypic_Selection"
  file.name <- paste0("Results_",Scenario,".csv")

  
  # Save burn-in to load later for subsequent scenarios
  save.image("population_burnin.RData")

}
