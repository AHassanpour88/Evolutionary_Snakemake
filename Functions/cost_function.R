## define the cost of your breeding programs based on different parameters
# for WCGALP the cost is as a constraint and depend on the variables
calc_cost <- function(x){
  
  cost <-
   x[1] * 4000 + # DH0 production
   x[2] * 3000
  return(unname((cost)))
}

# # DHO female
# x[3] * 25 + # DH0 production
#   x[3] * 0.4 * 3 + # Headrows in Nursery: loc 1 reps 1
#   x[3] * 0.06 * 5 + # MAS
#   # DH0 male
#   x[4] * 25 + # DH0 production
#   x[4] * 0.65 * 8 + # Headrows in Nursery
#   x[4] * 0.17 * 25 + # MAS
#   # P0 female
#   x[5] * 60 + # P0_female --> 1400 lines
#   x[5] * 1.2 * 20 + #TP0 --> 1680 lines: loc 1 reps 1
#   # P0 male
#   x[6] * 60 + # P0_male + --> 2000 lines
#   x[6] * 1.2 * 20 + #TP0 --> 2400 lines
#   # P1.1 female
#   x[7] * 30  + #Parent Multiplication (PM), microplot: loc 1 reps 1 + MAS
#   x[7] * 1.2 * 25 + # Nursery Disease trials: loc 1 reps 1.2
#   x[7] * 1.2 * 2 * 0.8 * 40 + # Hybrid yield trial: loc 2 rep 0.8
#   x[7] * 1.2 * 25 + # Nursery Disease trials: loc 1 reps 1
#   x[7] * 2.5 * 100 + # TP1: loc 1 reps 1
#   x[7] * 0.57 * 30 + # CMS initial crosses
#   # P1.2 female
#   x[7] * 0.8 * 4 * 25 * x[1] + # Parent Multiplication (PM), microplot: loc 1 reps 4
#   x[7] * 0.8 * 3 * 25 * x[1] + # Nursery Disease trial: loc 1 reps 3
#   x[7] * 2 * 3 * 1.2 * 40 * x[1] + # Hybrid yield trial: loc 3 rep 1.2
#   x[7] * 2 * 125 * x[1] + # Nursery Disease trial: loc 1 reps 3 + TP1
#   x[7] * 0.57 * 4 * 4 * 15 * x[1] + # Sterility test/cross B*A line: loc 4 reps 4
#   # P1 male
#   x[8] * 30  + #Parent Multiplication (PM), microplot: loc 1 reps 1 + MAS
#   x[8] * 1.2 * 25 + # Nursery Disease trials: loc 1 reps 1.2
#   x[8] * 1.2 * 2 * 0.8 * 35 + # Hybrid yield trial: loc 2 rep 0.8
#   x[8] * 1.2 * 25 + # Nursery Disease trials: loc 1 reps 1
#   x[8] * 2.5 * 100 + # TP1: loc 1 reps 1
#   # P2.1 female
#   x[9] * 4 * 25  + # Parent Multiplication (PM), microplot: loc 1 reps 4
#   x[9] * 4 * 2 * 25 + # Nursery Disease trials: loc 4 reps 2
#   x[9] * 3 * 80 + # Parent perse baking quality trial + Parent perse disease trial: loc 3 reps 1
#   x[9] * 15 + # HD-GBS, full trait panel (Parent Characterization)
#   x[9] * 2.5 * 6 * 2 * 40 + # Hybrid yield trial: loc 6 reps 2
#   x[9] * 2.5 * 2 * 40 + # Hybrid yield trial: loc 2 reps 1
#   x[9] * 2.5 * 4 * 2 * 25 + # Nursery Disease trials: loc 4 reps 2
#   x[9] * 7 * 200 + # TP2: loc 1 reps 1
#   x[9] * 200 + # CMS conversion BC1F1 + BC2F1
#   x[9] * 1000 + # 2 times MABC
#   # P2.2 female
#   x[9] * 0.6 * 4 * 25 * x[2] + # Parent Multiplication (PM), microplot: loc 1 reps 4
#   x[9] * 0.6 * 4 * 2 * 25 * x[2] + # Nursery Disease trial: loc 4 reps 2
#   x[9] * 0.6 * 3 * 1 * 80 * x[2] + # Parent perse baking quality trial + Parent perse disease trial: loc 3 rep 1
#   x[9] * 0.6 * 200 * x[2] + # Nursery Disease trial: loc 1 reps 3 + TP1
#   x[9] * 5 * 6 * 2 * 40 * x[2] + # Hybrid yield trial: loc 6 reps 2
#   x[9] * 5 * 2 * 40 * x[2] + # Hybrid yield trial: loc 2 reps 1
#   x[9] * 5 * 4 * 2 * 25 * x[2] + # Nursery Disease trial: loc 4 reps 2
#   x[9] * 5 * 200 * x[2] + # TP2: loc 1 reps 1
#   x[9] * 200 * x[2] + # CMS conversion BC3F1 + BC4F1
#   x[9] * 1000 * x[2] + # 2 times MABC
#   # P2 male
#   x[10] * 4 * 500 + # R G0 maintenance
#   x[10] * 4 * 2 * 25 + # Nursery Disease trials: loc 4 reps 2
#   x[10] * 3 * 70 + # Parent perse baking quality trial + Parent perse disease trial: loc 3 reps 1
#   x[10] * 15 + # HD-GBS, full trait panel (Parent Characterization)
#   x[10] * 2.5 * 6 * 2 * 35 + # Hybrid yield trial CHA/CMS: loc 6 reps 2
#   x[10] * 2.5 * 2 * 35 + # Hybrid yield trial CHA/CMS: loc 2 reps 1
#   x[10] * 2.5 * 4 * 2 * 25 + # Nursery Disease trials: loc 4 reps 2
#   x[10] * 12 * 200 + # TP2 CMS: loc 1 reps 1
#   # number of crosses female
#   x[11] * 60 +
#   # number of crosses male
#   x[12] * 60