generate_new <- function(){
  ### this can be as exact as the initial bound in sampleScript or new bounds
  ### As initial bounds exclude the area of optimum we add them here to check a broader range
  ### this will be activated only if n_off_random in iteration.csv in config folder will be > 0 
  x <- numeric(4)
  x[4] <- rbinom(1,1,0.5)
  x[1] <- sample(100:700,1)
  x[2] <- floor((10000000-3000*x[1])/(4000 + 1000*x[4]))
  x[3] <- sample(15:25, 1)

  x <- c(x[2],x[1],x[3],x[4])
  
  return(x)
}
