# This is a faster version of dnorm()
dnorm_approx <- function(distance, bw = 2){
  exp(-distance^2/(2*(bw^2))) / bw
}