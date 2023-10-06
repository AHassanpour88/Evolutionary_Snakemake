approx_density <- function(x, bw=1, results_smooth){
  
  distance <- numeric(nrow(results_smooth))
  
  for(index in 1:(length(param_cols))){ 
    distance <- distance + ((results_smooth[,param_cols[index]] - x[param_cols[index]])/bw[index])^2
  }
  
  weight <- dnorm_approx(sqrt(distance))
  
  target <-  sum(weight)
  
  return (target)
}