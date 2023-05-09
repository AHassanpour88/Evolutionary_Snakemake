## n = how many parameters
## the number that sum of all n should be = m

randomN = function(n, m){
  
  ns = round(sort(runif(n-1, min=0.5, max = m+0.5)))
  rand = c(ns[1], diff(ns))
  rand = c(rand, m - sum(rand))
}


