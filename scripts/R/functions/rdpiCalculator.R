#' Calculate "relative distance plasticity index" 
#' 
#' This index is taken from from Valladares F .. Zavala MA (2006) Journal of Ecology 94:1103
#' 
#' @param a trait value in environment 1
#' @param b trait value in environment 2
rdpiCalculator = function(a, b){  
  # First calculate the pairwise differences between the two vectors
  pair_dif = c()
  pair_sum = c()
  
  for(i in 1:length(a)){
    for(j in 1:length(b)){
      pair_dif = c(pair_dif, a[i] - b[j])
      pair_sum = c(pair_sum, a[i] + b[j])
    }
  }
  
  # Then sum and divide by number of pairwise comparisons
  rdpi = sum(pair_dif/pair_sum, na.rm=T)/length(pair_dif)
  
  return(rdpi)
}
