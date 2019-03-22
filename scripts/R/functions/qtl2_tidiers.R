#' QTL2 tidiers
#' 
if(!require(broom))  stop("Install broom package")

tidy.scan1coef <- function(x){
  
  # Convert coefficients to tibble
  coefs <- tibble::as_tibble(x, rownames = "marker") 
  coefs <- tidyr::gather(coefs, "coef", "estimate", -marker)
  
  # Convert SEs to tibble
  # and join together
  if("SE" %in% names(attributes(x))){
    SEs <- tibble::as_tibble(attr(x, "SE"), rownames = "marker")
    SEs <- tidyr::gather(SEs, "coef", "SE", -marker)
    
    coefs <- merge(coefs, SEs, by = c("marker", "coef"))
  }
  
  return(tibble::as_tibble(coefs))
  
}


