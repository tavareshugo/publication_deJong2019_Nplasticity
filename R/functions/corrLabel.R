#' Calculate correlation between two variables and return a pretty string 
#' 
#' @param x numeric vector
#' @param y numeric vector
#' @param ci whether or not to print confidence interval and p-value
#' @param ... further options passed to `cor.test()` function
corLabel <- function(x, y, ci = FALSE, ...){
  cor_test <- cor.test(x, y, ...)
  
  # Get correlation, confidence interval and p-value
  r <- cor_test$estimate %>% prettyNum(nsmall = 2, digits = 2)
  r_ci <- cor_test$conf.int %>% prettyNum(nsmall = 2, digits = 2) %>% paste(collapse = ", ")
  r_p <- cor_test$p.value %>% prettyNum(digits = 2)
  
  if(ci){
    cor_label <- paste0("r = ", r, " [", r_ci, "]\np = ", r_p)
  } else {
    cor_label <- paste0("r = ", r)
  }
  return(cor_label)
}
