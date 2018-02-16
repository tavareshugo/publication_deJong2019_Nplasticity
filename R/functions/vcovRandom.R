#' Obtain variance-covariance matrix from LMM produced by lme4
#'
#' @param model model to extract variance-covariance matrix
#' @param term random term
vcovRandomSlopes <- function(model, term = NULL){
  vcov_matrix <- VarCorr(model)
  
  if(length(vcov_matrix) != 1){
    stop("you must specificy a random term")
  } else {
    vcov_matrix <- vcov_matrix[[1]]
  }
  
  if(any(dim(vcov_matrix) != c(2, 2))){
    stop("this does not seem like a random slopes model")
  }
  
  # Variance-covariance matrix between intercept and slopes
  slopes_matrix <- vcov_matrix %>% as.data.frame() %>% as.matrix()
  
  return(slopes_matrix)
}


envCorr <- function(model){
  
  vcov_matrix <- vcovRandomSlopes(model)
  
  # design matrix
  phi <- matrix(c(1, 1, 0, 1), nrow = 2) 
  
  # Between-environment matrix
  between_matrix <- phi %*% vcov_matrix %*% t(phi)
  
  # Between-environment correlation
  between_matrix[1,2]/sqrt(prod(diag(between_matrix)))
}

intCorr <- function(model){
  vcov_matrix <- vcovRandomSlopes(model)
  
  vcov_matrix[1,2]/sqrt(prod(diag(vcov_matrix)))
}


