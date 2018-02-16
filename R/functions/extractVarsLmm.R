#' Title
#'
#' This function is used to extract variance components from a random slopes LMM.
#' It follows the methods in [Brommer (2013)](https://link.springer.com/article/10.1007/s00265-013-1603-9)
#'
#' @param model a lmerMod object
#' @param id_genotype the name of the variable that contains genotype IDs
#' @param overdisp the name of the variable with an overdispersion parameter 
#' (for GLMM models)
extractVarsLmm <- function(model, id_genotype = c("id_kover", "id_gwapp", "id_leyser"), 
                           overdisp = "id_unique"){
  
  # Get variance covariance matrices of random terms
  vcov_all <- VarCorr(model)
  
  # Get the variance-covariance matrix from line ID
  vcov_matrix <- vcov_all[[which(names(vcov_all) %in% id_genotype)]] %>% 
    as.data.frame %>% as.matrix()
  
  # Get the variance in slopes (plasticity)
  var_slopes <- vcov_matrix[2, 2]
  
  # Calculate correlation between LN (intercept) and the plasticity (slopes) 
  cor_ln_slopes <- vcov_matrix[1,2]/sqrt(prod(diag(vcov_matrix)))
  
  # Calculate the variance on HN from that variance-covariance matrix
  phi <- matrix(c(1, 1, 0, 1), nrow = 2) # design matrix
  
  # Between-environment matrix
  env_vcov_matrix <- phi %*% vcov_matrix %*% t(phi)
  
  # Extract the variance for LN, HN and the correlation between them
  var_ln <- env_vcov_matrix[1, 1]
  var_hn <- env_vcov_matrix[2, 2]
  
  cor_ln_hn <- env_vcov_matrix[1,2]/sqrt(prod(diag(env_vcov_matrix)))
  
  # Extract the residual variance
  if(family(model)$family == "gaussian"){
    var_res <- attr(VarCorr(model), "sc")^2
  } else if(family(model)$family == "poisson" & family(model)$link == "log"){
    
    # Calculate distribution-specific variance following table 2 of Nakagawa & Schielzeth 2010
    # Beta0 adapted from the function `r.squaredGLMM.merMod` from the MuMIn package
    beta0 <- model.matrix(model) %*% fixef(model) %>% mean
    var_dis <- log(1/exp(beta0)+1)
    
    # Add the over-dispersion variance to get total "residual" variance
    var_res <- var_dis + vcov_all[[overdisp]][1]
  }
  
  # Variance of the fixed effects
  var_fixed <- model.matrix(model) %*% fixef(model) %>% var
  
  # Extract the overall variance explained by the model (R^2)
  r2_m <- MuMIn::r.squaredGLMM(model)[1]
  r2_c <- MuMIn::r.squaredGLMM(model)[2]
  
  # Output a data.frame
  data.frame(var_ln = var_ln,
             var_hn = var_hn,
             var_gen = (var_ln + var_hn)/2,
             var_nitrate = var_fixed,
             var_plas = var_slopes,
             var_res = var_res,
             cor_ln_plas = cor_ln_slopes,
             cor_ln_hn = cor_ln_hn,
             r2_fixed = r2_m,
             r2_all = r2_c)
}
