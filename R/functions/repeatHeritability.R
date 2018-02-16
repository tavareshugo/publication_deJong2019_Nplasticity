#' Calculate broad sense heritability for a trait
#' 
#' This function calculates the broad-sense heritability for a trait, based 
#' on repeated measures of the same genotype (repeatability of individual 
#' measures). The variance is estimated from a linear mixed model using genotype 
#' as a random effect. It also provides an estimate of narrow sense heritability 
#' based on mean measurements of each individual.
#' 
#' @param phenotype a vector of phenotypic values.
#' @param genotype a vector with genotype IDs
#' @param plot_diag a logical defining whether to plot diagnostics of the linear 
#' model or not. Default: TRUE
#' @param ci a numeric between 0 and 1 defining the confidence interval to 
#' produce. The confidence interval is obtained by parametric bootstrap. 
#' Default: NULL
#' @param nboot number of bootstrap samples to compute CI
#' @param ncores number of cores to use when estimating the CI. Default: 1
#' 
#' @return a list with individual and average (per genotype) heritability
repeatHeritability <- function(phenotype, genotype,
                               plot_diag = TRUE, 
                               ci = NULL, nboot = 100, ncores = 1){
	
  # Check some things
  if(!require(lme4)) stop("Please install 'lme4' package")
  
	# Remove any NA values
  complete_obs <- !is.na(phenotype) & !is.na(genotype)
  phenotype <- phenotype[complete_obs]
  genotype <- genotype[complete_obs]
  
  # Make sure genotype is a factor
  genotype <- as.factor(genotype)
	
  # Make sure phenotype is numeric
  if(!is.numeric(phenotype)) stop("phenotype must be numeric. Is is: ", class(phenotype))
  
  # Fit linear mixed model with genotype as random effect
  lmm_fit <- lmer(phenotype ~ 1|genotype)
  
  # Plot diagnostics if requested
  if(plot_diag) .plotLmmDiag(lmm_fit)
  
  # Calculate heritabilities
  her <- .lmmHeritability(lmm_fit)

  # Calculate CI using parametric bootstrap
  if(!is.null(ci)){
    her <- bind_cols(her, 
                     .ciHeritability(lmm_fit, nboot = nboot, ncores = ncores, ci = ci))
  }
	
	# Result
  return(her)
	
}

.ciHeritability <- function(m, nboot = 100, ncores = 2, ci = 0.95){

  # Simulate phenotypes from the fitted model
  boot_simul <- simulate(m, nboot)
  
  # Estimate heritabilities for each simulated phenotype
  boot_fits <- parallel::mclapply(boot_simul, function(x){
    refit(m, newresp = x) %>% 
      .lmmHeritability()
  }, mc.cores = ncores)
  
  # Bind all estimates to a data.frame
  boot_fits <- bind_rows(boot_fits)
  
  # Return the Upper and Lower CI for each heritability
  boot_fits %>% 
    summarise_all(funs(lo = quantile(., c((1-ci)/2)),
                       hi = quantile(., ci + (1-ci)/2)))
}


.lmmHeritability <- function(lmm_fit){
  # Get genotype from model
  genotype <- lmm_fit@frame[,2]
  
  # Extract explained and non-explained variance
  model_variance <- VarCorr(lmm_fit) %>% as.data.frame %>% select(vcov)
  
  # Individual heritability, following Box et. al (2015) Curr Biol 25:194
  indiv_her <- model_variance[1,]/sum(model_variance)
  
  # Calculate genotype heritability, following Box et. al (2015) Curr Biol 25:194
  n_lines <- length(unique(genotype))
  n_reps <- table(genotype) %>% as.vector()
  
  if(n_lines != length(n_reps)) warning("Something's wrong! Check bug!")
  
  genotype_her <- vector(mode = "numeric", length = n_lines)
  
  for(i in seq_along(n_reps)){
    genotype_her[i] <- model_variance[1,]/(model_variance[1,] + model_variance[2,]/n_reps[i])
  }
  genotype_her <- sum(genotype_her)/n_lines
  
  return(data.frame(her_i = indiv_her, her_n = genotype_her))
}


.plotLmmDiag <- function(model){
  # "Augment" the model to obtain fitted and residual values
  model <- augment(model)
  
  # Setup graphical device
  par(mfrow = c(2, 2))
  
  # Q-Q plot
  qqnorm(model$.resid, main = "Residuals normal Q-Q plot"); 
  qqline(model$.resid, col = "red3", lty = 2)
  
  # Fitted vs Residuals
  scatter.smooth(model$.fitted, model$.resid, main = "Fitted vs. Residuals", 
                 xlab = "Fitted values", ylab = "Residuals",
                 lpars = list(col = "royalblue", lwd = 3))
  abline(h = 0, col = "red", lty = 2)
  
  # Scale-location plot
  scatter.smooth(model$.fitted, sqrt(abs(model$.resid)),
                 main = "Scale-location plot", 
                 xlab = "Fitted values", ylab = expression(sqrt(abs(Residuals))),
                 lpars = list(col = "royalblue", lwd = 3))
}
