#' Diagnostic plots from LMM object
#'
#' @param model model to produce diagnost plots for
plotLmmDiag <- function(model){
  if(!require(broom)) stop("Please install the package broom")
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


