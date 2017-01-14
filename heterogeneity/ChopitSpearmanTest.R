ChopitSpearmanTest <- function(chopit.results) {
  # Returns:
  #   Spearman's correlation coeficient between E(y.star.s / x) and E(y.s / x), equal to 1 without heteroskedastictiy, but possibly smaller when 
  #   heteroskedasticity
  #   Spearman's correlation coeficient between V(y.star.s / x) and V(y.s / x), equal to 1 when x.beta is a constant, but possibly smaller else 
    
  # Sanity Check
  if (!class(chopit.results) == "Chopit") cat("AttenuationAnalysis: argument chopit.results is not of class 'Chopit'")
  
  # Retrieving required objects
  fix.num <- 1 + 1 * fix  # equal to 1 if fix is FALSE, and to 2 if fix is TRUE
  heterosk <- chopit.results$constants$heterosk
  beta0 <- chopit.results$coef$beta0
  beta <- c(0, beta0)
  gamma <- chopit.results$coef$gamma
  kK <- chopit.results$constants$kK
  y.s <- OriginalData(chopit.results, what = "y.s")
  x.s <- OriginalData(chopit.results, what = "x.s")
  x.tau <- OriginalData(chopit.results, what = "x.tau")
  if (heterosk) {
    kappa0 <- chopit.results$coef$kappa0
    kappa <- c(0, kappa0)
    x.sigma <- OriginalData(chopit.results, what = "x.sigma")
  }
  if (!heterosk) {  # If heterosk was not allowed in chopit.results, we specify kappa0 = 0, and x.sigma = 0, so that both objects can be included w/out
  	                # loss of generality in all the calculations
    kappa0 <- kappa <- 0
    x.sigma <- matrix(0, ncol = 1, nrow = nrow(x.s))	
  }
  
  # Computing the thresholds, x.beta, x.kappa and exp.x.kappa
  tau0 <- as.matrix(x.tau) %*% gamma[, 1, drop = FALSE]
  if (kK >1) { 
    for (k in 2:(kK - 1)) tau0 <- cbind(tau0, tau0[, ncol(tau0), drop = FALSE] + exp(as.matrix(x.tau) %*% 
    gamma[, k, drop = FALSE]))
  }
  tau <- cbind(-Inf, tau0, + Inf)
  average.tau0 <- t(colMeans(tau0, na.rm = TRUE))
  average.tau0 <- matrix(1, nrow = length(y.s)) %x% average.tau0  # Replicating the average tau0 for each individual through Kronecker product
  average.tau <- cbind(-Inf, average.tau0, +Inf)
  x.beta <- x.s %*% beta
  x.kappa <- x.sigma %*% kappa
  exp.x.kappa <- exp(x.kappa)
 
  cor(x.beta, ) 
}