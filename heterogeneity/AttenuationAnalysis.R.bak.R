AttenuationAnalysis <- function(chopit.results, indiv.calculus = TRUE, fix = FALSE, plot.exp = FALSE) {
  # Args
  #   indiv.calculus: should each individual's threshold be used for a counterfactual calculation
  #   fix.heterosk: should the conditional sd of y.star.s be kept constant across individuals, taking the average sd as a value ?
  # Return
  #   spearman: a list containing
  #     Spearman's correlation coeficient between E(y.star.s / x) and E(y.s / x), equal to 1 without heteroskedastictiy, but possibly smaller when 
  #     heteroskedasticity
  #     Spearman's correlation coeficient between V(y.star.s / x) and V(y.s / x), equal to 1 when x.beta is a constant, but possibly smaller else 
    
  # Sanity Check
  if (!class(chopit.results) == "Chopit") cat("AttenuationAnalysis: argument chopit.results is not of class 'Chopit'")
  
  # If plot is TRUE, then indiv.calculus has to be FALSE
  if (plot.exp) {
  	indiv.calculus <- FALSE
  	par(mfrow = c(1, 2))
  }
  
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
  fix.exp.x.kappa <- exp.x.kappa * 0 + mean(exp.x.kappa, na.rm = TRUE)  # cf. 'fix.heterosk' in the args description
  fix.x.beta <- matrix(mean(x.beta, na.rm = TRUE), nrow = length(x.beta))

  
  # Correlation of conditional expectations
  #   Computing the correlation between E(y.star.s/x) = x.beta and E(y.s/x) when thresholds vary across respondents
  corr.exp.baseline <- CorrEYstarEY(x.beta = x.beta, exp.x.kappa = switch(fix.num, exp.x.kappa, fix.exp.x.kappa), tau = tau, kK = kK,
 								    plot.exp = plot.exp, round = 1)
  #   Computing the correlation between E(y.star.s/x) = x.beta and E(y.s/x) when, for all resp, the thresholds are the average thresholds
  corr.exp.average <- CorrEYstarEY(x.beta = x.beta, exp.x.kappa = switch(fix.num, exp.x.kappa, fix.exp.x.kappa), tau = average.tau, kK = kK, 
  								    plot.exp = plot.exp, round = 2)  
  #   We iterate this calculus N times, by taking the thresholds of respondent i as the benchmark thresholds (we impose resp i thresholds to all other 
  #   resp)
  if (indiv.calculus) corr.exp.individuals <- apply(tau[complete.cases(tau), ], MARGIN = 1, FUN = function(taui, x.beta, exp.x.kappa, kK, plot.exp) {
  	taui <- matrix(1, nrow = length(y.s)) %x% matrix(taui, ncol = length(taui))
  	CorrEYstarEY(x.beta = x.beta, exp.x.kappa = exp.x.kappa, tau = taui, kK = kK, plot.exp = plot.exp)
  }, x.beta = x.beta, exp.x.kappa = switch(fix.num, exp.x.kappa, fix.exp.x.kappa), kK = kK, plot.exp = FALSE)
  else corr.exp.individuals <- NULL

  
  # Correlation of conditional variances
  #   Computing the correlation between V(y.star.s/x) = x.beta and V(y.s/x) when thresholds vary across respondents
  corr.var.baseline <- CorrVYstarVY(x.beta = switch(fix.num, x.beta, fix.x.beta), exp.x.kappa = exp.x.kappa, tau = tau, kK = kK)
  #   Computing the correlation between V(y.star.s/x) = x.beta and V(y.s/x) when, for all resp, the thresholds are the average thresholds
  corr.var.average <- CorrVYstarVY(x.beta = switch(fix.num, x.beta, fix.x.beta), exp.x.kappa = exp.x.kappa, tau = average.tau, kK = kK)
  if (indiv.calculus) corr.var.individuals <- apply(tau[complete.cases(tau), ], MARGIN = 1, FUN = function(taui, x.beta, exp.x.kappa, kK) {
  	taui <- matrix(1, nrow = length(y.s)) %x% matrix(taui, ncol = length(taui))
  	CorrVYstarVY(x.beta = x.beta, exp.x.kappa = exp.x.kappa, tau = taui, kK = kK)
  }, x.beta = switch(fix.num, x.beta, fix.x.beta), exp.x.kappa = exp.x.kappa, kK = kK)
  else corr.var.individuals <- NULL
  
  
  # Spearman correlations (see 'return' in description)
  spearman.exp <- SpearmanEYstarEY(x.beta = x.beta, exp.x.kappa = exp.x.kappa, tau = tau, kK = kK)
  spearman.var <- SpearmanVYstarVY(x.beta = x.beta, exp.x.kappa = exp.x.kappa, tau = tau, kK = kK)  

  
  # Return
  return(list(corr.exp = list(corr.baseline = corr.exp.baseline, corr.average = corr.exp.average, corr.individuals = corr.exp.individuals),
              corr.var = list(corr.baseline = corr.var.baseline, corr.average = corr.var.average, corr.individuals = corr.var.individuals),
              spearman = list(spearman.exp = spearman.exp, spearman.var = spearman.var)))
}



CorrXbetaYs <- function(x.beta, exp.x.kappa, tau, kK) {
  # Calculates the correlation between E(y.star.s / x) and y.s for a given x.beta, sigma(x), tau and kK
  # Not used in the AttenuationANalysis.
  proba.i <- ProbaI(x.beta, exp.x.kappa, tau, kK)
  k.proba.i <- t(apply(X = proba.i, MARGIN = 1, FUN = function(x, kK) x * (1 : kK), kK = kK))  # each proba is multiplied by the associated response
  k.sq.proba.i <- t(apply(X = proba.i, MARGIN = 1, FUN = function(x, kK) x * (1 : kK)^2, kK = kK))  # each proba is multiplied by the associated response^2
  sum.k.proba.i <- rowSums(k.proba.i)  # sum over all the proba for each respondent

  cov.xbeta.ys <- mean(x.beta * sum.k.proba.i, na.rm = TRUE) - mean(x.beta, na.rm = TRUE) * mean(sum.k.proba.i, na.rm = TRUE)  #cov E(y.star.s/x) and y.s
  var.xbeta <- var(x.beta, na.rm = TRUE)  # var(y.star.s/x) = var(x.beta/x)
  var.ys <- mean(k.sq.proba.i, na.rm = TRUE) - (mean(k.proba.i, na.rm = TRUE))^2
  corr.xbeta.ys <- cov.xbeta.ys / sqrt(var.xbeta * var.ys)	

  return(corr.xbeta.ys)
}



CorrEYstarEY <- function(x.beta, exp.x.kappa, tau, kK, plot.exp, round = 0) {
  # Calculates the correlation between E(y.star.s / x) and E(y.s / x) for a given x.beta, sigma(x), tau and kK
  proba.i <- ProbaI(x.beta, exp.x.kappa, tau, kK)
  ey <- proba.i %*% matrix(1:kK, nrow = kK)  # sum over all the proba for each respondent (E(y.s / x))
  corr.eystar.ey <- cor(ey, x.beta, use = "complete.obs")
  if (plot.exp & round == 1) plot(density(ey, na.rm = TRUE))
  if (plot.exp & round == 2) lines(density(ey, na.rm = TRUE), col = "red")  
  
  return(corr.eystar.ey)
}



CorrVYstarVY <- function(x.beta, exp.x.kappa, tau, kK) {
  # Calculates the correlation between V(y.star.s / x) and V(y.s / x) for a given x.beta, sigma(x), tau and kK
  proba.i <- ProbaI(x.beta, exp.x.kappa, tau, kK)
  sum.k.proba.i <- proba.i %*% matrix(1 : kK, nrow = kK)  # sum over all the proba for each respondent (E(y.s / x))
  sum.k.sq.proba.i <- proba.i %*% matrix((1 : kK) ^ 2, nrow = kK)  # sum over all the proba for each respondent (E(y.s / x))

  var.y.star.s <- exp.x.kappa ^ 2
  var.y.s <- sum.k.sq.proba.i - (sum.k.proba.i) ^ 2
  corr.var <- cor(var.y.star.s, var.y.s, use = "complete.obs")

  return(corr.var)
}


SpearmanEYstarEY <- function(x.beta, exp.x.kappa, tau, kK) {
  # Calculates the correlation between E(y.star.s / x) and E(y.s / x) for a given x.beta, sigma(x), tau and kK
  proba.i <- ProbaI(x.beta, exp.x.kappa, tau, kK)
  ey <- proba.i %*% matrix(1:kK, nrow = kK)  # sum over all the proba for each respondent (E(y.s / x))
  spearman.eystar.ey <- cor(ey, x.beta, use = "complete.obs", method = "spearman")
#  plot(x.beta, ey)
  
  return(spearman.eystar.ey)
}



SpearmanVYstarVY <- function(x.beta, exp.x.kappa, tau, kK) {
  # Calculates the correlation between V(y.star.s / x) and V(y.s / x) for a given x.beta, sigma(x), tau and kK
  proba.i <- ProbaI(x.beta, exp.x.kappa, tau, kK)
  sum.k.proba.i <- proba.i %*% matrix(1 : kK, nrow = kK)  # sum over all the proba for each respondent (E(y.s / x))
  sum.k.sq.proba.i <- proba.i %*% matrix((1 : kK) ^ 2, nrow = kK)  # sum over all the proba for each respondent (E(y.s / x))

  var.y.star.s <- exp.x.kappa ^ 2
  var.y.s <- sum.k.sq.proba.i - (sum.k.proba.i) ^ 2
  spearman.var <- cor(var.y.star.s, var.y.s, use = "complete.obs", method = "spearman")
#  plot(var.y.star.s, var.y.s)
  
  return(spearman.var)
}



ProbaI <- function(x.beta, exp.x.kappa, tau, kK) {  # for all individuals, predicts the proba of each response cditional on X, given x.beta, exp.x.kappa, 
											        # and tau. Checked.
  proba.i <- matrix(NA, nrow = nrow(x.beta), ncol = kK)												  
  for (k in 1:kK) {
    proba.i[, k] <- pnorm(1 / exp.x.kappa * (tau[, k + 1] - x.beta)) - pnorm(1 / exp.x.kappa * (tau[, k] - x.beta))
  }
  return(proba.i)
}



table.attenuation <- function(..., individual.cf = TRUE, moment) {
  # ... should be a serie of attenuation analyses
  # moment: either "expectation" or "variance". If expectation, Corr[E(Y*/X), E(Y/X)] is calculated. If variance, Corr[V(Y*/X), V(Y/X)] is calculated. 
  if (!(moment %in% c("expectation", "variance"))) stop("table.attenuation: 'moment' wrongly specified")
  
  
  attenuation.analyses <- list(...)
  
  
  attenuation.table <- sapply(attenuation.analyses, FUN = function(analysis) {
  	if (moment == "expectation") {
  	  baseline.corr <- analysis$corr.exp$corr.baseline
  	  cf.average.thresh.corr <- analysis$corr.exp$corr.average
  	  if (individual.cf) {
  	  average.cf.corr <- mean(analysis$corr.exp$corr.individuals, na.rm = TRUE)  # Average correlation across individual correlations
  	  	cf.attenuations <- (baseline.corr - analysis$corr.exp$corr.individuals) / analysis$corr.exp$corr.individuals  # Attenuation in each counterfactual computation
  	    average.cf.attenuation <- - mean(cf.attenuations, na.rm = TRUE)  # Average attenuation among cf.attenuations (equal to minus the percentage change
  	  																   # from baseline to counterfactual, because we talk about attenuation)
  	    prop.attenuation <- sum(cf.attenuations < 0, na.rm = TRUE) / sum(!is.na(cf.attenuations)) * 100  # % of cf where attenuation > 0
  	  }
  	  else cf.attenuations <- - (baseline.corr - cf.average.thresh.corr) / cf.average.thresh.corr
  	}
  	
  	if (moment == "variance") {
  	  baseline.corr <- analysis$corr.var$corr.baseline
  	  cf.average.thresh.corr <- analysis$corr.var$corr.average
  	  if (individual.cf) {
  	  	average.cf.corr <- mean(analysis$corr.var$corr.individuals, na.rm = TRUE)  # Average correlation across individual correlations
   	    cf.attenuations <- (baseline.corr - analysis$corr.var$corr.individuals) / analysis$corr.var$corr.individuals  # Attenuation in each counterfactual computation
   	    average.cf.attenuation <- - mean(cf.attenuations, na.rm = TRUE)  # Average attenuation among cf.attenuations (equal to minus the percentage change
  	  																   # from baseline to counterfactual, because we talk about attenuation)
  	    prop.attenuation <- sum(cf.attenuations < 0, na.rm = TRUE) / sum(!is.na(cf.attenuations)) * 100  # % of cf where attenuation > 0
  	  }
  	  else cf.attenuations <- - (baseline.corr - cf.average.thresh.corr) / cf.average.thresh.corr
  	}
  	
    if (individual.cf) return(c(baseline.corr = baseline.corr, average.cf.corr = average.cf.corr, average.cf.attenuation  = average.cf.attenuation, 
           prop.attenuation  = prop.attenuation))
    else return(c(baseline.corr = baseline.corr, cf.average.thresh.corr = cf.average.thresh.corr, cf.attenuations = cf.attenuations))
  	}
  )


  attenuation.table <- t(attenuation.table)
  return(attenuation.table)
}