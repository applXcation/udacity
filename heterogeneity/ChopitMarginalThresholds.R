ChopitMarginalThresholds <- function(chopit.results, var.name) {
  # Args
  #   chopit.results: results from a Chopit() estimation 
  #   var.name: name of the variable of interest ; must belong to x.tau
  #
  # Return (as a list)
  #   Nota Bene, all the following objects contain two elements : marginal (the one described), and variance
  #   marginal.thresh:   N times (K-1) matrix where row i is the gradient of the thresholds with respect to var.name when x = x[i]. 
  #                      Column k thus is the first derivative of threshold k w.r.t var.name at point x = x[i].
  #   Remark: for discrete variables, a more appropriate approach would be to take the difference in lieu of the derivative. This should be
  #           implemented later.
  #   marginal.iaverage.thresh: N times 1 matrix, where each row i is the average is the thresholds for observation i
  #   marginal.idiff.thresh: same, but each row is the individual difference between last and first threshold
  #   average.marginal.iaverage.thresh: average in individual average thresholds
  #   average.marginal.idiff.thresh: average in individual space between thresholds 

  # Libraries
  library("numDeriv")

  # Sanity Check
  if (!class(chopit.results) == "Chopit") stop("ChopitMarginalThresholds: argument chopit.results is not of class 'Chopit'")

  # Retrieving the appropriate constants
  kK <- chopit.results$constants$kK
  kGammaNrow <- chopit.results$constants$kGammaNrow
  kBeta0Nrow <- chopit.results $constants$kBeta0Nrow
  beta0 <- chopit.results$coef$beta0
  gamma <- chopit.results$coef$gamma
  x.tau <- OriginalData(chopit.results, what = "x.tau")
  var.name.col <- match(var.name, colnames(x.tau))

  # Sanity Check
  if (!is.character(var.name)) stop("ChopitMarginal: var.name is not of type character")
  if (!(var.name %in% colnames(x.tau))) stop("ChopitMarginal: var.name does not belong to the regressor names")


  MarginalCalculations <- function(para, chopit.results) {
  	# Arg
  	#   para: c(beta0, gamma)
  	#   chopit.results, whose copy is going to be modified inside the function
  	#   All the other arguments obtained from the environment this function is defined in (i.e. ChopitMarginalThresholds)
  	#
  	# Return (ordered in a vector)
    #   marginal.thresh: N times (K-1) matrix where row i is the gradient of the thresholds with respect to var.name when x = x[i]. Column k thus
    #                      is the first derivative of threshold k w.r.t var.name at point x = x[i].
    #   Remark: for discrete variables, a more appropriate approach would be to take the difference in lieu of the derivative. This should be
    #           implemented later.
    #   marginal.iaverage.thresh: N times 1 matrix, where each row i is the average is the thresholds for observation i
    #   marginal.idiff.thresh: same, but each row is the individual difference between last and first threshold
    #   average.marginal.iaverage.thresh: average in individual average thresholds
    #   average.marginal.idiff.thresh: average in individual space between thresholds
    
    # Injecting para into chopit.results, in order to make the calculations with different values of para
      beta0 <- para[1:kBeta0Nrow]
      gamma <- matrix(para[(length(beta0) + 1) : (length(beta0) + (kK - 1) * kGammaNrow)], nrow = kGammaNrow)
      chopit.results$coef$beta0 <- beta0
      chopit.results$coef$gamma <- gamma

    # Computing marginal.thresh
    for (k in 1:(kK - 1)) {
    	if (k == 1) marginal.thresh <- matrix(gamma[var.name.col, 1], nrow = nrow(x.tau), ncol = 1)
    	else {
    	  deriv.tau <- marginal.thresh[, k - 1] + gamma[var.name.col, k] * exp(x.tau %*% gamma[, k, drop = FALSE])
    	  marginal.thresh <- cbind(marginal.thresh, deriv.tau)
    	}
    }

    # Individual marginal effects on the average threshold
    marginal.iaverage.thresh <- rowMeans(marginal.thresh)
    # Individual spaces between last and first thresholds
    marginal.idiff.thresh <- marginal.thresh[, kK - 1] - marginal.thresh[, 1]

    # Average in individual marginal effects on the average threshold
    average.marginal.iaverage.thresh <- mean(marginal.iaverage.thresh, na.rm = TRUE)
    # Average individual spaces
    average.marginal.idiff.thresh <- mean(marginal.idiff.thresh, na.rm = TRUE)
   
   # Return 
    return(c(marginal.thresh = marginal.thresh, 
             marginal.iaverage.thresh = marginal.iaverage.thresh, marginal.idiff.thresh = marginal.idiff.thresh, 
             average.marginal.iaverage.thresh = average.marginal.iaverage.thresh, average.marginal.idiff.thresh = average.marginal.idiff.thresh))
  }


  # Calculating the point-estimates of the above-mentionned quantitied, and parsing
  #   Parsing procedure: checked. Calculations: checked.
  marginal.calculations <- MarginalCalculations(c(beta0, gamma), chopit.results)
  #   count allows to know where we are in marginal.calculations
  marginal.thresh <- matrix(marginal.calculations[1 : (nrow(x.tau) * (kK - 1))], nrow = nrow(x.tau)) ; count <- length(marginal.thresh)
  marginal.iaverage.thresh <- marginal.calculations[(count + 1) : (count + nrow(x.tau))] ; count <- count + length(marginal.iaverage.thresh)
  marginal.idiff.thresh <- marginal.calculations[(count + 1) : (count + nrow(x.tau))] ; count <- count + length(marginal.idiff.thresh)
  average.marginal.iaverage.thresh <- marginal.calculations[count + 1] ; count <- count + 1
  average.marginal.idiff.thresh <- marginal.calculations[count + 1] ; count <- count + 1


  # Delta-Method
  # Calculating the jacobian of MarginalCalculations(para, chopit.results) wrt para, at the estimated parameters' value.
  # 'jacobian' is  a length(marginal.calculations) time number of parameters matrix
  jacobians <- jacobian(MarginalCalculations, x = c(beta0, gamma), method = "simple", chopit.results = chopit.results)
  # For each element of the vector returned by MarginalCalculations(estimated.para, chopit.results), calculating the variance through the delta
  # method. This requires performing a delta-method row by row.
  # A delta method could be performed in order to get the whole variance-covariane matrix for marginal.calculations, but the matrix would be
  # way too big.
  sigma <- chopit.results$var.cov[1 : length(c(beta0, gamma)), 1 : length(c(beta0, gamma))]
  variance <- function(jacobian, sigma) {
  	jacobian <- matrix(jacobian, ncol = length(jacobian))
  	variance <- jacobian %*% sigma %*% t(jacobian)
  }
  variances <- apply(X = jacobians, MARGIN = 1, FUN = variance, sigma = sigma)
  # Parsing the variance objects.
  var.marginal.thresh <- matrix(variances[1 : (nrow(x.tau) * (kK - 1))], nrow = nrow(x.tau)) ; count <- length(var.marginal.thresh)
  var.marginal.iaverage.thresh <- variances[(count + 1) : (count + nrow(x.tau))] ; count <- count + length(var.marginal.iaverage.thresh)
  var.marginal.idiff.thresh <- variances[(count + 1) : (count + nrow(x.tau))] ; count <- count + length(var.marginal.idiff.thresh)
  var.average.marginal.iaverage.thresh <- variances[count + 1] ; count <- count + 1
  var.average.marginal.idiff.thresh <- variances[count + 1] ; count <- count + 1
  
  
  # Return
  return(list(marginal.thresh = list(marginal = marginal.thresh, variance = var.marginal.thresh), 
         marginal.iaverage.thresh = list(marginal = marginal.iaverage.thresh, variance = var.marginal.iaverage.thresh), 
         marginal.idiff.thresh = list(marginal = marginal.idiff.thresh, variance = var.marginal.idiff.thresh), 
         average.marginal.iaverage.thresh = list(marginal = average.marginal.iaverage.thresh, variance = var.average.marginal.iaverage.thresh), 
         average.marginal.idiff.thresh = list(marginal = average.marginal.idiff.thresh, variance = var.average.marginal.idiff.thresh)))
}



# SYSTEMATIC ANALYSIS
SystematicImpactThresh <- function(chopit.results) {
  # SystematicImpactThresh: for each variable (except intercept) in chopit.results, displays (through print) the beta with its z-stat, together w/
  #                         the average of the variable's marginal impact on individual average threshold (and its z-stat), and on individual 
  #                         distance between last and first thresholds (and its z-stat). Does the same with the kappas if any.This object is of 
  #                         class  "systematic.marginal.thresh.work.dsbility"
  
  x.tau <- OriginalData(chopit.results, what = "x.tau")
  heterosk <- chopit.results$constants$heterosk
  # Requires the variables to be the same in x.s and x.tau (idem for x.sigma if heterosk has been allowed)
  analyses <- apply(X = matrix(colnames(x.tau[, -1]), ncol = 1),
                    MARGIN = 1,
                    FUN = function(x) {
                      marginal.thresholds <- ChopitMarginalThresholds(chopit.results, x)
                      
                      impact.average.thresh <- marginal.thresholds$average.marginal.iaverage.thresh$marginal
                      impact.average.thresh.var <- marginal.thresholds$average.marginal.iaverage.thresh$var
                      impact.average.thresh.t <- impact.average.thresh / sqrt(impact.average.thresh.var)
                      
                      impact.diff.thresh <- marginal.thresholds$average.marginal.idiff.thresh$marginal
                      impact.diff.thresh.var <- marginal.thresholds$average.marginal.idiff.thresh$variance
                      impact.diff.thresh.t <- impact.diff.thresh / sqrt(impact.diff.thresh.var)

                      latent.impact <- chopit.results$coef$beta0[x]
                      latent.impact.t <- latent.impact / sqrt(chopit.results$var.cov[x, x])

                      if (heterosk) {
                      	sd.impact <- chopit.results$coef$kappa0[paste("kappa", x)]
                        sd.impact.t <- sd.impact / sqrt(chopit.results$var.cov[paste("kappa", x), paste("kappa", x)])
                      }

                      if (!heterosk) return(c(latent.impact, latent.impact.t,
                                              impact.average.thresh, impact.average.thresh.t, 
                                              impact.diff.thresh, impact.diff.thresh.t))
                      if (heterosk)  return(c(latent.impact, latent.impact.t,
                                              impact.average.thresh, impact.average.thresh.t, 
                                              impact.diff.thresh, impact.diff.thresh.t,
                                              sd.impact, sd.impact.t))
                    }
                    )

  analyses <- t(analyses)
  analyses.results <- list(beta = list(value = analyses[, 1], z.stat = analyses[, 2]),
                           impact.av.thresh = list(value = analyses[, 3], z.stat = analyses[, 4]),
                           impact.diff.thresh = list(value = analyses[, 5], z.stat = analyses[, 6]))
  if (heterosk)  analyses.results$kappa <- list(value = analyses[, 7], z.stat = analyses[, 8])
  analyses.results$rownames <- colnames(x.tau[, -1])
  analyses.results$heterosk <- heterosk
  class(analyses.results) <- "systematic.impact.thresh"
  return(analyses.results)
}

print.systematic.impact.thresh <- function(systematic.marginal.thresh) {
  heterosk <- systematic.marginal.thresh$heterosk
  if (heterosk) analyses.to.print <- cbind(systematic.marginal.thresh$beta$value, systematic.marginal.thresh$beta$z.stat,
                                           systematic.marginal.thresh$kappa$value, systematic.marginal.thresh$kappa$z.stat,
                                           systematic.marginal.thresh$impact.av.thresh$value, systematic.marginal.thresh$impact.av.thresh$z.stat,
                                           systematic.marginal.thresh$impact.diff.thresh$value, 
                                           systematic.marginal.thresh$impact.diff.thresh$z.stat)
  if (!heterosk) analyses.to.print <- cbind(systematic.marginal.thresh$beta$value, systematic.marginal.thresh$beta$z.stat,
                                           systematic.marginal.thresh$impact.av.thresh$value, systematic.marginal.thresh$impact.av.thresh$z.stat,
                                           systematic.marginal.thresh$impact.diff.thresh$value, 
                                           systematic.marginal.thresh$impact.diff.thresh$z.stat)
  if (heterosk)  colnames(analyses.to.print) <- c("beta", "z-stat", "kappa", "z-stat", "impact av. thresh",  "z-stat", "impact diff. thresh",
                                                  "z-stat")
  if (!heterosk) colnames(analyses.to.print) <- c("beta", "z-stat", "impact av. thresh",  "z-stat", "impact diff. thresh", "z-stat")
  rownames(analyses.to.print) <- systematic.marginal.thresh$rownames
  print(analyses.to.print)
}


summary.systematic.impact.thresh <- function(systematic.marginal.thresh, alpha = 0.05) {
  # summary.systematic.impact.thresh: two columns table, that says
  #   1) among the variables whose impact on average threshold and on the expected value of the latent is significant, which percentage
  #      of the betas is significant and impact the latent in the same direction, which perc. is signif but impacts the latent to other way around
  #      Among the variables whose impact on average threshold is signif, which percentage of the betas is significant and impact the latent
  #      in the same direction, which perc. is signif but impacts the latent to other way around, and which percentage is insignificant.
  #      In a last row, says the percentage of variables that have signif impact on average thresh
  #   2) same thing, but with the impact on the difference between last and first threshold instead
  #   If heterosk (that if, if "kappa" is in the colnames of systematic.marginal.thresh), produces a second table, equivalent to this one, but 
  #   with the kappas replacing the betas
  # Args
  #   systematic.marginal.thresh: object of class 'systematic.marginal.thresh'
  #   alpha: level of the significance test (0.05 is the default)
  # Checked
  
  # Sanity Checks	
  if (class(systematic.marginal.thresh) != "systematic.impact.thresh") stop("summary.systematic.marginal.thresh: systematic.marginal.thresh is
                                      									       not of class 'systematic.impact.thresh'")
  
  kNBeta <- length(systematic.marginal.thresh$beta$value)
  
  signif.beta <- abs(systematic.marginal.thresh$beta$z.stat) > qnorm(1 - alpha / 2)  # TRUE if para affects significantly y.s
  signif.average.thresh <- abs(systematic.marginal.thresh$impact.av.thresh$z.stat) > qnorm(1 - alpha / 2)  # TRUE if para affects significantly average thresh
  signif.diff.thresh <- abs(systematic.marginal.thresh$impact.diff.thresh$z.stat) > qnorm(1 - alpha / 2)  # TRUE if para affects significantly diff btwn last and first thr.
  
  same.sign.average.thresh <- sign(systematic.marginal.thresh$beta$value) == sign(systematic.marginal.thresh$impact.av.thresh$value)
  same.sign.diff.thresh <- sign(systematic.marginal.thresh$beta$value) == sign(systematic.marginal.thresh$impact.diff.thresh$value)
  
  kSignifBeta <- sum(signif.beta, na.rm = TRUE)
  kSignifAverage <- sum(signif.average.thresh, na.rm = TRUE)  # Number of variables where impact on average thresh is signif
  kSignifBetaAverage <- sum(signif.average.thresh & signif.beta, na.rm = TRUE)
  kSignifDiff <- sum(signif.diff.thresh, na.rm = TRUE)  # Number of variables where impact on diff btwn last and first thresh is signif
  kSignifBetaDiff <- sum(signif.diff.thresh & signif.beta, na.rm = TRUE)
    
  # See description of these variables in the program description
  prop.all.signif.same.sign.average <- sum(signif.average.thresh & signif.beta  & same.sign.average.thresh, na.rm = TRUE) / kSignifBetaAverage
  prop.all.signif.contr.sign.average <- sum(signif.average.thresh & signif.beta & !same.sign.average.thresh, na.rm = TRUE) / kSignifBetaAverage  
  prop.impact.signif.same.sign.average <- sum(signif.average.thresh & signif.beta  & same.sign.average.thresh, na.rm = TRUE) / kSignifAverage
  prop.impact.signif.contr.sign.average <- sum(signif.average.thresh & signif.beta & !same.sign.average.thresh, na.rm = TRUE) / kSignifAverage  
  prop.not.signif.average <- sum(signif.average.thresh & !signif.beta, na.rm = TRUE) / kSignifAverage
  prop.signif.average <- kSignifAverage / kNBeta
  prop.all.signif.average <- sum(signif.average.thresh & signif.beta, na.rm = TRUE) / kNBeta
  card.all.signif.average <- sum(signif.average.thresh & signif.beta, na.rm = TRUE)

  prop.all.signif.same.sign.diff <- sum(signif.diff.thresh & signif.beta  & same.sign.diff.thresh, na.rm = TRUE) / kSignifBetaDiff
  prop.all.signif.contr.sign.diff <- sum(signif.diff.thresh & signif.beta & !same.sign.diff.thresh, na.rm = TRUE) / kSignifBetaDiff  
  prop.impact.signif.same.sign.diff <- sum(signif.diff.thresh & signif.beta  & same.sign.diff.thresh, na.rm = TRUE) / kSignifDiff
  prop.impact.signif.contr.sign.diff <- sum(signif.diff.thresh & signif.beta & !same.sign.diff.thresh, na.rm = TRUE) / kSignifDiff  
  prop.not.signif.diff <- sum(signif.diff.thresh & !signif.beta, na.rm = TRUE) / kSignifDiff
  prop.signif.diff <- kSignifDiff / kNBeta
  prop.all.signif.diff <- sum(signif.diff.thresh & signif.beta, na.rm = TRUE) / kNBeta
  card.all.signif.diff <- sum(signif.diff.thresh & signif.beta, na.rm = TRUE)

  prop.table.exp <- cbind(c(prop.all.signif.same.sign.average, prop.all.signif.contr.sign.average, prop.impact.signif.same.sign.average, 
                        prop.impact.signif.contr.sign.average, prop.not.signif.average, prop.signif.average, prop.all.signif.average, 
                        card.all.signif.average), 
  				      c(prop.all.signif.same.sign.diff, prop.all.signif.contr.sign.diff, prop.impact.signif.same.sign.diff, 
  				        prop.impact.signif.contr.sign.diff, prop.not.signif.diff, prop.signif.diff, prop.all.signif.diff, card.all.signif.diff))

  colnames(prop.table.exp) <- c("Average thresh.", "Diff. btwn last and 1st")
  rownames(prop.table.exp) <- c("Cditional impact and beta signif: same sign", "Cditional impact and beta signif: opp. sign", 
                            "Cditional impact signif: Impact and beta signif, same sign", "Cditional impact signif: Impact and beta signif, opp. sign", "Cditional impact signif: Impact signif, beta not signif", "Impact signif", "Prop. of impact and beta signif", "Number of impact and beta signif")
  print(prop.table.exp)
  
  if (systematic.marginal.thresh$heterosk == TRUE) {  # If heteroskedasticity is taken into account in the estimation
    
    signif.kappa <- abs(systematic.marginal.thresh$kappa$z.stat) > qnorm(1 - alpha / 2)
    
    kNKappa <- kNBeta
    kSignifKappaAverage <- sum(signif.average.thresh & signif.kappa, na.rm = TRUE)
    kSignifKappaDiff <- sum(signif.diff.thresh & signif.kappa, na.rm = TRUE)

    same.sign.average.thresh <- sign(systematic.marginal.thresh$kappa$value) == sign(systematic.marginal.thresh$impact.av.thresh$value)
    same.sign.diff.thresh <- sign(systematic.marginal.thresh$kappa$value) == sign(systematic.marginal.thresh$impact.diff.thresh$value)

    # See description of these variables in the program description
    prop.all.signif.same.sign.average <- sum(signif.average.thresh & signif.kappa  & same.sign.average.thresh, na.rm = TRUE) / kSignifKappaAverage
    prop.all.signif.contr.sign.average <- sum(signif.average.thresh & signif.kappa & !same.sign.average.thresh, na.rm = TRUE) / 
                                          kSignifKappaAverage  
    prop.impact.signif.same.sign.average <- sum(signif.average.thresh & signif.kappa  & same.sign.average.thresh, na.rm = TRUE) / kSignifAverage
    prop.impact.signif.contr.sign.average <- sum(signif.average.thresh & signif.kappa & !same.sign.average.thresh, na.rm = TRUE) / kSignifAverage  
    prop.not.signif.average <- sum(signif.average.thresh & !signif.kappa, na.rm = TRUE) / kSignifAverage
    prop.signif.average <- kSignifAverage / kNKappa
    prop.all.signif.average <- sum(signif.average.thresh & signif.kappa, na.rm = TRUE) / kNBeta
    card.all.signif.average <- sum(signif.average.thresh & signif.kappa, na.rm = TRUE)

    prop.all.signif.same.sign.diff <- sum(signif.diff.thresh & signif.kappa  & same.sign.diff.thresh, na.rm = TRUE) / kSignifKappaDiff
    prop.all.signif.contr.sign.diff <- sum(signif.diff.thresh & signif.kappa & !same.sign.diff.thresh, na.rm = TRUE) / kSignifKappaDiff  
    prop.impact.signif.same.sign.diff <- sum(signif.diff.thresh & signif.kappa  & same.sign.diff.thresh, na.rm = TRUE) / kSignifDiff
    prop.impact.signif.contr.sign.diff <- sum(signif.diff.thresh & signif.kappa & !same.sign.diff.thresh, na.rm = TRUE) / kSignifDiff  
    prop.not.signif.diff <- sum(signif.diff.thresh & !signif.kappa, na.rm = TRUE) / kSignifDiff
    prop.signif.diff <- kSignifDiff / kNKappa
    prop.all.signif.diff <- sum(signif.diff.thresh & signif.kappa, na.rm = TRUE) / kNBeta
    card.all.signif.diff <- sum(signif.diff.thresh & signif.kappa, na.rm = TRUE)
    
    prop.table.var <- cbind(c(prop.all.signif.same.sign.average, prop.all.signif.contr.sign.average, prop.impact.signif.same.sign.average, 
                        prop.impact.signif.contr.sign.average, prop.not.signif.average, prop.signif.average, prop.all.signif.average, 
                        card.all.signif.average), 
  				      c(prop.all.signif.same.sign.diff, prop.all.signif.contr.sign.diff, prop.impact.signif.same.sign.diff, 
  				        prop.impact.signif.contr.sign.diff, prop.not.signif.diff, prop.signif.diff, prop.all.signif.diff, card.all.signif.diff))

    colnames(prop.table.var) <- c("Average thresh.", "Diff. btwn last and 1st")
    rownames(prop.table.var) <- c("Cditional impact and kappa signif: same sign", "Cditional impact and kappa signif: opp. sign", 
                            "Cditional impact signif: Impact and kappa signif, same sign", "Cditional impact signif: Impact and kappa signif, opp. sign", "Cditional impact signif: Impact signif, kappa not signif", "Impact signif", "Prop. of impact and kappa signif", "Number of impact and kappa signif")
    print(prop.table.var)
  }
  else propr.table.var <- NULL
  
  return(list(prop.table.exp = prop.table.exp, prop.table.var = prop.table.var))
}

