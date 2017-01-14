IndividualAnalysisVar <- function(chopit.results, plot = FALSE, mfrow = TRUE, frac.displayed.points = 1, plot.bands = FALSE, boot.results = NULL, boot.thresh = 0.05, boot.p  = 10, boot.ci.method = "centiles") {
  # Function that regresses the individual conditional variances on the space between the last and the first thresholds, and then on each space between 
  # thresholds. Returns the correlation between these two types of objects as well. Plot the regression lines if required. Provides bootstrap confidence 
  # interval for the correlations, and confidence bands for the regression line, if boot.results is provided.
  # 
  #
  #
  # Args
  #   chopit.results: object of class Chopit, that contains the estimated parameters beta and gamma among others
  #   plot: TRUE if want the results to be plotted
  #   boot.results: indicates an object of class boot, that gives the results of a bootstrap procedure on the coef of IndividualAnalysisMean and
  #                 IndividualAnalysisVar
  #   boot.thresh: what is the tolerance threshold to be used if statistical inference is performed on the basis of boot.results ?
  #   boot.p: number on points on which the confidence band for the regression lines should be computed
  #
  # Returns (as a list)
  #   reglist: list containing the lm.summary() of each of the above-described regressions
  #   corrlist: correlations between the aboved-described objects
  #   ci.corrlist: confidence interval at boot.thresh for corrlist if boot.results has been provided
  #   cb.reglist: confidence intervals for the values cb.reglist$boot.p.values

    
  # Sanity Check
  if (!class(chopit.results) == "Chopit") cat("Warning: argument chopit.results is not of class 'Chopit'")
  if(!chopit.results$constants$heterosk) stop("GrouprAnalysis: chopit.results has not been estimated with sigma.s varying across respondents")
  if (frac.displayed.points > 1) stop("IndividualAnalysisMean: frac.displayed.points > 1")
  if (!(boot.ci.method %in% c("centiles", "pivotal"))) stop("IndividualAnalysisMean: boot.ci.method is none of 'centiles', 'pivotal'")
  
  # Calculating y.star.s.hat and tau0.hat. tau0.hat is defined as the matrix where each row is an individual, and the columns are the 
  # (non-inifinite) thresholds
  kappa0 <- chopit.results$coef$kappa0
  kappa <- c(0, kappa0)
  gamma <- chopit.results$coef$gamma
  kK <- chopit.results$constants$kK
  x.tau <- OriginalData(chopit.results, what = "x.tau")
  x.sigma <- OriginalData(chopit.results, what = "x.sigma")
  
  # Calculating sigma.sq and tau0.hat  
  sigma.sq <- (exp(x.sigma %*% t(t(kappa))))^2
  tau0.hat <- as.matrix(x.tau) %*% gamma[, 1, drop = FALSE]
  if (kK >1) { 
    for (k in 2:(kK - 1)) tau0.hat <- cbind(tau0.hat, tau0.hat[, ncol(tau0.hat), drop = FALSE] + exp(as.matrix(x.tau) %*% 
    gamma[, k, drop = FALSE]))
  }
  # cleaning sigma.sq and tau0.hat from infinite values (may happen in tau0.hat or in sigma.sq if some x*gamma or x*kappa goes wild in an exponential).
  inf.obs <- is.infinite(sigma.sq) | is.infinite(rowSums(tau0.hat, na.rm = TRUE))
  sigma.sq <- sigma.sq[!inf.obs]
  tau0.hat <- tau0.hat[!inf.obs, ]
  
  # Calculations
  last.to.first.tau0.hat <- tau0.hat[, kK - 1] - tau0.hat[, 1]
  diff.tau0.hat <- t(diff(t(tau0.hat)))
  reglist <- list()
  corrlist <- list()
  reglist$last.to.first <- lm(last.to.first.tau0.hat ~ sigma.sq)
  corrlist$last.to.first <- cor(last.to.first.tau0.hat, sigma.sq, use = "complete.obs")
  for (i in 1:(kK - 2)) {
    reglist[[i + 1]] <- lm(diff.tau0.hat[, i] ~ sigma.sq) ; names(reglist)[i + 1] <- paste("threshold", i+1, "-", "threshold", i)
    corrlist[[i + 1]] <- cor(diff.tau0.hat[, i], sigma.sq, use = "complete.obs")
  }
  summary.reglist <- lapply(reglist, summary)   

  
  # BOOTSTRAP
  # if the correlations and reg coef have been bootstrapped, then using these info to compute confidence interval and bands
  if (!is.null(boot.results)) {
  	B <- boot.results$R  # number of bootstrap iterations
    boot.values <- matrix(as.numeric(boot.results$t), ncol = 3*kK + 3*(kK - 1))  # as.numeric() is helpful in case something went wrong on a core (i.e. 
                                       						      # some elements in the original matrix could be characters, and so would then be the 
                                       						      # whole matrix
   boot.corr.individual.mean <- boot.values[, 1 : kK] ; count <- kK  # bootstraps on the correlations from IndividualAnalaysisMean
   boot.reg.coef.individual.mean <- boot.values[, (count + 1) : (count + 2 * kK)] ; count <- count + 2 * kK  # bootstraps on the reg coef from 
   																										# IndividualAnalaysisMean
   boot.corr.individual.var <- boot.values[, (count + 1) : (count + kK - 1)] ; count <- count + kK - 1  # bootstraps on the correlations from IndividualAnalaysisVar    
   boot.reg.coef.individual.var <- boot.values[, (count + 1) : (count + 2 * (kK - 1))]  # b. on the reg coef from IndividualAnalaysisVar
  }

  # Computing the bootstrapped confidence intervals and bands
  if (is.null(boot.results)) {
  	ci.corrlist <- NULL
  	confidence.bands <- NULL
  	cb.reglist <- NULL
  }
  if (!is.null(boot.results)) {
  	# Computing the confidence intervals for the correlation coefficients
  	ci.corrlist <- matrix(NA, nrow = kK - 1, ncol = 2)
  	for (i in 1:(kK - 1)) {  
  	ci.corrlist[i, ] <- BootCI(estimate = corrlist[[i]], bootstrap.sample = boot.corr.individual.var[, i], boot.thresh = boot.thresh,
  	         				   boot.ci.method = boot.ci.method)
  	}
  	# Computing the confidence bands for the regression line on boot.p points. These boot.p points are equally spaced on the 'x' axis, spanning from
  	# the first min to the max
  	cb.reglist <- list()
  	boot.p.values <- seq(from = min(sigma.sq, na.rm = TRUE), to = max(sigma.sq, na.rm = TRUE), length.out = boot.p)
  	for (i in 1:(kK - 1)) {  
  	  coef <- boot.reg.coef.individual.var[, c(2 * (i-1) + 1, 2 * (i-1) + 2)]  # Extracting the coefs of the i-th reg
      cb.reglist[[i]] <- sapply(boot.p.values, function(value, coef) {  # for each value in boot.p.values we calculate a + b*value with each coef and 
      																	# builds the CI form that
        predicted.y <- reglist[[i]]$coefficients[1] + reglist[[i]]$coefficients[2] * value  # predicted y with baseline coef 
        boot.predicted.y <- coef[, 1] + coef[, 2] * value
        cb <- BootCI(estimate = predicted.y, bootstrap.sample = boot.predicted.y, boot.thresh = boot.thresh,
  	         					 boot.ci.method = boot.ci.method)
        return(cb)
      }, coef = coef)
      names(cb.reglist)[i] <- paste("reg.", i, sep = "")
  	}  	
  cb.reglist$boot.p.values <- boot.p.values
  }

  
  
  # PLOT
  if (frac.displayed.points == 1) sample <- 1:length(sigma.sq)
  else sample <- sample(1:length(sigma.sq), size = trunc(frac.displayed.points * length(sigma.sq)), replace = FALSE)
  if (plot) {
  	if (mfrow) par(mfrow = c(1, kK - 1))
  	plot(sigma.sq, last.to.first.tau0.hat, main = "Space between last and 1st thresh", type = "n")
    if (!is.null(boot.results) & plot.bands) {  # if bootstrap results, we also plot the confidence bands obtained above
      polygon(c(cb.reglist$boot.p.values, rev(cb.reglist$boot.p.values)),
             c(cb.reglist[[1]][1, ],rev(cb.reglist[[1]][2, ])), col = "grey75", border = FALSE)
      lines(cb.reglist$boot.p.values, cb.reglist[[1]][1, ], col="red",lty=2)
      lines(cb.reglist$boot.p.values, cb.reglist[[1]][2, ], col="red",lty=2)
    }
    points(sigma.sq[sample], last.to.first.tau0.hat[sample], cex = .2)
    abline(reg = reglist$last.to.first)
    for (i in 1:(kK - 2)) {
      plot(sigma.sq, diff.tau0.hat[, i], main = paste("threshold", i+1, "-", "threshold", i), type = "n")
      if (!is.null(boot.results) & plot.bands) {  # if bootstrap results, we also plot the confidence bands obtained above
        polygon(c(cb.reglist$boot.p.values, rev(cb.reglist$boot.p.values)),
               c(cb.reglist[[i + 1]][1, ],rev(cb.reglist[[i + 1]][2, ])), col = "grey75", border = FALSE)
        lines(cb.reglist$boot.p.values, cb.reglist[[i + 1]][1, ], col="red",lty=2)
        lines(cb.reglist$boot.p.values, cb.reglist[[i + 1]][2, ], col="red",lty=2)
    }      
      points(sigma.sq[sample], diff.tau0.hat[sample, i], cex = .2)
      abline(reg = reglist[[i + 1]])
    }
  }


  # Return
  return(list(reglist = summary.reglist, corrlist = corrlist, ci.corrlist = ci.corrlist, cb.reglist = cb.reglist))
}


BootCI <- function(estimate, bootstrap.sample, boot.thresh, boot.ci.method) {
  if (boot.ci.method == "pivotal") boot.ci <- c(2 * estimate - quantile(bootstrap.sample, na.rm = TRUE, probs = 1 - boot.thresh / 2), 
                         			   2 * estimate - quantile(bootstrap.sample, na.rm = TRUE, probs =  boot.thresh / 2))
  if (boot.ci.method == "centiles") boot.ci <- c(quantile(bootstrap.sample, na.rm = TRUE, probs = boot.thresh / 2), 
   									    quantile(bootstrap.sample, na.rm = TRUE, probs = 1 - boot.thresh / 2))
  return(boot.ci)
}