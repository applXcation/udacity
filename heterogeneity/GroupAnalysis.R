GroupAnalysis <- function(chopit.results, grouping.var, both = FALSE, thresh.outliers = 0) {
  # Function that calculates:
  #   the average y.s for each group
  #   the variance of y.s for each group
  #   the average thresholds for each group
  #   the space between first and last threshold in each group
  #
  # Args
  #   chopit.results: object of class Chopit, that contains the estimated parameters beta and gamma among others
  #   grouping.var: data frame containing the grouping variables
  #   both: should both the mean of Y* and its variance be included in the regressions ?
  #   thresh.outliers: theshold for the number of individuals per group (strictly) below wich values calculated from a given group are dropped 
  #                    from the analysis (i.e. when regressing a given function of the thresholds on a given function of the latent).
  #                    0 (default) if no observation dropped
  #
  # Returns an object of class "group.analysis", basically a list containing
  #   group.mean.y.star: estimated average for y.star.s in each group
  #   group.var.y.star: estimated variance for y.star.s in each group
  #   average.thresholds: number of groups * (K-1) matrix, whose i-th row is the vector containing each average thresholds in the i-th group
  #   average.average.thresholds: for each group, gives the average value of all the thresholds. This captures the shift in average threshold
  #                               from a group to another
  #   diff.average.thresholds: for each group, calculates the difference between the last and the first (non-infinites) thresholds. This captures
  #                            the shift in spaces between thresholds across groups
  #   group.count: the j-th row gives, for the j-th group, the number of individuals for which it was possible to estimate y.star.s[i],
  #                and the number of individuals for which it was possible to estimate the full set of tau[i, ]
  #   above.thresh: logical vector ; TRUE if obs in the group >=  thresh.outliers ; FALSE if not, NA cannot say
  
  # Sanity Check
  cat("GroupAnalysis: This program does not allow for sigma to vary across groups \n")
  if (!class(chopit.results) == "Chopit") cat("Warning: argument chopit.results is not of class 'Chopit'")
  if (!is.data.frame(grouping.var)) stop("GroupAnalysis: grouping.var should be a data.frame")
  if (length(chopit.results$coef$beta) != ncol(chopit.results$original.x$original.x.s) -1) cat("GroupAnalysis: in chopit.results, x.s et beta0 are  
                                                                                                not compatible")
  if (nrow(chopit.results$coef$gamma) != ncol(chopit.results$original.x$original.x.tau)) cat("GroupAnalysis: in chopit.results, x.s et beta0 are 
                                                                                                not  compatible")
    
  # Calculating y.star.s.hat and tau0.hat. tau0.hat is defined as the matrix where each row is an individual, and the columns are the 
  # (non-inifinite) thresholds
  beta0 <- chopit.results$coef$beta0
  gamma <- chopit.results$coef$gamma
  kK <- chopit.results$constants$kK
  x.s <- chopit.results$original.x$original.x.s[, -1]
  x.tau <- chopit.results$original.x$original.x.tau
  y.star.s.hat <- as.matrix(x.s) %*% t(t(beta0))
  tau0.hat <- as.matrix(x.tau) %*% gamma[, 1, drop = FALSE]
  if (kK >1) { 
    for (k in 2:(kK - 1)) tau0.hat <- cbind(tau0.hat, tau0.hat[, ncol(tau0.hat), drop = FALSE] + exp(as.matrix(x.tau) %*% 
    gamma[, k, drop = FALSE]))
  }
  
  # Estimating the first two moments for y* in each group
  group.mean.y.star <- tapply(y.star.s.hat, interaction(grouping.var), mean, na.rm = TRUE)
  group.var.y.star <- tapply(y.star.s.hat^2, interaction(grouping.var), mean, na.rm = TRUE) + 1 - group.mean.y.star^2

  # Average of each threshold by groups. On each column, we apply a function that calculates, for each group, the average thresholds.
  # These calculations result in a number of groups * (K-1) matrix, whose i-th row gives the average thresholds for the i-th group
  AverageThreshold <- function(tau.k, grouping.var) {
  	tapply(tau.k, interaction(grouping.var), mean, na.rm = TRUE)
  }
  average.thresholds <- apply(tau0.hat, 2, AverageThreshold, grouping.var = interaction(grouping.var))
  rownames(average.thresholds) <- paste("group", rownames(average.thresholds)) ; colnames(average.thresholds) <- paste("tau", 1:(kK-1))
  
  
  # Recovering the mean value of all the thresholds in each group, as well as the difference between the first and the last threshold
  average.average.thresholds <- apply(average.thresholds, 1, mean, na.rm = TRUE)
  diff.average.thresholds <- average.thresholds[, kK - 1] - average.thresholds[, 1]
  
  # Number of observations used to calculate these moments in each group
  group.count.y.star <- tapply(complete.cases(y.star.s.hat), interaction(grouping.var), sum)
  group.count.tau <- tapply(complete.cases(tau0.hat), interaction(grouping.var), sum)
  group.count <- cbind(group.count.y.star, group.count.tau) 
  rownames(group.count) <- paste("group", rownames(group.count)) ; colnames(group.count) <- c("Obs. for y.star", "Obs. for tau")
 
  # Regressions
  above.thresh <- c((group.count[, 1] >= thresh.outliers) & (group.count[, 2] >= thresh.outliers))  # logical vector ; TRUE if (weakly) more 
 																									# obs in the group than thresh.outliers ; 
 																									# FALSE if not, NA cannot say
  if (!both) {
  	lm.expectation <- lm(average.average.thresholds ~ group.mean.y.star, subset = above.thresh)
    lm.variance <- lm(diff.average.thresholds ~ group.var.y.star, subset = above.thresh)  
  }
  else {
  	lm.expectation <- lm(average.average.thresholds ~ group.mean.y.star + group.var.y.star, subset = above.thresh)
    lm.variance <- lm(diff.average.thresholds ~ group.var.y.star + group.mean.y.star, subset = above.thresh)  

  }
  
  # Returns
  group.analysis <- list(group.mean.y.star = group.mean.y.star, group.var.y.star = group.var.y.star, average.thresholds = average.thresholds, 
  average.average.thresholds = average.average.thresholds, diff.average.thresholds = diff.average.thresholds,
  group.count = group.count, lm.expectation = lm.expectation, lm.variance = lm.variance, thresh.outliers = thresh.outliers , 
  above.thresh = above.thresh)
  class(group.analysis) <- "group.analysis"
  return(group.analysis)
}



summary.group.analysis <- function(group.analysis, no.print = FALSE) {
  # Functions that plots the regression results for the two main analyses performed by GroupAnalysis:
  #   - regression of the groups' estimated expected values in y.star.s on the estimated expected value of the groups' average thresholds values
  #   - regression of the groups' estimated variance of y.star.s on the estimated expected value of the groups' space btwn last and first threshold
  # The data and regression lines are plotted, together with N, the number of observations, and the p-value for the slope of the reg line
  #
  # Args:
  #   group.analysis: object returner by the function GroupAnalysis()
  #   no.print: if true, no print in the summary (only graph)
  #
  # Return:
  #   Noting
  
  # Sanity Checks
  if (class(group.analysis) != "group.analysis") stop("plot.GroupAnalysis: group.analysis is not of class group.analysis")

  # Defining variables
  lm.expectation <- group.analysis$lm.expectation
  lm.variance <- group.analysis$lm.variance
  above.thresh <- group.analysis$above.thresh
  thresh.outliers <- group.analysis$thresh.outliers
  summary.expectation <- summary(lm.expectation)
  summary.variance <- summary(lm.variance)


  # Printing results
  if (!no.print) {
    cat("1°) Regression of the groups' estimated expected values in y.star.s on the estimated expected value of the groups' average thresholds 
         values")
    print(summary.expectation)
    cat("2°) Regression of the groups' estimated variance of y.star.s on the estimated expected value of the groups' space btwn last and first 
        threshold")
    print(summary.variance)
    cat("3°) Number of observations that can be used in each group in order to make the computations \n")
    print(group.analysis$group.count) ; cat("\n")
    cat("4°) Groups are dropped from the analysis if they contain less than", thresh.outliers, "observations \n")
  }
  
  # Plotting
  par(mfrow = c(1,2))
  #   First Graph
  plot(x = group.analysis$group.mean.y.star[above.thresh], y = group.analysis$average.average.thresholds[above.thresh], 
        xlab = "Group's average y*", ylab = "Average group's threshold",
       sub = paste("p-value:", round(summary.expectation$coefficients[2, 4], digits = 3), "; G =", length(summary.expectation$residuals)))
  abline(reg = lm.expectation)
  #   Second Graph
  plot(x = group.analysis$group.var.y.star[above.thresh], y = group.analysis$diff.average.thresholds[above.thresh], 
        xlab = "Group's variance in y*", ylab = "Group's diff. between first and last thresh.",
       sub = paste("p-value:", round(summary.variance$coefficients[2, 4], digits = 3), "; G =", length(summary.variance$residuals)))
  abline(reg = lm.variance)
}


