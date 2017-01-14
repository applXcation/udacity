# NOTA BENE: en l'état, ne fonctionne qu'avec les variables continues.
# TODO: intégrer une gestion des variables discrètes
library("numDeriv")

ChopitAverageMarginalSpan <- function(chopit.results, var.name, kNPoints) {
  # ChopitMarginalSpan: calculates the EXPECTED (w.r.t. k) derivative of P(y.s >= k / bar_x[-var.name], var.name.value, beta, gamma) 
  #                     w.r.t several values of 
  #                     var.name.value (with bar_x being the average values of the regressors, and bar_x[-var.name] being the same but w/out
  #                     regressor var.name). 
  #                     i.e. calculates sum_k P(y.s == k) P(y.s >= k / bar_x[-var.name], var.name.value, beta, gamma) for each var.name.value.
  #                     var.name.value takes kNPoints values, spanning from the first percentile of var.name to its last
  #                     percentile, and being equally spaced.
  #
  # Args
  #   chopit.results: object of class "Chopit". The x.tau and x.s should be the same for the marginal effects to be calculated.
  #   var.name: name of the variable for which marginal impacts should be calculated
  #   response: marginal effects are calculated on P(y.s >= response / x). Should be between 2 (proba equal to 1 if response = 1) and the highest
  #             response category
  #   kNPoints: number of points on which the marginal effects are going to be calculated
  #
  # Returns
  #   marginal.effects: list containing
  #     total.marginal: a list containing
  #       marginal: a vector of of which each element s is expected (w.r.t. response) the expected marginal impacts of variable in column 
  #               var.name.col in values on  P(y.s >= / values[i, s], beta, gamma),
  #               i.e. sum_response P(y.s == k) * P(y.s >= k / values[i, s], beta, gamma)
  #       variance: variance for these v alues
  #       CI: CI for these values
  #   latent.marginal: same thing, but considering only the marginal impact of the variable through the change in y.star.s
  #   tau.marginal: same thing, but considering only the marginal impact of the variable through the change in thresholds
  #   values.var: values of var.name on which the marginal effects are calculated  
  #   var.name: name of the variable w.r.t. which the marginal effects are calculated
  
  # Sanity Checks
  if (!class(chopit.results) == "Chopit") cat("ChopitAverageMarginalSpan: argument chopit.results is not of class 'Chopit' \n")
  if (!is.character(var.name)) stop("ChopitMarginal: var.name is not of type character")
  if (!(var.name %in% colnames(chopit.results$original.x$original.x.s))) stop("ChopitMarginal: var.name does not belong to the regressor names")
  if (!all(chopit.results$original.x$original.x.tau == chopit.results$original.x$original.x.s, na.rm = TRUE)) stop("ChopitMarginal: in
     chopit.results, x.tau != x.s")


  # Making the calculations: each row contain the results for a given value of x[var.name], with each column being the result for a given 
  # response, going from k = 2 to k = kK. Checked.
  kK <- chopit.results$constants$kK
  kGammaNrow <- chopit.results$constants$kGammaNrow
  marginal <- sapply(X = 2:(kK), FUN = ChopitMarginalSpanTotal, chopit.results, var.name, kNPoints)
  latent.marginal <- sapply(X = 2:(kK), FUN = ChopitMarginalSpanLatent, chopit.results, var.name, kNPoints)
  tau.marginal <- sapply(X = 2:(kK), FUN = ChopitMarginalSpanTau, chopit.results, var.name, kNPoints)
  values.var <- ChopitMarginalSpanValuesVar(2, chopit.results, var.name, kNPoints)  # The values are the same for all k
  
    
  # Calculating the expected value in each row. Checked.
  prop <- prop.table(table(chopit.results$original.y.s))
  prop <- as.matrix(prop, nrow = ncol(marginal))
  average.marginal <- marginal %*% prop[-1]
  average.latent.marginal <- latent.marginal %*% prop[-1]
  average.tau.marginal <- tau.marginal %*% prop[-1]


  # Variance for each element in marginal.effect (calculation through numerical delta method)
  #   First, we defined a function used to calculate the variance through the delta method
  MarginalEffectsJacobians <- function(para, chopit.results, var.name, kNPoints) {
  	  # MarginalEffectsJacobians : Calculates the Average Marginal Impact of var.name on a span of kNPoints, just like ChopitAverageMarginalSpan.
  	  #                            The main difference is that this function use para as the vector of paramters. To be used to calculate
  	  #                            jacobians. Checked.
  	  # Args and returns
  	  #   Same as ChopitAverageMarginalSpan, except that it takes para as a first argument. Vector of parameters at which the marginal effects
  	  #   should be calculated.

	  kK <- chopit.results$constants$kK
      kGammaNrow <- chopit.results$constants$kGammaNrow 
  	  
      # Sanity Checks
      if (length(para) != (ncol(chopit.results$original.x$original.x.s) - 1 + (kK - 1) * chopit.results$constants$kGammaNrow)) 
          stop("marginal.effects.jacobians: lentgh of parameter not conform")

      # Injecting para into chopit.results, in order to make the calculations with different values of para

      beta0 <- para[1:(ncol(chopit.results$original.x$original.x.s) - 1)]
      gamma <- matrix(para[(length(beta0) + 1) : (length(beta0) + (kK - 1) * kGammaNrow)], nrow = kGammaNrow)
      
      chopit.results$coef$beta0 <- beta0
      chopit.results$coef$gamma <- gamma
      
  	  # Making the calculations: each row contain the results for a given value of x[var.name], with each column being the result for a given 
      # response, going from k = 2 to k = kK      
	  marginal <- sapply(X = 2:(kK), FUN = ChopitMarginalSpanTotal, chopit.results, var.name, kNPoints)
	  latent.marginal <- sapply(X = 2:(kK), FUN = ChopitMarginalSpanLatent, chopit.results, var.name, kNPoints)
	  tau.marginal <- sapply(X = 2:(kK), FUN = ChopitMarginalSpanTau, chopit.results, var.name, kNPoints)
	  values.var <- ChopitMarginalSpanValuesVar(2, chopit.results, var.name, kNPoints)  # The values are the same for all k
  
	  # Calculating the expected value in each row
	  prop <- prop.table(table(chopit.results$original.y.s))
	  prop <- as.matrix(prop, nrow = ncol(marginal))
	  average.marginal <- marginal %*% prop[-1]
	  average.latent.marginal <- latent.marginal %*% prop[-1]
	  average.tau.marginal <- tau.marginal %*% prop[-1]
	  
	  return(c(average.marginal, average.latent.marginal, average.tau.marginal))
  }
  
  marginal.effects.jacobians <- jacobian(MarginalEffectsJacobians, c(chopit.results$coef$beta0, chopit.results$coef$gamma), method = "simple",
                                         chopit.results = chopit.results, var.name = var.name, kNPoints = kNPoints)
  
  #   Extracting the variance-covariance matrix associated to the vector c(beta0, gamma)
  kNpara <- length(chopit.results$coef$beta0) + (kK - 1) * kGammaNrow
  sigma <- chopit.results$var.cov[1 : kNpara, 1 : kNpara]

  #   Applying the delta method on each lign of marginal.effects.jacobians
  variances <- apply(marginal.effects.jacobians, MARGIN = 1, FUN = function(marginal.effects.jacobian, sigma) {
  	       jacob <- matrix(marginal.effects.jacobian, ncol = length(marginal.effects.jacobian))
  	       variance <- jacob %*% sigma %*% t(jacob)
  	       return(variance)
        }, sigma = sigma)
  marginal.effects.variances <- list(marginal.var = variances[1: kNPoints], latent.marginal.var = variances[(kNPoints + 1): (2 * kNPoints)],
                                     tau.marginal.var = variances[(2 * kNPoints + 1): (3 * kNPoints)])
  
  
  # Calculating the confidence intervals
  marginal.effects.CI <- list(
    marginal.CI = cbind(average.marginal - qnorm(0.975) * sqrt( marginal.effects.variances$marginal.var), 
                        average.marginal + qnorm(0.975) * sqrt( marginal.effects.variances$marginal.var)),
    latent.marginal.CI = cbind(average.latent.marginal - qnorm(0.975) * sqrt( marginal.effects.variances$latent.marginal.var), 
                              average.latent.marginal + qnorm(0.975) * sqrt( marginal.effects.variances$latent.marginal.var)),
    tau.marginal.CI = cbind(average.tau.marginal - qnorm(0.975) * sqrt( marginal.effects.variances$tau.marginal.var), 
                            average.tau.marginal + qnorm(0.975) * sqrt( marginal.effects.variances$tau.marginal.var))
  )
  
  
  # Return
  marginal.effects <- list(total.marginal = list(marginal = average.marginal, variance = marginal.effects.variances$marginal.var, 
                                                 CI = marginal.effects.CI$marginal.CI), 
                           latent.marginal = list(marginal = average.latent.marginal, variance = marginal.effects.variances$latent.marginal.var, 
                                                  CI = marginal.effects.CI$latent.marginal.CI), 
                           tau.marginal = list(marginal = average.tau.marginal, variance = marginal.effects.variances$tau.marginal.var, 
                                               CI = marginal.effects.CI$tau.marginal.CI),
                           values.var = values.var,
                           var.name = var.name)

  class(marginal.effects) <- "chopit.marginal.span"
  return(marginal.effects)
}



plot.chopit.marginal.span <- function(chopit.marginal.span) {
  # Plot marginal effects (total, through the latent, and through the thresholds), together with the CI for the last two types of effects
  
  # Sanity Checks
  if (class(chopit.marginal.span) != "chopit.marginal.span") stop("plot.chopit.marginal.span: object chopit.marginal.span is not of class
                                                                   chopit.marginal.span")

  # Creating shortcut objects (helps to make the code easier to read)
  total.marginal <- chopit.marginal.span$total.marginal
  latent.marginal <- chopit.marginal.span$latent.marginal
  tau.marginal <- chopit.marginal.span$tau.marginal
  values.var <- chopit.marginal.span$values.var
  var.name <- chopit.marginal.span$var.name
  
  # Putting together all marginal impacts and CIs to define the range of the plots y-axes                                                                 
  all.values <- c(total.marginal$marginal, total.marginal$CI, 
                  latent.marginal$marginal, latent.marginal$CI, 
                  tau.marginal$marginal, tau.marginal$CI)
                  
  # Plotting
  #   Plotting the total marginal impact
  plot(values.var, total.marginal$marginal, type = "l", 
       ylim = c(min(all.values, 0), max(all.values, 0)),
       main = var.name)
  #   Plotting the partial, through-the-latent, impact
  lines(values.var, latent.marginal$marginal, col = "blue")  # lty = "longdash"
  lines(values.var, latent.marginal$CI[, 1], col = "blue", lty = "dotdash")
  lines(values.var, latent.marginal$CI[, 2], col = "blue", lty = "dotdash")
      
  #   Plotting the partial, through-the-thresholds, impact  
  lines(values.var, tau.marginal$marginal, col = "red")  # lty = "longdash"
  lines(values.var, tau.marginal$CI[, 1], col = "red", lty = "dotdash")
  lines(values.var, tau.marginal$CI[, 2], col = "red", lty = "dotdash")

  # lines(chopit.marginal.span$marginal.effects$values.var, chopit.marginal.span$marginal.effects$latent.marginal, lty = "longdash")
  # lines(chopit.marginal.span$marginal.effects$values.var, chopit.marginal.span$marginal.effects$tau.marginal, lty = "dotted")
  abline(h = 0)
}



ChopitMarginalSpan <- function(response , chopit.results, var.name, kNPoints) {
  # ChopitMarginalSpan: calculates the derivative of P(y.s >= response / bar_x[-var.name], var.name.value, beta, gamma) w.r.t several values of 
  #                     var.name.value (with bar_x being the average values of the regressors, and bar_x[-var.name] being the same but w/out
  #                     regressor var.name). var.name.value takes kNPoints values, spanning from the first percentile of var.name to its last
  #                     percentile, and being equally spaced.
  #
  # Args
  #   chopit.results: object of class "Chopit". The x.tau and x.s should be the same for the marginal effects to be calculated.
  #   var.name: name of the variable for which marginal impacts should be calculated
  #   response: marginal effects are calculated on P(y.s >= response / x). Should be between 2 (proba equal to 1 if response = 1) and the highest
  #             response category
  #   kNPoints: number of points on which the marginal effects are going to be calculated
  #
  # Returns
  #   marginal.effects: list containing
  #     marginal: a vector of of which each element s is the marginal impacts of variable in column var.name.col in values on 
  #               P(y.s >= response / values[i, s], beta, gamma)
  #     latent.marginal: same thing, but considering only the marginal impact of the variable through the change in y.star.s
  #     tau.marginal: same thing, but considering only the marginal impact of the variable through the change in thresholds
  #     values.var: values of var.name on which the marginal effects are calculated
  
  # Sanity Check
  if (!class(chopit.results) == "Chopit") cat("ChopitMarginalSpan: argument chopit.results is not of class 'Chopit' \n")
  if (!is.character(var.name)) stop("ChopitMarginalSpan: var.name is not of type character")
  if (!(var.name %in% colnames(chopit.results$original.x$original.x.s))) stop("ChopitMarginalSpan: var.name does not belong to the regressor 
      names")
  if (!all(chopit.results$original.x$original.x.tau == chopit.results$original.x$original.x.s, na.rm = TRUE)) stop("ChopitMarginalSpan: in
     chopit.results, x.tau != x.s")
  if (response < 2) stop("ChopitMarginalSpan: response smaller than 2")   
  if (response > chopit.results$constants$kK) stop("ChopitMarginalSpan: response biggher than kK from chopit.results")
  if (kNPoints < 2) stop("ChopitMarginalSpan: kNPoints smaller than 2")
  
  # Locating the column of var.name
  var.name.col <- match(var.name, colnames(chopit.results$original.x$original.x.s))
  
  # Calculating the average values of all the regressors, including the intercept.
  x.average <- colMeans(chopit.results$original.x$original.x.s, na.rm = TRUE)
  
  # Calculating the set of values on which the marginal probabilities are going to be calculated
  var.name.quantiles <- quantile(chopit.results$original.x$original.x.s[, var.name.col], probs = c(0.01, 0.99), na.rm = TRUE)
  var.name.values <- seq(from = var.name.quantiles[1], to = var.name.quantiles[2], length.out = kNPoints)
  # each row of values is a value of x on which marginal impact should be calculated. All the variables are equal to their average, except of 
  # variable var.name, which takes one value of var.name.values in each row
  values <- t(as.matrix(unlist(replicate(length(var.name.values), x.average)), nrow = length(var.name.values), byrow = TRUE))
  values[, var.name.col] <- var.name.values 

  # Retrieving the right coefficients
  beta0 <- chopit.results$coef$beta0
  beta <- c(0, chopit.results$coef$beta0)
  gamma <- chopit.results$coef$gamma
  kK <- chopit.results$constants$kK

  # Calculate the marginal effects
  marginal.effects <- ChopitMarginal(values, var.name.col, response, beta, gamma)  
  
  # Return
  class(marginal.effects) <- "chopit.marginal.span"
  return(marginal.effects)
}



ChopitMarginalSpanTotal <- function(response, chopit.results, var.name, kNPoints) {
  # Same as ChopitMarginalSpan, but returns only the total marginal effect ("marginal")
  marginal.total <- ChopitMarginalSpan(response, chopit.results, var.name, kNPoints)$marginal
  return(c(marginal.total))
}
ChopitMarginalSpanLatent <- function(response , chopit.results, var.name, kNPoints) {
  # Same as ChopitMarginalSpan, but returns only the through-the-latent marginal effect ("latent.marginal")
  latent.marginal <- ChopitMarginalSpan(response , chopit.results, var.name, kNPoints)$latent.marginal
  return(c(latent.marginal))
}
ChopitMarginalSpanTau <- function(response , chopit.results, var.name, kNPoints) {
  # Same as ChopitMarginalSpan, but returns only the through-the-thresholds marginal effect ("tau.marginal")
  tau.marginal <- ChopitMarginalSpan(response , chopit.results, var.name, kNPoints)$tau.marginal
  return(c(tau.marginal))
}
ChopitMarginalSpanValuesVar <- function(response , chopit.results, var.name, kNPoints) {
  # Same as ChopitMarginalSpan, but returns only values
  values <- ChopitMarginalSpan(response , chopit.results, var.name, kNPoints)$values
  return(c(values))
}


ChopitMarginal <- function(values, var.name.col, response, beta, gamma) {
  # ChopitMarginal: for each row i of values, calculates the marginal impacts of the variable associated to the var.name.col column of values,
  #                 i.e. the derivative of P(y.s >= response / values[i, ], beta, gamma) w.r.t values[i, var.name.col].
  #                 Remark: "values" is used both for the equation of y.star.s and for the threshold equations, that is:
  #                         y.star.s = values %*% beta
  #                         threshold_1 = values %*% gamma etc.
  #                         Thus, given the usual normalization of y.star.s not containing an intercept, beta[1] is expected to be equal to 0.
  #
  # Args
  #   values: matrix of which each row is a vector of values x = x.s = x.tau on which the marginal effect should be calculated.
  #           Remark: it is assumed that the variables entering the y.star.s equation are the same as the variables entering the threshold 
  #           equations, except for the constant.
  #   var.name.col: column number of the variable for which the marginal effect is calculated
  #   response: response as in P(y.s >= response / values[i, ], beta, gamma)
  #   beta: vector beta (beta[1] is expected to be 0, cf. above)
  #   gamma: matrix that cbinds gamma.1, ..., gamma.K-1, with gamma.k a column vector whose number of rows matches the number of columns in x
  #
  # Returns
  #   marginal.effects: list containing
  #     marginal: a vector of of which each element s is the marginal impacts of variable in column var.name.col in values on 
  #               P(y.s >= response / values[i, s], beta, gamma)
  #     latent.marginal: same thing, but considering only the marginal impact of the variable through the change in y.star.s
  #     tau.marginal: same thing, but considering only the marginal impact of the variable through the change in thresholds
  #     values.var: values of var.name on which the marginal effects are calculated
  
  # Retrieving the right constants
  kK <- ncol(gamma) + 1
  kNPoints <- nrow(values)
  
  # Sanity Check
  if (!is.matrix(values)) stop("ChopitMarginal: values is not a matrix")
  if (response < 2) stop("ChopitMarginal: response smaller than 2")   
  if (response > kK) stop("ChopitMarginal: response biggher than kK from chopit.results")
  if (beta[1] != 0) cat("WARNING: ChopitMarginal: the first term of the vector beta is not a 0. \n")
  if (!is.matrix(gamma)) stop("ChopitMarginal: gamma is not a matrix")
  if (nrow(gamma) != ncol(values)) stop("ChopitMarginal: number of rows in gamma does not match the number of columns in values")

  # Calculating the tau for each row of values  
  tau0 <- values %*% gamma[, 1, drop = FALSE]
  if (kK >1) { 
    for (k in 2:(kK - 1)) tau0 <- cbind(tau0, tau0[, ncol(tau0), drop = FALSE] + exp(values %*% gamma[, k, drop = FALSE]))
  }
    
  # Calculating the through-y.star.s variable's marginal impact on P(y.s >= response) for each row in values
  latent.marginal <- beta[var.name.col] * dnorm(values %*% t(t(beta)) - tau0[, response - 1])
  
  # Calculating the through-tau variables marginal impact on P(y.s >= response) for each row in values
  #   First, calculating deriv_tau.rm1, a row vector of which row i contains the marginal variations in tau_{response - 1} for row i of values
  if (response == 2) deriv.tau.rm1 <- seq(from = gamma[var.name.col, 1], to = gamma[var.name.col, 1], length.out = kNPoints)
  if (response > 2) {
  	deriv.tau.rm1 <- 0
  	for (i in (response - 1):(2)) {
  	  deriv.tau.rm1 <- deriv.tau.rm1 + gamma[var.name.col, i] * exp(values %*% gamma[, i, drop = FALSE])
  	  }
  	  deriv.tau.rm1 <- deriv.tau.rm1 + gamma[var.name.col, 1]
  }
  #   Then, obtain the row of throug-tau contributions to the marginal impact
  tau.marginal <- - deriv.tau.rm1 * dnorm(values %*% beta - tau0[, response - 1])
  
  # Eventually, calculating the vector of marginal effects of variable var.name on the probability
  # Checked: so far, the algorithm is numerically validated
  marginal <- latent.marginal + tau.marginal

  return(marginal.effect <- list(marginal = marginal, latent.marginal = latent.marginal, tau.marginal = tau.marginal, 
                                 values.var = values[, var.name.col]))	
}



ChopitMarginalDisrete <- function(values, var.name.col, response, beta, gamma) {
  # Args
  #   values: row-vector x = x.s =x.tau on which the marginal effect should be calculated (only value[var.name.col] is going to be modified by the
  #           algo, to take values 1 and 0)
  #           Remark: it is assumed that the variables entering the y.star.s equation are the same as the variables entering the threshold 
  #           equations, except for the constant.
  #   var.name.col: column number of the variable for which the marginal effect is calculated
  #   response: response as in P(y.s >= response / values[i, ], beta, gamma)
  #   beta: vector beta (beta[1] is expected to be 0, cf. above)
  #   gamma: matrix that cbinds gamma.1, ..., gamma.K-1, with gamma.k a column vector whose number of rows matches the number of columns in x
  #
  # Returns
  #   marginal.effects: list containing
  #     marginal: marginal impacts of variable in column var.name.col in values on 
  #               P(y.s >= response / values[i, s], beta, gamma)
  #     latent.marginal: same thing, but considering only the marginal impact of the variable through the change in y.star.s
  #     tau.marginal: same thing, but considering only the marginal impact of the variable through the change in thresholds
  
  # Retrieving the right constants
  kK <- ncol(gamma) + 1
  
  # Sanity Check
  if (!is.matrix(values)) stop("ChopitMarginalDisrete: values should be a matrix")
  if (nrow(values) != 1) stop("ChopitMarginalDisrete: values should be a row-vector")
  if (response < 2) stop("ChopitMarginalDisrete: response smaller than 2")   
  if (response > kK) stop("ChopitMarginalDisrete: response biggher than kK from chopit.results")
  if (beta[1] != 0) cat("WARNING: ChopitMarginalDisrete: the first term of the vector beta is not a 0. \n")
  if (!is.matrix(gamma)) stop("ChopitMarginalDisrete: gamma is not a matrix")
  if (nrow(gamma) != ncol(values)) stop("ChopitMarginalDisrete: number of rows in gamma does not match the number of columns in values")

  # Values when var.name is equal to 0 (resp. 1)
  values.0 <- values ; values.0[var.name.col] <- 0
  values.1 <- values ; values.1[var.name.col] <- 1

  # Calculating the tau for each row of values  
  tau0.0 <- values.0 %*% gamma[, 1, drop = FALSE]
  if (kK >1) { 
    for (k in 2:(kK - 1)) tau0.0 <- cbind(tau0.0, tau0.0[, ncol(tau0.0), drop = FALSE] + exp(values.0 %*% gamma[, k, drop = FALSE]))
  }
  tau0.1 <- values.1 %*% gamma[, 1, drop = FALSE]
  if (kK >1) { 
    for (k in 2:(kK - 1)) tau0.1 <- cbind(tau0.1, tau0.1[, ncol(tau0.1), drop = FALSE] + exp(values.1 %*% gamma[, k, drop = FALSE]))
  }
  
  # Calculating the through-y.star.s variable's marginal impact on P(y.s >= response)
  latent.marginal <-  pnorm(values.1 %*% t(t(beta)) - tau0.0[, response - 1]) - pnorm(values.0 %*% t(t(beta)) - tau0.0[, response - 1])
  
  # Calculating the through-tau variables marginal impact on P(y.s >= response) for each row in values
  tau.marginal <- pnorm(values.0 %*% t(t(beta)) - tau0.1[, response - 1]) - pnorm(values.0 %*% t(t(beta)) - tau0.0[, response - 1])
    
  # Eventually, calculating the vector of marginal effects of variable var.name on the probability
  marginal <- pnorm(values.1 %*% t(t(beta)) - tau0.1[, response - 1]) - pnorm(values.0 %*% t(t(beta)) - tau0.0[, response - 1])

  return(marginal.effect <- list(marginal = marginal, latent.marginal = latent.marginal, tau.marginal = tau.marginal))
}



ChopitAverageMarginalDiscrete <- function(chopit.results, var.name) {
  # Same description as for ChopitAverageMarginal, but 
  # Has no method for plot()
  
  # Sanity Checks
  if (!class(chopit.results) == "Chopit") cat("ChopitAverageMarginalSpan: argument chopit.results is not of class 'Chopit' \n")
  if (!is.character(var.name)) stop("ChopitMarginal: var.name is not of type character")
  if (!(var.name %in% colnames(chopit.results$original.x$original.x.s))) stop("ChopitMarginal: var.name does not belong to the regressor names")
  if (!all(chopit.results$original.x$original.x.tau == chopit.results$original.x$original.x.s, na.rm = TRUE)) stop("ChopitMarginal: in
     chopit.results, x.tau != x.s")
  
  # Getting back the right variables
  kK <- chopit.results$constants$kK
  
  ChopitAverageMarginalDiscreteCalculus <- function(chopit.results, var.name) {
  
    # Getting back the right variables
    kK <- chopit.results$constants$kK
    beta0 <- chopit.results$coef$beta0
    beta <- c(0, chopit.results$coef$beta0)
    gamma <- chopit.results$coef$gamma
  
    # Locating the column of var.name
    var.name.col <- match(var.name, colnames(chopit.results$original.x$original.x.s))
  
    # Calculating the average values of all the regressors, including the intercept.
    x.average <- colMeans(chopit.results$original.x$original.x.s, na.rm = TRUE)
  
    # generating values
    values <- matrix(x.average, nrow = 1)

    # Making the calculations: each row contain the results for a given value of x[var.name], with each column being the result for a given 
    # response, going from k = 2 to k = kK.
    TempChopitMarginalDisrete <-  function(response, values, var.name.col, beta, gamma) {  # same as ChopitMarginalDiscrete, but takes resp as 
    	                   																   # first argument
    	marginal.results <- ChopitMarginalDisrete(values = values, var.name.col = var.name.col, response = response, beta = beta, gamma = gamma)
    	return(marginal.results)
    }  
    marginal.results <- lapply(X = 2:(kK), FUN = TempChopitMarginalDisrete, values = values, var.name.col = var.name.col, beta = beta, 
                               gamma = gamma)
    
    # Parsing the results in the appropriate vectors
    marginal <- numeric(0) ; latent.marginal <- numeric(0) ; tau.marginal <- numeric(0)
    for (k in (1:(kK-1))) {
    	marginal <- c(marginal, marginal.results[[k]][[1]])
    	latent.marginal <- c(latent.marginal, marginal.results[[k]][[2]])
    	tau.marginal <- c(tau.marginal, marginal.results[[k]][[3]])
    }
    marginal <- matrix(marginal, nrow = 1) ; latent.marginal <- matrix(latent.marginal, nrow = 1) ; tau.marginal <- matrix(tau.marginal, nrow = 1)  
    
    # Calculating the expected value in each row. Checked.
    prop <- prop.table(table(chopit.results$original.y.s))
    prop <- as.matrix(prop, nrow = ncol(marginal))
    average.marginal <- marginal %*% prop[-1]
    average.latent.marginal <- latent.marginal %*% prop[-1]
    average.tau.marginal <- tau.marginal %*% prop[-1]
  
    return(c(average.marginal = average.marginal, average.latent.marginal = average.latent.marginal, 
               average.tau.marginal = average.tau.marginal))
  }
  
  # Calculating the marginal effects
  average.marginal <- ChopitAverageMarginalDiscreteCalculus(chopit.results = chopit.results, var.name  = var.name)
    
  # First, we defined a function used to calculate the variance through the delta method
  MarginalEffectsJacobiansDiscrete <- function(para, chopit.results, var.name) {
	  kK <- chopit.results$constants$kK
      kGammaNrow <- chopit.results$constants$kGammaNrow 
      # Sanity Checks
      if (length(para) != (ncol(chopit.results$original.x$original.x.s) - 1 + (kK - 1) * chopit.results$constants$kGammaNrow)) 
          stop("marginal.effects.jacobians: lentgh of parameter not conform")
      # Injecting para into chopit.results, in order to make the calculations with different values of para
      beta0 <- para[1:(ncol(chopit.results$original.x$original.x.s) - 1)]
      gamma <- matrix(para[(length(beta0) + 1) : (length(beta0) + (kK - 1) * kGammaNrow)], nrow = kGammaNrow)
      chopit.results$coef$beta0 <- beta0
      chopit.results$coef$gamma <- gamma
      marginal.effects <- ChopitAverageMarginalDiscreteCalculus(chopit.results = chopit.results, var.name  = var.name)
      return(marginal.effects)
  }
  
  # marginal.effects.jacobians: 3 times length(beta0,gamma) matrix, for which the k column contains the first derivative of the 3 marginal effects
  #                             with respect to the k element of c(beta0, gamma)
  marginal.effects.jacobians <- jacobian(MarginalEffectsJacobiansDiscrete, c(chopit.results$coef$beta0, chopit.results$coef$gamma), method = "simple",
                                         chopit.results = chopit.results, var.name = var.name)
  
  #   Extracting the variance-covariance matrix associated to the vector c(beta0, gamma)
  kNpara <- length(chopit.results$coef$beta0) + (kK - 1) * nrow(chopit.results$coef$gamma)
  sigma <- chopit.results$var.cov[1 : kNpara, 1 : kNpara]

  #  Applying the delta method on each lign of marginal.effects.jacobians
  variances <- apply(marginal.effects.jacobians, MARGIN = 1, FUN = function(marginal.effects.jacobian, sigma) {
  	       jacob <- matrix(marginal.effects.jacobian, ncol = length(marginal.effects.jacobian))
  	       variance <- jacob %*% sigma %*% t(jacob)
  	       return(variance)
        }, sigma = sigma)
  marginal.effects.variances <- list(marginal.var = variances[1], latent.marginal.var = variances[2],
                                     tau.marginal.var = variances[3])
  
  
  # Calculating the confidence intervals
  marginal.effects.CI <- list(
    marginal.CI = cbind(average.marginal[1] - qnorm(0.975) * sqrt( marginal.effects.variances$marginal.var), 
                        average.marginal[1] + qnorm(0.975) * sqrt( marginal.effects.variances$marginal.var)),
    latent.marginal.CI = cbind(average.marginal[2] - qnorm(0.975) * sqrt( marginal.effects.variances$latent.marginal.var), 
                              average.marginal[2] + qnorm(0.975) * sqrt( marginal.effects.variances$latent.marginal.var)),
    tau.marginal.CI = cbind(average.marginal[3] - qnorm(0.975) * sqrt( marginal.effects.variances$tau.marginal.var), 
                            average.marginal[3] + qnorm(0.975) * sqrt( marginal.effects.variances$tau.marginal.var))
  )
  
  return(list(marginal.effects = average.marginal, marginal.effects.variances = marginal.effects.variances, 
              marginal.effects.CI = marginal.effects.CI))
}









DoNotExecute <- function() {# VALIDATING NUMERICALLY THE PREVIOUS ALGO
ChopitCumulative <- function(value.l, value.ml, chopit.results, response, var.name) {
  # ChopitCumulative: finds proba(y.s >= respone / value.l, value.ml ; chopit.results$coef)
  
  var.name.col <- match(var.name, colnames(chopit.results$original.x$original.x.s))
  values <- c(value.ml[1:(var.name.col - 1)], value.l, value.ml[(var.name.col):length(value.ml)])  # Works only of variable l is not the last 
                                                                                                   # regressor
  
  # Retrieving the right coefficients
  beta0 <- chopit.results$coef$beta0
  beta <- c(0, chopit.results$coef$beta0)
  gamma <- chopit.results$coef$gamma
  kK <- chopit.results$constants$kK

  # Calculating the tau for each row of values
  tau0 <- as.matrix(values %*% gamma[, 1, drop = FALSE])
  if (kK > 1) { 
    for (k in 2:(kK - 1)) tau0 <- cbind(tau0, tau0[, ncol(tau0), drop = FALSE] + exp(values %*% gamma[, k, drop = FALSE]))
  }
  
  # Calculating the above-described proba
  proba <- pnorm(values %*% beta - tau0[, response - 1])
  
  #Return
  return(proba)
}




library(numDeriv)
l <- 3  # Variable number
kNPoints <- 10
response <- 2
chopit.results <- chopit.pain
marginal <- ChopitMarginal(chopit.pain, "age", response, kNPoints)

for (s in 1:kNPoints) {
print(grad(ChopitCumulative, marginal$values[s, l], value.ml = marginal$values[s, -l], chopit.results = chopit.results, 
     response = response, var.name = "age"))
}
marginal$marginal
}