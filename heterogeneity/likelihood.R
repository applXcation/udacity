# Calculating the LL
# Plan:
#   First, we compute the tau
#   Then we compute the log-likelihood for self-evaluations
#   Then we compte the  log-likelihood for the evaluations of each vignette
#   Eventually, we sum the whole


ChopitTau <- function(x.tau, gamma, kK) {
  # computes the tau
  #
  # Args:
  #   x.tau: matrix with kN rows of variables explaining the thresholds
  #   gamma: matrix that cbinds gamma.1, ..., gamma.K-1, with gamma.k a column vector whose number of rows matches the number of columns in x.tau
  #   kK: number of statements categories
  #
  # Returns:
  #   tau: matrix that cbinds tau.0, tau.1, ..., tau.K, with tau.k a column vector of size kN that contains the value for tau.k for every 
  #   individuals
  #
  # Checked on simulated data: recover the right taus given x.tau, gamma and kK
  
  # Sanity Checks
  if (ncol(x.tau) != nrow(gamma)) stop("ChopitTau: ncol(x.tau) != nrow(gamma)")
  
  # Constants
  kN <- nrow(x.tau)
  
  # Calculations
  tau <- matrix(-Inf, nrow = kN)  # Adding tau.0
  tau <- cbind(tau, x.tau %*% gamma[, 1, drop = FALSE])  # Adding tau.1
  if (kK >1) { 
    for (k in 2:(kK - 1)) tau <- cbind(tau, tau[, ncol(tau), drop = FALSE] + exp(x.tau %*% gamma[, k, drop = FALSE]))  # if K > 1, adding tau.k
  }
  tau <- cbind(tau, matrix(+Inf, nrow = kN))
  
  # Return
  return(tau)
}


ChopitLlSelfEval <- function(y.s, x.s, tau, beta0, x.sigma = NULL, kappa0 = NULL, heterosk = FALSE) {
  # Args:
  #   y.s: reported self-evaluation
  #   x.s: x for the self-evaluation equation (including the first row made of 1)
  #   tau: matrix that cbinds tau.0, tau.1, ..., tau.K, with tau.k a column vector of size kN that contains the value for tau.k for every 
  #   individuals
  #   beta0: corresponding coefficients except the first line (i.e. without the intercept)
  #
  # Return:
  #   chopit.ll.self.eval, the vector of all the individuals' self-evaluation contributions to the log-likelihood
  
  # Defining a function that calculates, for every individual, the contribution of his/her self-eval to the log-likelihood
  kK <- ncol(tau) - 1
  
  if (!heterosk) {
    beta <- rbind(0, beta0)  # proper value of beta
    x.beta <- x.s %*% beta
    chopit.l.self.eval <- (y.s == 1) * pnorm(tau[cbind(1:nrow(tau), y.s + 1)] - x.beta) +
    					  (y.s > 1 & y.s < kK) * (pnorm(tau[cbind(1:nrow(tau), y.s + 1)] - x.beta) - pnorm(tau[cbind(1:nrow(tau), y.s)] - x.beta)) +
    					  (y.s == kK) * (1 - pnorm(tau[cbind(1:nrow(tau), y.s)]  - x.beta))
    chopit.ll.self.eval <- log(chopit.l.self.eval)
    chopit.ll.self.eval[chopit.ll.self.eval == -Inf] <- log(5e-324) # In the case the diff between proba is too small, log() = -Inf. In 
  	# this case, the individual likelihood is coerced to take the closest to 0 value R knows
   
    # Return
    return(chopit.ll.self.eval)
  }

  if (heterosk) {
    beta <- rbind(0, beta0)  # proper value of beta
    kappa <- rbind(0, kappa0)  # proper value of kappa
    x.beta <- x.s %*% beta
    x.kappa <- x.sigma %*% kappa
    chopit.l.self.eval <- (y.s == 1) * pnorm(1 / exp(x.kappa) * (tau[cbind(1:nrow(tau), y.s + 1)] - x.beta)) +
    					  (y.s > 1 & y.s < kK) * (pnorm(1 / exp(x.kappa) * (tau[cbind(1:nrow(tau), y.s + 1)] - x.beta))
    					                          - pnorm(1 / exp(x.kappa) * (tau[cbind(1:nrow(tau), y.s)] - x.beta))) +
    					  (y.s == kK) * (1 - pnorm(1 / exp(x.kappa) * (tau[cbind(1:nrow(tau), y.s)]  - x.beta)))
    chopit.ll.self.eval <- log(chopit.l.self.eval)
    chopit.ll.self.eval[chopit.ll.self.eval == -Inf] <- log(5e-324) # In the case the diff between proba is too small, log() = -Inf. In 
  	# this case, the individual likelihood is coerced to take the closest to 0 value R knows
   
    # Return
    return(chopit.ll.self.eval)
  }
}


ChopitLlVignEvalBak <- function(y.v, tau, theta, sigma.tilde.v) {
  # NOT INTENDED FOR USE. Kept for keeping an example of bad practice in term of use of loops.
  # Args (all args should be vectors or row-matrix)
  #   y.v: one reported vignette evaluation
  #   tau: matrix that cbinds tau.0, tau.1, ..., tau.K, with tau.k a column vector of size kN that contains the value for tau.k for every 
  #   individuals
  #   theta: vignette-specific intercepts
  #   sigma.tilde.v: ln(1/sigma.v) with sigma.v the standard-deviation
  #
  # Return:
  #   chopit.ll.vignette.eval, the vector of all the individuals' vignette evaluation contributions to the log-likelihood (for one vignette)

  # Sanity Check
  if (length(theta) > 1) stop("ChopitLlVignEval: theta contains more than one row")

  # Defining a function that calculates, for every individual, this contribution
  ChopitLlVignEvali <- function(y.vi, taui, theta, sigma.tilde.v) {
  	# For a given statement y.si, calculation of log-likelihood
  	chopit.ll.vign.evali <- log(pnorm(exp(sigma.tilde.v) * (taui[y.vi + 1] - theta)) - pnorm(exp(sigma.tilde.v) * (taui[y.vi] - theta)))
  	return(chopit.ll.vign.evali)
  }
  
  # Evaluation the partial log-likelihood
  chopit.ll.vign.eval <- mapply(FUN = ChopitLlVignEvali, y.v, data.frame(t(tau)), theta, sigma.tilde.v)
 
  # Return
  return(chopit.ll.vign.eval)
}



ChopitLlVignEval <- function(y.v, tau, theta, sigma.tilde.v) {
  # Args (all args should be vectors or row-matrix)
  #   y.v: one reported vignette evaluation
  #   tau: matrix that cbinds tau.0, tau.1, ..., tau.K, with tau.k a column vector of size kN that contains the value for tau.k for every 
  #   individuals
  #   theta: vignette-specific intercepts
  #   sigma.tilde.v: ln(1/sigma.v) with sigma.v the standard-deviation
  #
  # Return:
  #   chopit.ll.vignette.eval, the vector of all the individuals' vignette evaluation contributions to the log-likelihood (for one vignette)

  # Sanity Check
  if (length(theta) > 1) stop("ChopitLlVignEval: theta contains more than one row")
 
  # Evaluation the partial log-likelihood
  kK <- ncol(tau) - 1
  chopit.l.vign.eval <- (y.v == 1) * pnorm(exp(sigma.tilde.v) * (tau[cbind(1:nrow(tau), y.v + 1)] - theta)) +
    					  (y.v > 1 & y.v < kK) * (pnorm(exp(sigma.tilde.v) * (tau[cbind(1:nrow(tau), y.v + 1)] - theta))
    					                          - pnorm(exp(sigma.tilde.v) * (tau[cbind(1:nrow(tau), y.v)] - theta))) +
    					  (y.v == kK) * (1 - pnorm(exp(sigma.tilde.v) * (tau[cbind(1:nrow(tau), y.v)]  - theta)))
  chopit.ll.vign.eval <- log(chopit.l.vign.eval)
  # chopit.ll.vign.eval <- log(pnorm(exp(sigma.tilde.v) * (tau[cbind(1:nrow(tau), y.v + 1)] - theta))
                             # - pnorm(exp(sigma.tilde.v) * (tau[cbind(1:nrow(tau), y.v)] - theta)))
  chopit.ll.vign.eval[chopit.ll.vign.eval == -Inf] <- log(5e-324) # In the case the diff between proba is too small, log() = -Inf. In 
  	# this case, the individual likelihood is coerced to take the closest to 0 value R knows
 
  # Return
  return(chopit.ll.vign.eval)
}




ChopitLl <- function(y.s, y.v, x.s, x.tau, beta0, gamma, theta, sigma.tilde.v, kK, chopit.envir,
                     x.sigma = NULL, kappa0 = NULL, heterosk = FALSE) {
  # Args:
  #   y.s: reported self-evaluation
  #   y.v: cbind() on all the vignettes (each column is the observations for one vignette)
  #   x.s: x for the self-evaluation equation (including the first row made of 1)
  #   tau: matrix that cbinds tau.0, tau.1, ..., tau.K, with tau.k a column vector of size kN that contains the value for tau.k for every 
  #   individuals
  #   beta0: corresponding coefficients except the first line (i.e. without the intercept)
  #   gamma: matrix that cbinds gamma.1, ..., gamma.K-1, with gamma.k a column vector whose number of rows matches the number of columns in x.tau
  #   theta: row vector in which row j is theta_j
  #   sigma.tilde.v: row vector in which row j is log(1/sigma.vj) for vignette j
  #   kK: number of statements categories
  #   chopit.envir: environment created by the call to the Chopit function
  # NB: if adding one parameter, should be also taken into account in ChopitLlCompound()
  #
  # Returns:
  #   A vector of all the contributions (self-assessment, vignette assessments) to the log-likelihood
  #   Assigns the number of contributions to the likelihood in each sub-part of the likelihood to an object in the main Chopit frame 
  #   (contrib.number)
  
  tau <- ChopitTau(x.tau= x.tau, gamma = gamma, kK = kK) 
  if (!heterosk) chopit.ll.self.eval <- ChopitLlSelfEval(y.s = y.s, x.s = x.s, tau = tau, beta0 = beta0)  # Atomic vector containing the individual contribution to LL in self-evaluation
  if (heterosk) chopit.ll.self.eval <- ChopitLlSelfEval(y.s = y.s, x.s = x.s, tau = tau, beta0 = beta0, x.sigma = x.sigma,
                                                        kappa0 = kappa0, heterosk = heterosk)
  chopit.ll.vign.eval <- mapply(FUN = ChopitLlVignEval, data.frame(y.v), theta, sigma.tilde.v, MoreArgs = list(tau = tau), SIMPLIFY = FALSE) 
                         # In the case there is only one vignette variable, chopit.ll.vign.eval is a list with a single item, containing the 
                         # contribution to the LL for this vignettes. In the case there are several vignettes, each element of the list contains 
                         # the contributions to the LL for one vignette.
                           
  chopit.ll.self.eval <- na.exclude(chopit.ll.self.eval) ; chopit.ll.vign.eval <- lapply(chopit.ll.vign.eval, na.exclude)  # Getting rid of NA
  																											        # as explained in Chopit.R
  contrib.number <- c(length(chopit.ll.self.eval), lapply(chopit.ll.vign.eval, length))
  assign("contrib.number", contrib.number, envir = chopit.envir)
  ll.contrib <- c(chopit.ll.self.eval, unlist(chopit.ll.vign.eval))  # Log-likelihood
  
  # Return
  return(ll.contrib)
}




ChopitLlCompound <- function(par, y.s, y.v, x.s, x.tau, kK, kBeta0Nrow, kGammaNrow, kV, chopit.envir,
                             heterosk, x.sigma, kKappa0Nrow) {
  # Function that makes ChopitLl usable by optim(). It takes the parameters in compound form, and then parses them to ChopitLl()
  #
  # Args:
  #   par = (beta0', bar.gamma, theta, sigma.tilde.v) with bar.gamma being equal to (gamma^1', ..., gamma^K')
  #   kK: number of statement categories
  #   kBeta0Nrow0: Number of parameters in beta0
  #   kGammaNrow: Number of gamma parameters for each threshold equation
  #   kV: Number of vignettes
  #   chopit.envir: environment created by the call to the Chopit function
  #   Other arguments as described in ChopitLl
  #
  # Returns:
  #   Minus the log-likelihood
    
  # Sanity Checks:
  #   par should be a simple vector (no dimension)
  if (!is.vector(par)) stop("ChopitLlCompound: par is not a simple vector")
  
  # Parsing the arguments
  # A special attention has been dedicated to the dimensions
  beta0 <- t(t(par[1 : kBeta0Nrow]))
  gamma <- matrix(par[(length(beta0) + 1) : (length(beta0) + (kK - 1) * kGammaNrow)], nrow = kGammaNrow)
  theta <- t(par[(length(beta0) + length(gamma) + 1) : (length(beta0) + length(gamma) + kV)])
  sigma.tilde.v <- t(par[(length(beta0) + length(gamma) + length(theta) + 1) : (length(beta0) + length(gamma) + length(theta) + kV)])
  if (heterosk)  kappa0 <- t(t(par[(length(beta0) + length(gamma) + length(theta) + length(sigma.tilde.v) + 1) : (length(beta0) + length(gamma) + length(theta) + length(sigma.tilde.v) + kKappa0Nrow)]))
  
  # Running ChopitLl
  if (!heterosk) ll.contrib <- ChopitLl(y.s = y.s, y.v = y.v, x.s = x.s, x.tau = x.tau, beta0 = beta0, gamma = gamma, theta = theta, 
                         sigma.tilde.v = sigma.tilde.v, kK = kK, chopit.envir = chopit.envir)
  if (heterosk) ll.contrib <- ChopitLl(y.s = y.s, y.v = y.v, x.s = x.s, x.tau = x.tau, beta0 = beta0, gamma = gamma, theta = theta, 
                         sigma.tilde.v = sigma.tilde.v, kK = kK, chopit.envir = chopit.envir,
                         x.sigma = x.sigma, kappa0 = kappa0, heterosk = heterosk)
  
  # Return
  return(ll.contrib)
}


