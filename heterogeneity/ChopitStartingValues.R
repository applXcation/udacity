ChopitStartingValues <- function(y.s, x.s, kK, kBeta0Nrow, kGammaNrow, kV) {
  # ChopitStartingValues()
  # Description:
  #   Functions that looks for good starting values, given the data. First, it calculates, for each equation (self-assessment, each vignette 
  #   assessment) the value of the coefficients when then Xs have no impact on the thresholds.
  #   The beta thesholds of the self-assessment equation are taken as starting values.
  #   Then all the thetas are set as the average value of hat(y.star.s), and the sigmas as the standard error of hat(y.star.s) multiplied by 1.5. 
  #   The fact I multiply this s.d. by  1.5 helps to make sure individual likelihood are not getting too high or too low, but are staying close to
  #   0.5 (the higher sigma, the smaller 1/sigma*(threshold[j] - theta), the closer pnorm(1/sigma*(threshold[j] - theta)) gets to 0)
  # Args: 
  #   cf. description of arguments in ChopitLlCompound
  # Value:
  #   a vector of proper initial parameters

   # Initial beta parameters
  if (kBeta0Nrow > 0) polr.y.s <- polr(as.factor(y.s) ~ x.s[, -1], method = "probit")
  else polr.y.s <- polr(as.factor(y.s) ~ 1, method = "probit")
  beta0.init <- polr.y.s$coefficients

  # Initial theta and sigmas
  lp.y.s <- polr.y.s$lp  # linear prediction of y.star.s, say hat(y.star.s)
  theta.init <- numeric(kV) + median(lp.y.s)
  sigma.v <- 1.5*sd(lp.y.s) ;   sigma.tilde.v.init <- numeric(kV) + log(1/sigma.v)
  
  # Initial gammas
  # First, calculcating the first column, i.e. the constant of each threshold equation. Work also when only one tau (kK = 2).
  tau.y.s <- polr.y.s$zeta
  gamma.col.1 <- c(tau.y.s[1], log(diff(tau.y.s)))  
  # Then we add columns of 0 if there are more than one gammas per equation
  if (kGammaNrow > 1) {
  	gamma.col.m1 <- matrix(0 ,nrow = kGammaNrow - 1, ncol = kK - 1)
  	gamma.init <- rbind(gamma.col.1, gamma.col.m1)
  }
  else gamma.init <- gamma.col.1
  
  # Returns
  par.init <- c(beta0.init, gamma.init, theta.init, sigma.tilde.v.init)
  return(par.init)
}