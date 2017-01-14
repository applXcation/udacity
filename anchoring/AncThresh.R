# TO DO: check on simulated data

AncThresholds <- function(ancreg.results) {
  # space: for each row k, the g-th column gives the space between thresholds k+1 and k in group g
  # space.sd: associated sd
  # tstat.list: list of which the k-th element contains a matrix, whose cell (g, h) tests whether for the equality in the diff between tau_{k+1} and tau_{k} in group g, h
  # rand.eff: if TRUE, uses the results obtained with random effects. FALSE otherwise.
  
  if (ancreg.results$constants$H == 1) rand.eff <- FALSE
  else rand.eff <- TRUE
  
  # Useful matrices
  if (rand.eff == FALSE) {
    para <- ancreg.results$structural.para$structural.para.no.quad
    var.cov <- ancreg.results$structural.var.cov$structural.var.cov.no.quad
  }
  if (rand.eff == TRUE) {
    para <- ancreg.results$structural.para$structural.para.quad
    var.cov <- ancreg.results$structural.var.cov$structural.var.cov.quad
  }
  # Useful constants
  kK <- ancreg.results$constants$kK
  kG <- ancreg.results$constants$kG
  kL <- length(para)  # length of the parameter vector

  space <- matrix(0, nrow = kK - 2, ncol = kG)
  space.sd <- matrix(0, nrow = kK - 2, ncol = kG)
  tstat.list <- list()
  for (k in 2 : (kK - 1)) {

  	tstat.k <- matrix(0, nrow = kG, ncol = kG)
  	for (g in 1 : kG) {
  	  pos.k.g <- kL - (kK - 1) * (kG - g + 1) + 1 + (k - 1) - kG * rand.eff # position of tau.k.g among the parameters
  	  pos.km1.g <- kL - (kK - 1) * (kG - g + 1) + 1 + (k - 2) - kG * rand.eff  # position of tau.km1.g among the parameters
      
      # Space between thresholds
      space[k - 1, g] <- para[pos.k.g] - para[pos.km1.g]
      space.sd[k - 1, g] <- sqrt(var.cov[pos.k.g, pos.k.g] + var.cov[pos.km1.g, pos.km1.g] - 2 * var.cov[pos.k.g, pos.km1.g] )
      
      for (h in 1 : kG) if (g != h) {
  	    pos.k.h <- kL - (kK - 1) * (kG - h + 1) + 1 + (k - 1) - kG * rand.eff  # position of tau.k.h among the parameters
  	    pos.km1.h <- kL - (kK - 1) * (kG - h + 1) + 1 + (k - 2) - kG * rand.eff  # position of tau.km1.h among the parameters
		pos.vector <- c(pos.k.g, pos.km1.g, pos.k.h, pos.km1.h)
        
        # Variance-covariance of the parameters
		var.cov.small <- matrix(0, nrow = 4, ncol = 4)
		for (l in 1 : 4) {
		  for (m in 1 : 4) {
		    var.cov.small[l, m] <- var.cov[pos.vector[l], pos.vector[m]]	
		  }
		}
        
        # Calculating the space between thresholds
        
        
        # Calculating the sd for the spaces between thresholds
        
        # Calculating the t-stat
		num <- (para[pos.k.g] - para[pos.km1.g]) - (para[pos.k.h] - para[pos.km1.h])
		den <- sqrt(sum(diag(var.cov.small)) - 2 * var.cov[pos.k.g, pos.km1.g] - 2 * var.cov[pos.k.g, pos.k.h] 
		+ 2 * var.cov[pos.k.g, pos.km1.h] + 2 * var.cov[pos.km1.g, pos.k.h] 
		- 2 * var.cov[pos.km1.g, pos.km1.h] - 2 * var.cov[pos.k.h, pos.km1.h])
		t.stat <- num / den
		tstat.k[g, h] <- t.stat

  	  }
  	}
  tstat.list[[k - 1]] <- tstat.k
  }
  
  # Return
  return(list(space = space, space.sd = space.sd, tstat.list = tstat.list))
}