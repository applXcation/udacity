AncRegLl <- function(y, x.s, x.bar, group.index, indiv.index,
					 beta.s = NULL, beta.bar.star = NULL, sigma.eps = NULL, tau.star = NULL, sigma.u.star = NULL, 
					 Gamma.d.star = NULL, Gamma.bar.star = NULL, sigma.eta = NULL, tau = NULL, sigma.v.star = NULL,
                     v.z, kG, kK, kN, H, n.quad.mat.sum, list.mat.quad.sum, Gsplit, list.mat.quad.prod, test.rf) {
  # Sanity checks
  # if (nrow(y) != nrow(x.s) | nrow(y) != nrow(x.bar) | nrow(y) != length(group.index) | nrow(y) != length(indiv.index)) stop("AncRegLl: number of row incompatibles among the data matrices")
  # if (nrow(beta.s) != ncol(x.s)) stop("AncRegLl: number of rows in beta.s not equal to number of cols in x.s")
  # if (nrow(beta.bar.star) != ncol(x.bar)) stop("AncRegLl: number of rows in beta.bar.star not equal to number of cols in x.bar")
  # if (length(sigma.eps) != kG) stop("AncRegLl: length of sigma.eps not equal to kG")   
  # if (nrow(tau.star) != kK + 1 | ncol(tau.star) != kG) stop("AncRegLl: dimension in tau.star incompatible with kK and/or kG") 
  
  # Calculation of p.yi = P(Y_{i} / x.s, x) for each respondent
  # If H == 1
  if (H == 1) {
  	if (!test.rf) {
      sigma.eps.bar <- SigmaBar(sigma = sigma.eps, group.index = group.index)
      latent <- x.s %*% beta.s / sigma.eps.bar + x.bar %*% beta.bar.star  # part of the proba that comes into both normal distribution function
      p.yit.u <- pnorm(tau.star[cbind(y + 1, group.index)] - latent) - pnorm(tau.star[cbind(y, group.index)] - latent)
    }
    if (test.rf) {
      sigma.eta.bar <- SigmaBar(sigma = sigma.eta, group.index = group.index)
      latent <- Gamma.d.star[group.index] + x.bar %*% Gamma.bar.star
      p.yit.u <- (y == 1) * (pnorm(tau[y + 1] / sigma.eta.bar - latent)) +
      			 (y > 1 & y < kK) * (pnorm(tau[y + 1] / sigma.eta.bar - latent) - pnorm(tau[y] / sigma.eta.bar - latent)) +
      			 (y == kK) * (1 - pnorm(tau[y] / sigma.eta.bar - latent))
    }
    lli <- log(p.yit.u)
    return(lli)
  }
  
  
  # IF H > 1
  if (H > 1) {
  	if (!test.rf) {
  	  sigma.eps.bar <- SigmaBar(sigma = sigma.eps, group.index = group.index)
  	  sigma.u.bar.star <- SigmaBar(sigma = sigma.u.star, group.index = group.index)
      quad.index <- matrix(1 : H, nrow = H) %x% (numeric(length(y)) + 1)  # index for the quadratures
      p.yit.u.all <- numeric()
  	  for (h in 1 : H) {  # for each h, calculates p.yit.u.h ; concatenates it to p.yit.u.all
        latent <- x.s %*% beta.s / sigma.eps.bar + x.bar %*% beta.bar.star + sqrt(2) * sigma.u.bar.star * v.z[h, 1]
        p.yit.u.h <- pnorm(tau.star[cbind(y + 1, group.index)] - latent) - pnorm(tau.star[cbind(y, group.index)] - latent)
        p.yit.u.all <- c(p.yit.u.all, p.yit.u.h)
  	  }
  	}
  	if (test.rf) {  		
  	  sigma.eta.bar <- SigmaBar(sigma = sigma.eta, group.index = group.index)
  	  sigma.v.bar.star <- SigmaBar(sigma = sigma.v.star, group.index = group.index)
      quad.index <- matrix(1 : H, nrow = H) %x% (numeric(length(y)) + 1)  # index for the quadratures
      p.yit.u.all <- numeric()
  	  for (h in 1 : H) {  # for each h, calculates p.yit.u.h ; concatenates it to p.yit.u.all
        latent <- Gamma.d.star[group.index] + x.bar %*% Gamma.bar.star + sqrt(2) * sigma.v.bar.star * v.z[h, 1]
        p.yit.u.h <- (y == 1) * (pnorm(tau[y + 1] / sigma.eta.bar - latent)) +
      			 (y > 1 & y < kK) * (pnorm(tau[y + 1] / sigma.eta.bar - latent) - pnorm(tau[y] / sigma.eta.bar - latent)) +
      			 (y == kK) * (1 - pnorm(tau[y] / sigma.eta.bar - latent))
        p.yit.u.all <- c(p.yit.u.all, p.yit.u.h)
  	  }
  	}

    # For each individual in indiv.index, calculating p.yi
	# First, for each h, product of each respondent's probabilities over time
	log.p.yit.u.all <- cbind(log(p.yit.u.all))
	p.yi.h <- numeric()
	for (h in 1 : H) {
	  base.row <- (h - 1) * length(indiv.index) + 1 # Row at which the computation should start for this sequence
	  for (gs in 1 : Gsplit) {
	  	row.start <- list.mat.quad.prod[[Gsplit + 1]][gs, 1] + base.row - 1
	  	row.end <- list.mat.quad.prod[[Gsplit + 1]][gs, 2] + base.row - 1
	  	# cat("bla", c(row.start, row.end), "\n")
	  	# print(list.mat.quad.prod[[Gsplit + 1]])
	    p.yi.h  <- c(p.yi.h, exp(log.p.yit.u.all[row.start : row.end, ] %*% list.mat.quad.prod[[gs]]))
	  }
	}
	p.yi.h <- p.yi.h[order(rep.int(unique(indiv.index), times = H), (1 : H) %x% (numeric(kN) + 1))]
	p.yi.h <- rbind(p.yi.h)
    
	# p.yi.h <- rbind(tapply(p.yit.u.all, INDEX = interaction(quad.index, indiv.index), FUN = prod))  # the outcome is vector ordered by respondent, and for 

	# Then for each respondent, summing over h
    p.yi <- numeric()
    for (gs in 1 : Gsplit) {
      p.yi <- c(p.yi, p.yi.h[, list.mat.quad.sum[[Gsplit+1]][gs, 1] : list.mat.quad.sum[[Gsplit+1]][gs, 2], drop = FALSE] %*% list.mat.quad.sum[[gs]])
    }
    return(log(p.yi))
    
    	# Older code (without Gsplit)
    # new.indiv.index <- (1 : kN) %x% rep.int(1, times = H)  # New index for the respondents
    # p.yi <- tapply(X = p.yi.h, INDEX = new.indiv.index, FUN = function(p.yi.all) {
      # return(sum(p.yi.all * v.z[, 2]))  # For each individual, each p(y.i, h) is multiplied by the h-th weight. Then the terms are summed
    # })
    # p.yi <- 1 / sqrt(pi) * p.yi
    # return(log(p.yi))
  # }
  }
 
}



AncRegLlCompound <- function(coef, y, x.s, x.bar, group.index, indiv.index, v.z, kG, kK, kN, H, n.quad.mat.sum, list.mat.quad.sum, Gsplit, 
							 list.mat.quad.prod, test.rf) {
  # AncRegLlCompound: Transforms 'coef' into the coefficient matrices, and put them into AncRegLl to get the log-likelihood
  # Args
  #   coef: composed of
  #     tau.star.0: same as tau.star, without first and last rows (i.e. without the infinites)
  
  # count stands for the latest element of coef used
  if (!test.rf) {
  beta.s <- coef[1 : ncol(x.s)] ; beta.s <- cbind(beta.s) ; count <- ncol(x.s)
  beta.bar.star <- coef[(count + 1) : (count + ncol(x.bar))] ; beta.bar.star <- cbind(beta.bar.star) ; count <- count + ncol(x.bar)
  kappa.eps.0 <- coef[(count + 1) : (count + kG - 1)] ; kappa.eps.0 <- cbind(kappa.eps.0) ; count <- count + kG - 1
  sigma.eps.0 <- exp(kappa.eps.0)
  sigma.eps <- cbind(c(1, sigma.eps.0))
  gamma.star.0 <- coef[(count + 1) : (count + (kK - 1) * kG)] ; gamma.star.0 <- matrix(gamma.star.0, nrow = kK - 1, ncol = kG) ; count <- count + kG + kK - 1
  # Calculating tau.star from gamma.star.0
  tau.star.0 <- rbind(gamma.star.0[1, ])
  for (k in 2:(kK - 1)) {
    tau.star.0 <- rbind(tau.star.0,
    					tau.star.0[k - 1, ]+ exp(gamma.star.0[k, ]))
  }
  tau.star <- rbind(-Inf, tau.star.0, +Inf)
  if (H == 1) sigma.u.star <- numeric(kG)
  if (H > 1) {
  	kappa.u.star <- coef[(length(coef) - kG + 1) : length(coef)]
  	sigma.u.star <- exp(kappa.u.star)
  }
  sigma.u.star <- cbind(sigma.u.star)
  
  lli <- AncRegLl(y = y, x.s = x.s, x.bar = x.bar, group.index = group.index, indiv.index = indiv.index, beta.s = beta.s, 
  				  beta.bar.star =  beta.bar.star, sigma.eps = sigma.eps, tau.star = tau.star, sigma.u.star = sigma.u.star, v.z = v.z, kG = kG, kK = kK, 
  				  kN = kN, H = H, n.quad.mat.sum = n.quad.mat.sum, list.mat.quad.sum = list.mat.quad.sum, Gsplit = Gsplit, 
  				  list.mat.quad.prod = list.mat.quad.prod, test.rf = test.rf)
  }

  if (test.rf) {
  	# count stands for the latest element of coef used
    Gamma.d.star.0 <- coef[1 : (kG - 1)] ; count <- kG - 1
    Gamma.bar.star <- coef[(count + 1) : (count + ncol(x.bar))] ; count <- count + ncol(x.bar)
    kappa.eta.0 <- coef[(count + 1) : (count + kG - 1)] ; count <- count + kG - 1
    gamma.0 <- coef[(count + 1) : (count + kK - 1)] ; count <- count + kK - 1
   
    if (H > 1) {
      kappa.v.star <- coef[(count + 1) : (count + kG)] ; count <- count + kG
      sigma.v.star <- exp(kappa.v.star)
    }
    else sigma.v.star <- numeric(kG)
   
    tau.0 <- c(gamma.0[1])
    for (k in 2:(kK - 1)) {
      tau.0 <- c(tau.0, tau.0[k - 1] + exp(gamma.0[k]))
    }
    
    Gamma.d.star <- c(0, Gamma.d.star.0)
    sigma.eta.0 <- exp(kappa.eta.0)
    sigma.eta  <- c(1, sigma.eta.0)
    
    tau <- c(-Inf, tau.0, +Inf)    
    
    
    lli <- AncRegLl(y = y, x.s = NULL, x.bar = x.bar, group.index = group.index, indiv.index = indiv.index,
     			    Gamma.d.star = Gamma.d.star, Gamma.bar.star = Gamma.bar.star, sigma.eta = sigma.eta, tau = tau, sigma.v.star = sigma.v.star,
  				    v.z = v.z, kG = kG, kK = kK, kN = kN, H = H, n.quad.mat.sum = n.quad.mat.sum, list.mat.quad.sum = list.mat.quad.sum, 
  				    Gsplit = Gsplit, list.mat.quad.prod = list.mat.quad.prod, test.rf = test.rf)

  }
    
  # Return
  return(lli)
}



SigmaBar <- function(sigma, group.index) {
  # SigmaUBar: calculates a col-vector for wich each observation corresponds to an individual observation (i,t) of sigma.g, depending on the group
  # 		   the respondent belongs to.
  # Remark: group.index has to be sorted (and obviously, the ordering of the rows in x has to be conistent with this index). See XBar for
  #   		more details.
  
  # Sanity checks (desactivited because the two conditions are guaranteed in AncReg AND because it slows the algo down a lot)
#  if (any(group.index != group.index[order(group.index)])) stop("SigmaUBar: group.index is not sorted (meaning that x and index are not sorted by the levels of group.index either")
#  if (nrow(sigma) != length(levels(as.factor(group.index)))) stop("SigmaUBar: the number of levels in group.index is not equal to the length of sigma")

  # Constants
  nt.g <- as.vector(table(group.index))  # Matrix that contains the number of observations in each group
  kG <- length(nt.g)  # number of groups
  
  # Calculation
  sigma.bar <- numeric(0)
  for (g in 1 : kG) {
    sigma.bar <- c(sigma.bar, rep.int(sigma[g], nt.g[g]))
  }
  return(sigma.bar)
}



StructuralPara <- function(coef, kG, kK, kS, l.g, quad, test.rf = FALSE) {
  # count stands for the latest element of coef used
  # Getting back the reduced-form 'star' parameters
  
  if (!test.rf) {
  beta.s <- coef[1 : kS] ; beta.s <- cbind(beta.s) ; count <- kS
  beta.bar.star <- coef[(count + 1) : (count + sum(l.g))] ; beta.bar.star <- cbind(beta.bar.star) ; count <- count + sum(l.g)
  kappa.eps.0 <- coef[(count + 1) : (count + kG - 1)] ; kappa.eps.0 <- cbind(kappa.eps.0) ; count <- count + kG - 1
  sigma.eps.0 <- exp(kappa.eps.0)
  gamma.star.0 <- coef[(count + 1) : (count + (kK - 1) * kG)] ; gamma.star.0 <- matrix(gamma.star.0, nrow = kK - 1, ncol = kG) ; count <- count + kG * (kK - 1)
  if (quad == TRUE) {
  	kappa.u.star <- coef[(count + 1) : (count + kG)] ; count <- count + kG
  	sigma.u.star <- exp(kappa.u.star)
  }
  # Calculating tau.star from gamma.star.0
  tau.star.0 <- rbind(gamma.star.0[1, ])
  for (k in 2:(kK - 1)) {
    tau.star.0 <- rbind(tau.star.0,
    					tau.star.0[k - 1, ] + exp(gamma.star.0[k, ]))
  }
  
  # Getting back beta.bar
  extended.sigma.eps.0 <- cbind(c(rep.int(1, l.g[1]) ,mapply(rep.int, sigma.eps.0, l.g[-1])))
  beta.bar <- beta.bar.star * extended.sigma.eps.0
  
  # Getting back tau
  extended.sigma.eps.0 <- cbind(rep.int(1, kK - 1) ,mapply(rep.int, sigma.eps.0, kK - 1))  
  tau.0 <- tau.star.0 * extended.sigma.eps.0
  
  # Getting back sigma.u
  if (quad == TRUE) sigma.u <- sigma.u.star * c(1, sigma.eps.0)
  
  # Structural para: dividing everything by beta^1_s
  # structural.para <- c(beta.s, beta.bar, sigma.eps.0, tau.0)  / beta.s[1]
  if (quad == FALSE) structural.para <- c(beta.s, beta.bar, sigma.eps.0, tau.0)
  if (quad == TRUE) structural.para <- c(beta.s, beta.bar, sigma.eps.0, tau.0, sigma.u)
  }
  
  if (test.rf) {
  	# Getting back the reduced-form rf coefs
    Gamma.d.star.0 <- coef[1 : (kG - 1)] ; count <- kG - 1
    Gamma.bar.star <- coef[(count + 1) : (count + sum(l.g))] ; count <- count + sum(l.g)
    kappa.eta.0 <- coef[(count + 1) : (count + kG - 1)] ; count <- count + kG - 1
    gamma.0 <- coef[(count + 1) : (count + kK - 1)] ; count <- count + kK - 1
    kappa.v.star <- coef[(count + 1) : (count + kG)] ; count <- count + kG
    
    # Getting back the structural rf coefs
    tau.0 <- rbind(gamma.0[1])
    for (k in 2:(kK - 1)) {
      tau.0 <- rbind(tau.0,
       				  	  tau.0[k - 1] + exp(gamma.0[k]))
    }
    sigma.eta.0 <- exp(kappa.eta.0) ; sigma.eta <- c(1, sigma.eta.0)
    sigma.v <- exp(kappa.v.star) * sigma.eta
    Gamma.d.0 <- Gamma.d.star.0 * sigma.eta.0
    extended.sigma.eta <- cbind(c(rep.int(1, l.g[1]) , mapply(rep.int, sigma.eta.0, l.g[-1])))
    Gamma.bar <- Gamma.bar.star * extended.sigma.eta

    structural.para <- c(Gamma.d.0, Gamma.bar, sigma.eta.0, tau.0, sigma.v)
  }
  
  return(structural.para)
}

