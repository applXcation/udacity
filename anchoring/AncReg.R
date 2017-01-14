# LOADING SUB-ROUTINES
if (Sys.info()["sysname"] == "Darwin" & Sys.info()["user"] == "xavier") {
  ancreg.routines.directory = "~/Documents/Travail/Boulot Fac/Doctorat/1_ Productions personnelles/28_Anchoring Regressors/3_programmes/"
}

if (Sys.info()["sysname"] == "Linux" & Sys.info()["user"] == "x.fontaine") {
  ancreg.routines.directory = "~/U/Travail/AncReg/0_Routines/"
}

source(paste(ancreg.routines.directory,"AncRegLL.R", sep = ""), chdir = TRUE)

AncReg <- function(formula, data, indiv.index, group.index, H = 1, Gsplit = 160, multico.check = TRUE, init.para, test.rf = FALSE) {
  # Args
  #   Formula: list that contains
  #     y: formula of type 'y ~ x', with y the response and x the variables whose coef are allowed to change across individuals
  #     x.s: formula of type ' ~ x.s', with x.s the variables whose impact is assumed to be the same for all individuals
  #   indiv.index: index of repondents. A same individual may belong to several groups at the same time. In which case, the algo
  #                will "split" the respondents's observations and consider him as several individuals (one per group).
  #   group.index: index of the groups
  #   H: number of points for the quadrature ; 1 by default
  #   Gsplit: number of groups into which the observations are to be splitted to make the algorithm faster with panel data.
  #   multico.check: should a multicolinearity check be performed ? Remark: uses polr(), which sometimes won't find proper starting values, which 
  #     		     leads the program to stop its executions
  #   init.para: if specified and H == 1, used to estimate the model without RE. If specified and H > 1, this first round is skipped, and used to
  #     		   estimate directly the model with RE.
  #   test.rf: logical. TRUE if the alternative model (reduced form) with same thresholds for all groups, but group dummies and differenciated impact
  #   	       of parent death should be estimated instead of the main model. The results are displayed in a different form as compared to
  #            what they are when test.rf is not specified. Besides
  #
  # Return
  
  # LIBRARIES
  library(maxLik)
  library(MASS)
  library(numDeriv)
  library(gaussquad)
  
  # SANITY CHECKS
  if (nrow(data) != length(indiv.index) | nrow(data) != length(group.index)) stop("AncReg: Number of rows in data not compatible with length of indiv.index and / or group.index")
  if (!is.numeric(group.index)) stop("AncReg: group.index is not numeric")
  if (any(as.numeric(names(table(group.index))) != (1 : max(group.index, na.rm = TRUE)))) stop("AncReg: group.index is not ordered from 1 to its maximum (should be a sequence of step 1)")
  if (test.rf) cat("Important : since test.rf is TRUE, the variables in x.s are expected to be included in x, at the BEGINNING of the formula. \n")
  if (test.rf & H == 1) stop("Ancreg: if test.rf is true, H is expected to be higher than 1 (results for reduced-form are displayed only with random effects")
 
 
  # INDEXES
  # creating a new individual index vector. Individuals belonging to two groups (depending on the period) will be "splitted" into two.
  # The unused levels are dropped
  old.index <- indiv.index
  indiv.index <- interaction(group.index, indiv.index, drop = TRUE, lex.order = TRUE)
  
  # ordering the observations by group.index, indiv.index
  ordering <- order(group.index, indiv.index)
  data <- data[ordering, ]
  indiv.index <- indiv.index[ordering]
  group.index <- group.index[ordering]
  old.index <- old.index[ordering]
  
  
  # DATA MATRICES (including NA dropping)
  # Parsing the formula
  f.y <- formula[["y"]]
  f.x.s <- formula[["x.s"]]
  
  # Getting back y 
  mf.y <- model.frame(formula = f.y, data = data, na.action = NULL, drop.unused.levels = TRUE)
  y <- model.response(data = mf.y) ; y <- cbind(y)
  # Sanity checks on y
  if (any(as.numeric(names(table(y))) != (1 : max(y, na.rm = TRUE)))) stop("AncReg: y not ordered from 1 to its maximum (should be a sequence of step 1)")
  if (any(table(y, group.index) == 0)) stop("AncReg: all the groups are not using the same values of y")
  
  # Getting back x.s
  mf.x.s <- model.frame(formula = f.x.s, data = data, na.action = NULL, drop.unused.levels = TRUE)
  x.s <- model.matrix(object = f.x.s, data = mf.x.s) ; x.s <- cbind(x.s)  
  x.s <- x.s[, -1, drop = FALSE]  # not forgetting to drop the constant
  
  # Dropping NA in y, x.s, index, group.index and in data (in order to generate x.list hereafter) when there are NA in y, x.s, x, index, group.index 
  mf.x <- model.frame(formula = f.y, data = data, na.action = NULL, drop.unused.levels = TRUE)
  x <- model.matrix(object = f.y, data = mf.x) ; x <- cbind(x)
  complete <- complete.cases(y, x.s, x, indiv.index, group.index)
  y <- y[complete, , drop = FALSE]
  x.s <- x.s[complete, , drop = FALSE]
  data <- data[complete, ]
  indiv.index <- indiv.index[complete] ; indiv.index <- droplevels(indiv.index)
  group.index <- group.index[complete]
  old.index <- old.index[complete]
  cat("Number of respondents:", nlevels(as.factor(old.index)), "\n")
  cat("Number of observations:", nrow(y), "\n")  
  # for (i in 1:ncol(x.s)) print(table(x.s[, i], group.index))
  for (i in 1:ncol(x.s)) print(table(x.s[, i]))
  rm(x)  
  
  # For each group, creating a data matrix x. Putting all the G matrices in a list x.kist.
  #   The levels not used in a given group are dropped from the dataframe for this group, meaning that no dummy is going to be created 
  #   for these levels.
  #   This should be enough to take care of most of the multi-co problems
  x.list <- list()
  for (g in 1:max(group.index)) {
  	mf.x.g <- model.frame(formula = f.y, data = data[group.index == g, ], na.action = NULL, drop.unused.levels = TRUE)
  	x.g <- model.matrix(object = f.y, data = mf.x.g) ; x.g <- cbind(x.g)
  	x.list[[g]] <- x.g[, -1, drop = FALSE]  # not forgetting to drop the constant
  	# Checking there is no constant column in the x.g matrices. Otherwise, no convergence is possible, due to multi-co between thresh and this col
  	col.var.g <- apply(x.list[[g]], MARGIN = 2, FUN = var, na.rm = TRUE)
  	if (any(col.var.g == 0)) {
  	  cat("AncReg: the variable(s) ", colnames(x.list[[g]])[col.var.g == 0], " are constants for group ", g, ". Dropping them for this group. \n", 
  	      sep = "")
  	  x.list[[g]] <- x.list[[g]][, col.var.g != 0, drop = FALSE]
  	}
  }

  # DEALING WITH PERFECT MULTI-CO
  # For each group, dropping variables for which we have multi-colinearity. It should be noted that x.s is not included here, so if
  # its inclusion causes multicolinearity (in ill cases)  
  if (multico.check) {
    kG <- length(table(group.index))
    for (g in 1 : kG) {
      if (ncol(x.list[[g]]) > 2) {  # Dealing with multico matters only if there is more than the constant in x.s
    	x.temp <- x.list[[g]]
        temp.polr <- polr(as.factor(y)[group.index == g] ~ x.temp , method = "probit")
        names <- gsub("x.temp", "", names(temp.polr$coefficients))
        col.to.drop.x <- !(colnames(x.temp) %in% names)  # Index of the columns in x.s to be dropped
        if (any(col.to.drop.x)) {
          cat("Ancreg: x is not full-rank in group", g, ". Dropping the following columns:", colnames(x.temp)[col.to.drop.x], "\n")
          x.list[[g]] <- x.list[[g]][, !col.to.drop.x, drop = FALSE]
        }
        rm(temp.polr, x.temp)
      }
    }
  }
  
  
  # CONSTANTS
  kG <- length(table(group.index))
  kS <- ncol(x.s)
  kK <- max(y)
  kN <- length(table(indiv.index))
  re.only <- (H > 1 & !missing(init.para))


  # L.G, X.BAR
  # Creating l.g, matrix that contains the number of variables in x for each group, L_g. If test.rf, then x contains the variables in x.s as well.
  l.g <- sapply(x.list, ncol)

  # Creating x.bar
  for (g in 1 : kG) {
    if (g == 1) {
      x.g <- x.list[[g]]
      x.bar <- cbind(x.g, matrix(0, nrow = nrow(x.g), ncol = sum(l.g[-g])))
    }
    if (g > 1 & g < kG) {
      x.g <- x.list[[g]]
      x.bar <- rbind(x.bar,
     				 cbind(matrix(0, nrow = nrow(x.g), ncol = sum(l.g[1 : (g-1)])), x.g, matrix(0, nrow = nrow(x.g), ncol = sum(l.g[(g + 1) : kG]))))
    }
    if (g == kG) {
      x.g <- x.list[[g]]
      x.bar <- rbind(x.bar,
                     cbind(matrix(0, nrow = nrow(x.g), ncol = sum(l.g[-g])), x.g))
    }
  }
  
  
  # INITIAL PARAMETERS
  if (!test.rf) {
    beta.s <- numeric(length = ncol(x.s))
    beta.bar.star <- numeric(length = ncol(x.bar))
    kappa.eps.0 <- numeric(kG - 1)
    gamma.star.0 <- rep(-5, kG * (kK - 1))
    kappa.u.star <- numeric(length = kG)
    starting.values.no.quad <- c(beta.s = beta.s, beta.bar.star = beta.bar.star, kappa.eps.0 = kappa.eps.0, gamma.star.0 = gamma.star.0)
  }
  
  if (test.rf) {
  	Gamma.d.star.0 <- numeric(length = kG - 1)
  	Gamma.bar.star <- numeric(length = ncol(x.bar))
  	kappa.eta.0 <- numeric(length = kG - 1)
  	gamma.0 <- rep(0.1, kK - 1)
  	kappa.v.star <- numeric(length = kG) - 0.5
  	starting.values.no.quad <- c(Gamma.d.star.0, Gamma.bar.star, kappa.eta.0, gamma.0)
  }

  
  # POINTS FOR THE QUADRATURE
  if (H > 1) {
    v.z.first <- hermite.h.quadrature.rules(10) ; v.z.first <- v.z.first[[length(v.z.first)]]  # points for the first round of quadrature
    v.z <- hermite.h.quadrature.rules(H) ; v.z <- v.z[[length(v.z)]]  # points for the last round of quadrature
  }
  
  
  # CREATING MATRICES FOR LL MAX WITH QUADRATURE (list.mat.quad.sum, list.mat.quad.prod)
  if (H > 1) {
  	# list.mat.quad.sum
  	n.quad.mat.sum <- numeric(Gsplit) + floor(kN / Gsplit)  # vector that says what is the number of individuals to be expected in each sub.matrix

  	n.quad.mat.sum[Gsplit] <- n.quad.mat.sum[Gsplit]  + kN - sum(n.quad.mat.sum)
  	list.mat.quad.sum <- list()
  	for (gs in 1 : Gsplit) {
  	  list.mat.quad.sum[[gs]] <- 1 / sqrt(pi) * diag(n.quad.mat.sum[gs]) %x% matrix(c(v.z[, 2]), nrow = H)  # See "quadrature matrix" sheet
  	}
  	list.mat.quad.sum[[Gsplit + 1]] <- cbind(c(1, 1 + cumsum(n.quad.mat.sum)[-Gsplit] * H), cumsum(n.quad.mat.sum) * H)  # Two-col matrix of indexes
	
	# list.mat.quad.prod
	numeric.indiv.index <- indiv.index ; levels(numeric.indiv.index) <- (1 : nlevels(indiv.index))
	numeric.indiv.index <- as.numeric(numeric.indiv.index)  # Creating a numeric version of indiv.index, that goes from 1 to N
    list.mat.quad.prod <- list() ; list.mat.quad.prod[[Gsplit + 1]] <- numeric()  ; count.indiv <- 1 ; count.row <- 1
  	for (gs in 1 : Gsplit) {
  	  kNTGs <- sum(table(numeric.indiv.index)[count.indiv : (count.indiv + n.quad.mat.sum[gs] - 1)])  # Number of observations treated here
      list.mat.quad.prod[[gs]] <- matrix(0, nrow = kNTGs, ncol = n.quad.mat.sum[gs])
      list.mat.quad.prod[[gs]][cbind(1 : kNTGs, numeric.indiv.index[count.row : (count.row + kNTGs - 1)] - count.indiv + 1)] <- 1
      list.mat.quad.prod[[Gsplit + 1]] <- rbind(list.mat.quad.prod[[Gsplit + 1]], cbind(count.row, count.row + kNTGs - 1)) # Index matrix
      count.indiv <- count.indiv + n.quad.mat.sum[gs]  # Next indiv for the computation
      count.row <- count.row + kNTGs  # Next row in indiv.index
    }
  }
  
  # LIKELIHOOD MAXIMIZATION
  # Proceeding in several round. First round is without . If H = 1 (no quadrature), then this is the last round as well.
  # Else, second round with round(H/3) points for the quadrature, starting from the estimates obtained in the first round.
  # Last round is with all the H points, starting from the estimates of the previous round.
  if (!re.only) {  # If H > 1 and init.para have been given, skip this step
  if (H == 1 & !missing(init.para)) starting.values.no.quad <- init.para
    optim.results.no.quad <- maxLik(logLik = AncRegLlCompound, grad = NULL, hess = NULL, start = starting.values.no.quad, finalHessian = FALSE, 
                            iterlim = 2e3, method = "BFGS", print.level = 2,
                            y = y, x.s = x.s, x.bar = x.bar, group.index = group.index, indiv.index = indiv.index, kG = kG, kK = kK, kN = kN, H = 1, 
                            v.z = v.z, n.quad.mat.sum = NULL, list.mat.quad.sum = NULL, Gsplit = Gsplit, list.mat.quad.prod = NULL, test.rf = test.rf)
                              # assign("test.grad", t(optim.results.no.quad$gradientObs) %*% optim.results.no.quad$gradientObs, envir = .GlobalEnv)
                              # assign("test.OPG", optim.results.no.quad$gradientObs, envir = .GlobalEnv)
  } 
  if (H > 1) {
  	if (missing(init.para) & !test.rf) starting.values <- c(optim.results.no.quad$estimate, kappa.u.star)  # if no init para, use para from previous
  																										   # estimation
  	if (missing(init.para) & test.rf) starting.values <- c(optim.results.no.quad$estimate, kappa.v.star)
  	if (re.only) starting.values <- init.para  # if init para, and if not meant to be used for first reg (H !=1), then use them
    # optim.results.2nd.round <- maxLik(logLik = AncRegLlCompound, grad = NULL, hess = NULL, start = starting.values, finalHessian = FALSE, 
                            # iterlim = 2e+3, method = "BFGS", print.level = 2,
                            # y = y, x.s = x.s, x.bar = x.bar, group.index = group.index, indiv.index = indiv.index, kG = kG, kK = kK, kN = kN,  H = 10, 
                            # v.z = v.z.first)
  	# starting.values <- optim.results.2nd.round$estimate                            
    optim.results.quad <- maxLik(logLik = AncRegLlCompound, grad = NULL, hess = NULL, start = starting.values, finalHessian = FALSE, 
                            iterlim = 2e+3, method = "BFGS", print.level = 2,
                            y = y, x.s = x.s, x.bar = x.bar, group.index = group.index, indiv.index = indiv.index, kG = kG, kN = kN, kK = kK, H = H, 
                            v.z = v.z, n.quad.mat.sum = n.quad.mat.sum, list.mat.quad.sum = list.mat.quad.sum, Gsplit = Gsplit, 
                            list.mat.quad.prod = list.mat.quad.prod, test.rf = test.rf)
  }
  
  if (test.rf) {
  	# Getting back the structural parameters and var-cov matrix for the rf model. H expected to be > 1
  	structural.para.rf <- StructuralPara(coef = optim.results.quad$estimate, kG = kG, kK = kK, kS = kS, l.g = l.g, quad = TRUE, test.rf = TRUE)
  	var.cov.star.rf <- try(solve(t(optim.results.quad$gradientObs) %*% optim.results.quad$gradientObs))
  	if ((class(var.cov.star.rf) == "try-error")) structural.var.cov.rf <- matrix(NA, nrow = length(optim.results.quad$estimate), ncol = length(optim.results.quad$estimate))
  	if (!(class(var.cov.star.rf) == "try-error")) {
  	  gradient <- jacobian(func = StructuralPara, x = optim.results.quad$estimate,
					       kG = kG, kK = kK, kS = kS, l.g = l.g, quad = TRUE, test.rf = TRUE)
      structural.var.cov.rf <- gradient %*% var.cov.star.rf %*% t(gradient)
    }
  	
  	return(list(optim.results.rf = optim.results.quad, structural.para.rf = structural.para.rf, structural.var.cov.rf = structural.var.cov.rf, list(kK = kK, kG = kG, kN = kN, l.g = l.g, kS = kS, H = H)))
  }

  
  # STRUCTURAL PARAMETERS AND VAR-COV MATRICES
  if (!re.only) structural.para.no.quad <- StructuralPara(coef = optim.results.no.quad$estimate, kG = kG, kK = kK, kS = kS, l.g = l.g, quad = FALSE)
  if (H > 1) structural.para.quad <- StructuralPara(coef = optim.results.quad$estimate, kG = kG, kK = kK, kS = kS, l.g = l.g, quad = TRUE)
  
  # Calculating the var-cov matrices (star, i.e. reduced formed, and then normal, i.e. structural)
  if (!re.only) var.cov.star.no.quad <- try(solve(t(optim.results.no.quad$gradientObs) %*% optim.results.no.quad$gradientObs))
  if (H > 1) var.cov.star.quad <- try(solve(t(optim.results.quad$gradientObs) %*% optim.results.quad$gradientObs))
  # Var-cov for the case without quadrature
  if (!re.only) {
  	if ((class(var.cov.star.no.quad) == "try-error")) structural.var.cov.no.quad <- matrix(NA, nrow = length(optim.results.no.quad$estimate), ncol = length(optim.results.no.quad$estimate))
  	if (!(class(var.cov.star.no.quad) == "try-error")) {
  	  gradient <- jacobian(func = StructuralPara, x = optim.results.no.quad$estimate,
					     kG = kG, kK = kK, kS = kS, l.g = l.g, quad = FALSE)
      structural.var.cov.no.quad <- gradient %*% var.cov.star.no.quad %*% t(gradient)
    }
  }
  # With quadrature
  if (H > 1) {
  	if ((class(var.cov.star.quad) == "try-error")) structural.var.cov.no.quad <- matrix(NA, nrow = length(optim.results.quad$estimate), ncol = length(optim.results.quad$estimate))
  	if (!(class(var.cov.star.quad) == "try-error")) {
  	  gradient <- jacobian(func = StructuralPara, x = optim.results.quad$estimate,
					       kG = kG, kK = kK, kS = kS, l.g = l.g, quad = TRUE)
      structural.var.cov.quad <- gradient %*% var.cov.star.quad %*% t(gradient)
    }
  }

  # NAMING
  # Retrieving the names of the structural parameteters
  coef.names <- colnames(x.s)
  for (g in 1 : kG) {
  	coef.names <- c(coef.names, paste("Group", g, colnames(x.list[[g]])))
  }
  for (g in 2 : kG) {
  	coef.names <- c(coef.names, paste("Group", g, "sigma.eps"))
  }  
  for (g in 1 : kG) {
  	coef.names <- c(coef.names, paste("Group", g, 
  									  paste("tau", 1 : (kK - 1))))
  }
  # Naming the structural coef vector and the variance-covariance matrix
  if (!re.only) {
  	names(structural.para.no.quad) <- coef.names
    rownames(structural.var.cov.no.quad) <- colnames(structural.var.cov.no.quad) <- coef.names
  }
  if (H > 1) {
    for (g in 1 : kG) {  	
  	  coef.names <- c(coef.names, paste("Group", g, "sigma.u"))
    }
    names(structural.para.quad) <- coef.names
    rownames(structural.var.cov.quad) <- colnames(structural.var.cov.quad) <- coef.names
  }
  

  # TABLE
  if (H == 1 & !re.only) {
  	table <- cbind(structural.para.no.quad, sqrt(diag(structural.var.cov.no.quad)), structural.para.no.quad / sqrt(diag(structural.var.cov.no.quad)))
    rownames(table) <- coef.names
    colnames(table) <- c("coef", "sd", "t-stat")
  }
  if (H > 1 & !re.only) {
  	table.no.quad <- cbind(structural.para.no.quad, sqrt(diag(structural.var.cov.no.quad)), 
  				   structural.para.no.quad / sqrt(diag(structural.var.cov.no.quad)))
    table.quad <- cbind(structural.para.quad, sqrt(diag(structural.var.cov.quad)), 
  				   structural.para.quad / sqrt(diag(structural.var.cov.quad)))
  	table <- cbind(rbind(table.no.quad, matrix(NA, nrow = kG, ncol = 3)),
  				   table.quad)
    rownames(table) <- coef.names
    colnames(table) <- c("no quad : coef", "sd", "t-stat", "quad : coef", "sd", "t-stat")    
  }
  if (re.only) {
    table.quad <- cbind(structural.para.quad, sqrt(diag(structural.var.cov.quad)), 
  				   structural.para.quad / sqrt(diag(structural.var.cov.quad)))
  	table <- table.quad
    rownames(table) <- coef.names
    colnames(table) <- c("quad : coef", "sd", "t-stat")    
  }
  
  
  # RETURN
  print(table)
  if (H == 1 & !re.only) return(list(optim.results = list(optim.results.no.quad = optim.results.no.quad), 
                         structural.para = list(structural.para.no.quad = structural.para.no.quad),
                         structural.var.cov = list(structural.var.cov.no.quad = structural.var.cov.no.quad),
                         table = table, constants = list(kK = kK, kG = kG, kN = kN, l.g = l.g, kS = kS, H = H)))
  if (H > 1 & !re.only) return(list(optim.results = list(optim.results.no.quad = optim.results.no.quad, optim.results.quad = optim.results.quad), 
                         structural.para = list(structural.para.no.quad = structural.para.no.quad, structural.para.quad = structural.para.quad),
                         structural.var.cov = list(structural.var.cov.no.quad = structural.var.cov.no.quad, 
                                                   structural.var.cov.quad =  structural.var.cov.quad), 
                         table = table, constants = list(kK = kK, kG = kG, kN = kN, l.g = l.g, kS = kS, H = H)))
  if (re.only) return(list(optim.results = list(optim.results.quad = optim.results.quad), 
                         structural.para = list(structural.para.quad = structural.para.quad),
                         structural.var.cov = list(structural.var.cov.quad =  structural.var.cov.quad), 
                         table = table, constants = list(kK = kK, kG = kG, kN = kN, l.g = l.g, kS = kS, H = H)))

}
