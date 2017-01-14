if (Sys.info()["sysname"] == "Darwin" & Sys.info()["user"] == "xavier") {
  chopit.routines.directory = "~/Documents/Travail/Boulot Fac/Doctorat/1_ Productions personnelles/24_Range-Frequency/2_programmes/"
  share.data.directory = "~/Documents/Travail/Boulot Fac/Doctorat/Data bases/SHARE/0_merged_data/"
}

if (Sys.info()["sysname"] == "Linux" & Sys.info()["user"] == "x.fontaine") {
  chopit.routines.directory = "~/U/Travail/Range_frequency/routines/"
  share.data.directory = "~/U/Travail/Range_frequency/data/"
}


source(paste(chopit.routines.directory,"likelihood.R", sep = ""), chdir = TRUE)
source(paste(chopit.routines.directory,"ChopitStartingValues.R", sep = ""), chdir = TRUE)
source(paste(chopit.routines.directory,"summary.Chopit.R", sep = ""), chdir = TRUE)



Chopit <- function(formula, data, heterosk = FALSE, naive = FALSE, par.init = NULL, optim.method = "BFGS", varcov.method = "OPG") {
  # CHOPIT()
  # Args:
  #   formula: list containing 3 formulae arguments ; the first one must be names self, the second one vign, the last one tau. Remark that a
  #            a constant is expected to be specified in all the equations, even though the program will perform the necessary normalization
  #            beta["cste"] = 0. If heterosk == TRUE, a 4th argument (sigma) is required.
  #   data: dataset associated to the formulae
  #   heterosk: Should the model account for potential heteroskedasticity ? If so, sd(Y*_s / x_s, x_sigma) = sigma_s(x_sigma) = exp(x_sigma %*% kappa)
  #   par.init: initial parameters to be passed to optim
  #   naive: TRUE if initial values for parameters should be set the "naive" way (i.e. all parameters are 0, except intercept of each threshold
  #          equation which is set to 0.1. Else, a more sophisticated way is used, through ChopitStartingValues
  #   varcov.method: method to be used to estimate the variance-covariane matrix for the parameters: none, hessian or OPG
  #
  # Returns (as a list)
  #   optim.results: result from optimization (object of classes  maxLik and maxim) ; for sake of sparsity, gradientObs is delete
  #   coef: list of coefficients. beta0 is the beta vector, with a 0 as a first argument by normalization since there is an intercept in x.s. Same holds
  #         for kappa0
  #   var.cov: variance covariance matrix (if varcov.method != "none")
  #   constants: list of important constants (defined throughout the code)
  #   contrib.number: number of contribution to partial likelihood
  #   call.arg: argument used when the function has been called. 'data.name' and parent.frame replaces the original dataset, by giving both the name of
  #             this dataframe and its location. By doing so, the object returned remains quite small in size.
  #   col.to.drop: if Chopit() finds any colinearity problem, it drops the incriminated collumns. The arguments in col.to.drop contain the number of the
  #                columns to be dropped


  # Libraries
  library(MASS)  # used for finding proper starting values, using polr
  library(maxLik)
  
  # Start the clock
  ptm <- proc.time()
  
    
  # GLOBAL SANITY CHECKS
  # data should be a data.frame
  if (!is.data.frame(data)) stop("Chopit: the argument data is not a data.frame")
  # varcov.method should be either "none", "hessian" or OPG
  if (!varcov.method %in% c("none", "hessian", "OPG")) stop("Chopit: varcov.method should be either 'none', 'hessian' or 'OPG'")  
  
  # PARSING THE FORMULAE
  # Sanity checks: checking the "formula" list is properly defined
  #   Is "formula" a list ?
  if (!is.list(formula)) stop("Chopit: formula is not a list")
  #   Is each element a formula ?
  IsFormula <- function(x) {
  	is.formula <- (class(x) == "formula")
  	return(is.formula)
  }
  if (any(!sapply(X = formula, FUN = IsFormula))) stop("Chopit: at least one element in the list formula is not a formula")  
  #   Are the names apropriately chosen ?
  if (any(!c("self", "tau", "vign") %in% names(formula))) stop("Chopit: names in formula badly specified (one of self, tau, vign is missing)")
  if  (heterosk & !("sigma" %in% names(formula))) stop("Chopit: 'heterosk == TRUE' but no 'sigma' in the formula")
  if  (!heterosk & ("sigma" %in% names(formula))) stop("Chopit: heterosk = FALSE while there is 'sigma' in the formula")

  
  # Parsing process
  f.self <- formula$self
  f.tau <- formula$tau
  f.vign <- formula$vign
  if (heterosk) f.sigma <- formula$sigma  

  # PRODUCING THE DATA MATRICES
  # Sanity checks:
  #   A constant is expected to be included in the self-assessment equation (no + 0 or - 1 should take place in this equation). If no constant
  #   specified, we should force a constant to exist (cf. terms.object so see what to modify), and indicate we do so.
  if (attr(terms(f.self), "intercept") == 0) stop("No constant in the self-assessment equation formula (one is expected, event though we normalize 
  the associated coef to 0)")
  
  # Getting the name of the provided argument 'data' before data is evaluated
  data.name <- deparse(substitute(data))
  # Dropping unused levels in data if any remaining
  data <- droplevels(data)


# #   # Self-assessment. cbind() is used to make sure each object is at least a column-vector (else a matrix). original.objects are created to
  # # keep the objects as they were with NA, so that it can be returned by Chopit() and eventually passed to other functions (GroupAnalysis).
  # # Indeed, whenever some values in x.tau is missing, then observation in x.s is deleted here ; but these observations can be used in 
  # # GroupAnalysis(). Remark that variables dropped due to multi-co here are also going to be dropped in the "original" objects.
  # mf.self <- model.frame(formula = f.self, data = data, na.action = NULL)
  # y.s <- model.response(data = mf.self) ; y.s <- cbind(y.s)
  
  # # Vignettes
  # mf.vign <- model.frame(formula = f.vign, data = data, na.action = NULL)
  # y.v <- model.response(data = mf.vign) ; y.v <- cbind(y.v)
  
  # # Checking again for missing levels, but now when all y.s and y.v are NA. Otherwise, some levels may not be missing for the whole dataset, but
  # # missing when considering only the observations where we have at least one assessment avaible.
  # # Then getting back y.s and y.v again for this smaller dataset
  # data <- data[!is.na(y.s) | rowSums(!is.na(y.v)), ]
  # data <- droplevels(data)
  # mf.self <- model.frame(formula = f.self, data = data, na.action = NULL)
  # y.s <- model.response(data = mf.self) ; y.s <- cbind(y.s)
  # mf.vign <- model.frame(formula = f.vign, data = data, na.action = NULL)
  # y.v <- model.response(data = mf.vign) ; y.v <- cbind(y.v)
  
  # # Self-assessment: X
  # mf.self <- model.frame(formula = f.self, data = data, na.action = NULL)
  # x.s <- model.matrix(object = f.self, data = mf.self) ; x.s <- cbind(x.s)  

  # # Tau
  # mf.tau <- model.frame(formula = f.tau, data = data, na.action = NULL) ; 
  # x.tau <- model.matrix(object = f.tau, data = mf.tau) ; x.tau <- cbind(x.tau)


  
  # Self-assessment. cbind() is used to make sure each object is at least a column-vector (else a matrix). original.objects are created to
  # keep the objects as they were with NA, so that it can be returned by Chopit() and eventually passed to other functions (GroupAnalysis).
  # Indeed, whenever some values in x.tau is missing, then observation in x.s is deleted here ; but these observations can be used in 
  # GroupAnalysis(). Remark that variables dropped due to multi-co here are also going to be dropped in the "original" objects.
  mf.self <- model.frame(formula = f.self, data = data, na.action = NULL)
  y.s <- model.response(data = mf.self) ; y.s <- cbind(y.s)
  x.s <- model.matrix(object = f.self, data = mf.self) ; x.s <- cbind(x.s)
  
  # Tau
  mf.tau <- model.frame(formula = f.tau, data = data, na.action = NULL) ; 
  x.tau <- model.matrix(object = f.tau, data = mf.tau) ; x.tau <- cbind(x.tau)
  
  # Vignettes
  mf.vign <- model.frame(formula = f.vign, data = data, na.action = NULL)
  y.v <- model.response(data = mf.vign) ; y.v <- cbind(y.v)
  
  # Heteroskedasticity
  if (heterosk) {
  	mf.sigma <- model.frame(formula = f.sigma, data = data, na.action = NULL)
    x.sigma <- model.matrix(object = f.sigma, data = mf.sigma) ; x.sigma <- cbind(x.sigma)
  }
  else {
  	x.sigma <- NULL 
  	original.x.sigma <- NULL
  }
  
  # DEALING WITH NA
  # Observations that cannot be used at all:
  #   - Rows for which there is no self-assessment AND vignette information
  #   - Rows for which one of the x.tau is missing (impossible to calculate the thresholds)
  # Observations that cannot be used to calculate individual contribution to the self-assessment question likelihood (could be used for vignettes)
  #   - Row for which one of the x.s is missing
  #   - Rows for which y.s is missing
  # Observations that cannot be used to calculate indiv contrib to ONE vignette v question likelihood
  #   - Rows for which the vignette statement is missing
  # Here, I drop all the data for which one or both of the first two conditions are not met. This avoids calculating likelihood contributions for 
  # observations that are going to be NA. The reste of the NA are handled through the use of na.rm = TRUE in the sum functions of the ll function.
  # Alternatives could have been used, but this way of handling NA offers a good trade-off between speed and generality.
  #
  # Deleting rows corresponding to the first two cases
  #   First, detecting the incomplete cases
  any.non.na <- function(x) {  # Function that takes a vector, and says whether there is any non-missing value
  	any(!is.na(x))
  }
  incomplete <- ((is.na(y.s)) & !(apply(y.v, 1, any.non.na))) | (!complete.cases(x.tau)) 
  #   Then, deleting the incriminated rows
  y.s <- y.s[!incomplete, , drop = FALSE]
  x.s <- x.s[!incomplete, , drop = FALSE]
  x.tau <- x.tau[!incomplete, , drop = FALSE]
  y.v <- y.v[!incomplete, , drop = FALSE]
  if (heterosk) x.sigma <- x.sigma[!incomplete, , drop = FALSE]
  
  
  # DEALING WITH PERFECT COLINEARITY
  # Detecting perfect colinearity through polr, displaying the name of the incriminating variables, and dropping the corresponding variables
  # We focus on observations for which we observe at the same time the self-assessment variable and any vignette assessment. This allows us to check 
  # that we observe each variable when, and especially that none of these variables is a constant (e.g. dummies not-used)
  s.e.and.v.e <- !is.na(y.s) & rowSums(1 * !is.na(y.v))  # Obervations for which we observe both a self-eval and a vignette eval

  #   x.s
  if (ncol(x.s) > 2) {  # Dealing with multico matters only if there is more than the constant in x.s
    temp.polr <- polr(as.factor(y.s[s.e.and.v.e]) ~ x.s[s.e.and.v.e, -1] , method = "probit")
    proper.names <- gsub("x.s\\[s.e.and.v.e, -1\\]", "", names(temp.polr$coefficients))
    col.to.drop.x.s <- !(colnames(x.s) %in% proper.names) ; col.to.drop.x.s[1] <- FALSE  # Index of the columns in x.s to be dropped
    if (any(col.to.drop.x.s)) {
      cat("Chopit: x.s is not full-rank. Dropping the following columns:", colnames(x.s)[col.to.drop.x.s], "\n")
      x.s <- x.s[, !col.to.drop.x.s, drop = FALSE]
    }
    rm(temp.polr)
  }
  #  x.tau
  col.to.drop.x.tau <- NULL
  if (ncol(x.tau) > 2) {  # Dealing with multico matters only if there is more than the constant in x.tau
    temp.polr <- polr(as.factor(y.s[s.e.and.v.e]) ~ x.tau[s.e.and.v.e, -1] , method = "probit")
    proper.names <- gsub("x.tau\\[s.e.and.v.e, -1\\]", "", names(temp.polr$coefficients))
    col.to.drop.x.tau <- !(colnames(x.tau) %in% proper.names) ; col.to.drop.x.tau[1] <- FALSE  # Index of the columns in x.tau to be dropped
    if (any(col.to.drop.x.tau)) {
      cat("Chopit: x.tau is not full-rank. Dropping the following columns:", colnames(x.tau)[col.to.drop.x.tau], "\n")
      x.tau <- x.tau[, !col.to.drop.x.tau, drop = FALSE]
    }
    rm(temp.polr)
  }
  #  x.sigma
  col.to.drop.x.sigma <- NULL
  if (heterosk) if (ncol(x.sigma) > 2) {  # Dealing with multico matters only if there is more than the constant in x.sigma
    temp.polr <- polr(as.factor(y.s[s.e.and.v.e]) ~ x.sigma[s.e.and.v.e, -1] , method = "probit")
    proper.names <- gsub("x.sigma\\[s.e.and.v.e, -1\\]", "", names(temp.polr$coefficients))
    col.to.drop.x.sigma <- !(colnames(x.sigma) %in% proper.names) ; col.to.drop.x.sigma[1] <- FALSE  # Index of the columns in x.sigma to be dropped
    if (any(col.to.drop.x.sigma)) {
      cat("Chopit: x.sigma is not full-rank. Dropping the following columns:", colnames(x.sigma)[col.to.drop.x.sigma], "\n")
      x.sigma <- x.sigma[, ! col.to.drop.x.sigma, drop = FALSE]
    }
    rm(temp.polr)
  }  
   

   
    
  # CONSTANTS
  # Calculating the number of statement categories, by taking the maximum among the self-assessements and the vignettes
 length.levels.as.factor <- function(x) {  # Function that calculate the length of the levels of an objet passed to as.factor
  	length(levels(as.factor(x)))
  }
  kK <- max(sapply(X = data.frame(cbind(y.s, y.v)), FUN = length.levels.as.factor))  # Number of statement categories
  kBeta0Nrow <- ncol(x.s) - 1  # Number of parameters in beta0
  kGammaNrow <- ncol(x.tau)  # Number of gamma parameters for each threshold equation
  kV <- ncol(y.v)  # Number of vignettes
  if (heterosk) {
  	kKappa0Nrow <- ncol(x.sigma) - 1  # Number of row in kappa (except the 0 for the intercept)
  }
  else {
    kKappa0Nrow <- NULL  
  }
  

  # GENERATING STARTING VALUES
  if (naive == TRUE & is.null(par.init)) {
    beta0.init <- numeric(kBeta0Nrow)
    gamma.init <- matrix(0, nrow = kGammaNrow, ncol = kK - 1) ; gamma.init[1, ] <- 0.1
    theta.init = numeric(kV)
    sigma.tilde.v.init <- numeric(kV)  # equivalent to setting sigma. v = (1, ..., 1)
    if (heterosk) {
      kappa0.init <- numeric(kKappa0Nrow)  # parameters kappa, except the first 0
    }
    if (heterosk) par.init <- c(beta0.init, gamma.init, theta.init, sigma.tilde.v.init, kappa0.init)
    else par.init <- c(beta0.init, gamma.init, theta.init, sigma.tilde.v.init)
  }
  if (naive == FALSE & is.null(par.init)) {
    par.init <- ChopitStartingValues(y.s = y.s, x.s = x.s, kK = kK, kBeta0Nrow = kBeta0Nrow, kGammaNrow, kV)
    if (heterosk) {             
      par.init <- c(par.init, numeric(kKappa0Nrow))
      cat("Chopit: The non-naive parameter initilalization is not optimized for being used when heteroskedasticticy is allowed. Naive could be 
           prefered. \n")
    }
  }  
  

  # LIKELIHOOD MAXIMIZATION
  chopit.envir <- environment()
#  optim.results <- optim(par = par.init, fn = ChopitLlCompound, y.s = y.s, y.v = y.v, x.s = x.s, x.tau = x.tau, kK = kK, kBeta0Nrow = kBeta0Nrow, 
#        kGammaNrow = kGammaNrow, kV = kV, method = optim.method, control = list(trace = 2), hessian = varcov.calculus)
  optim.results <- maxLik(logLik = ChopitLlCompound, grad = NULL, hess = NULL, start = par.init, finalHessian = (varcov.method == "hessian"), 
                          iterlim = 2e+3, method = optim.method, print.level = 2, y.s = y.s, y.v = y.v, x.s = x.s, x.tau = x.tau, kK = kK, 
                          kBeta0Nrow = kBeta0Nrow, kGammaNrow = kGammaNrow, kV = kV, chopit.envir = chopit.envir,
                          heterosk = heterosk, x.sigma = x.sigma, kKappa0Nrow = kKappa0Nrow)

  # VAR-COV MATRIX
  if (varcov.method == "none") {
  	var.cov <- matrix(NA, nrow = length(par.init), ncol = length(par.init))
  }
  if (varcov.method == "hessian") {
  	var.cov <- - solve(optim.results$hessian)
  }  
  if (varcov.method == "OPG") {
  	var.cov <- solve(t(optim.results$gradientObs) %*% optim.results$gradientObs)
  }    
  
  
  # NAMING
  # Naming the rows (and cols) of the estimated parameters and of the varcov matrix
  #   Creating a vector of names
  beta0.names <- colnames(x.s[, -1])
  gamma.names <- vector("numeric")
  for (i in 1:(kK - 1)) {
  	gamma.names <- cbind(gamma.names, paste(paste("gamma", i, sep = ""), colnames(x.tau)))
  }
  theta.names <- paste("theta", 1:kV, sep = "")
  sigma.tilde.v.names <- paste("sigma.tilde.v", 1:kV, sep = "")
  if (heterosk) kappa.names <- paste("kappa", colnames(x.sigma[, -1]))
  else kappa.names <- NULL
  names <- c(beta0.names, gamma.names, theta.names, sigma.tilde.v.names, kappa.names)
  # Renaming the appropriate objects
  names(optim.results$estimate) <- names
  rownames(var.cov) <- colnames(var.cov) <- names  
  
  
  # DEPARSING COEF
  # Deparsing the coefficients in sub-categories to facilitate further uses
  beta0 <- optim.results$estimate[1:kBeta0Nrow]
  gamma <- matrix(optim.results$estimate[(length(beta0) + 1) : (length(beta0) + kGammaNrow * (kK - 1))], ncol = kK - 1)
  theta <- optim.results$estimate[(length(beta0) + length(gamma) + 1) : (length(beta0) + length(gamma) + kV)]
  sigma.tilde.v <- optim.results$estimate[(length(beta0) + length(gamma) + length(theta) + 1) : (length(beta0) + length(gamma) + length(theta)
                                          + kV)]  
  if (heterosk) kappa0 <- optim.results$estimate[(length(beta0) + length(gamma) + length(theta) + length(sigma.tilde.v) + 1) : (length(beta0) + length(gamma) + length(theta) + length(sigma.tilde.v) + kKappa0Nrow)] 
  else kappa0 <- NULL
  
  
  # Switching the clock off
  elapsed.time <- proc.time() - ptm
  
  
  # RETURN
  optim.results$gradientObs <- NULL  # Dropping gradientObs for the sake of sparsity
  results <- list(optim.results = optim.results, coef = list(beta0 = beta0, gamma = gamma, theta = theta, sigma.tilde.v =sigma.tilde.v, 
                 									         kappa0 = kappa0), 
                 var.cov = var.cov, constants = list(kK = kK, kBeta0Nrow = kBeta0Nrow, kGammaNrow = kGammaNrow, kV = kV, 
                                                     kKappa0Nrow = kKappa0Nrow,  heterosk = heterosk), 
                 contrib.number = contrib.number,
                 call.arg = list(formula = formula, data.name = data.name, heterosk = heterosk, naive = naive, par.init = par.init, 
                                 optim.method = optim.method, varcov.method = varcov.method, parent.frame = parent.frame()),
                 col.to.drop = list(x.s = col.to.drop.x.s, x.tau = col.to.drop.x.tau, x.sigma = col.to.drop.x.sigma),
                 elapsed.time = elapsed.time)
  class(results) <- "Chopit"
  return(results)
}



# "Compiling" the code to make it faster
library(compiler)
ChopitTau <- cmpfun(ChopitTau)
ChopitLlSelfEval <- cmpfun(ChopitLlSelfEval)
ChopitLlVignEval <- cmpfun(ChopitLlVignEval)
ChopitLl <- cmpfun(ChopitLl)
ChopitLlCompound <- cmpfun(ChopitLlCompound)
Chopit <- cmpfun(Chopit)


