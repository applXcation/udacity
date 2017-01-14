BootIndividualAnalysis <- function(data, formula, heterosk, initial.estimation = NULL, B = 200, parallel = "multicore", ncpus) {
  # Bootstrap over BootIteration (see this function's description). Checked
  # Args
  #   initial.estimation: if the formula has been estimated on these data before, then initial.estimation of class Chopit can be provided. This allows
  #                       to speed up the first step of boot(), where the model is re-estimated on the same data, by providing the appropriate initial 
  #						  para.
  
  # Library
  library(parallel)
  # Sanity Checks
  if (!is.null(initial.estimation)) {
  	if (class(initial.estimation) != "Chopit") stop("BootIndividualAnalysis: initial.estimation is not of class 'Chopit'")
  	if (initial.estimation$call.arg$heterosk != heterosk) stop("BootIndividualAnalysis: initial.estimation$call.arg$heterosk != heterosk")
  }  
  
  # Creating a list of index vectors, each of them being the index for one interation of the bootstrap
  set.seed(123456789)
  indexes <- replicate(n = B, expr = sample(x = 1 : nrow(data), size = nrow(data), replace = TRUE),
  	        simplify = FALSE)
  
  # Bootstrap
  boot.results <- list()
  boot.results$t <- mclapply(X = indexes, FUN = BootIteration, mc.cores = ncpus,
  						     data = data, formula = formula, heterosk = heterosk, initial.estimation = initial.estimation)
  
  # Proportion of iteration for which a problem occured
  perc.wrong <- prop.table(table(boot.results$t == "BootIteration: error on this iteration"))["TRUE"]
  if (!is.na(perc.wrong)) if (perc.wrong != 0) cat("BootIndividualAnalysis: beware,", perc.wrong, "% of the calculations went wrong (bug in IndividualAnalysisMean and IndividualAnalysisVar)")
  
  # Reformatting the results so that it fits IndividualAnalysisMean/Var
  boot.results$t <- t(as.matrix(data.frame(boot.results$t)))
  boot.results$R <- B

  # Return
  return(boot.results)
}



BootIteration <- function(data, index, formula, heterosk, initial.estimation) {
  # One bootstrap iteration: reestimate the chopit as in chopit.results, except that
  #   - uses data[index, ]
  #   - do NOT uses chopit.results$optim.results$estimate as starting values, because these values may not be acceptable for some resampling of the data
  #   - set NAIVE to true (best way to avoid troubles with initial parameters)
  #   - does not calculate the var.cov matrix
  # Then passes the result through 
  # Returns a vector containing :
  #   the (kK) bootstraped coef of corr between hat.E(y.star.s/x.s, theta.hat) and the average threshold, and the other thresholds
  #   the 2*(kK) associated regression coef
  #   the (kK-1) bootstraped coef of corr between V(y.star.s/x.s, theta.hat) and the average threshold, and the other thresholds
  #   the 2*(kK-1) associated regression coef
 
  indexed.data <- data[index, ]
  # In the case the data are the same as the initial data (i.e. first step of the boot procedure), we can use the estimated parameters obtained in a
  # previous estimation, 'initial.estimation'
  if (all(indexed.data == data, na.rm = TRUE) & !is.null(initial.estimation)) par.init <- initial.estimation$optim.results$estimate
  else par.init <- NULL
  
  chopit.boot <- Chopit(formula = formula, data = indexed.data, naive = TRUE, par.init = par.init, varcov.method = "none", heterosk = heterosk)
  
  # If the estimation went wrong (tremendous coef vector), then impossible to do IndividualAnalysisMean and IndividualAnalysisVar.
  # An handling error is thus proposed for now (to be deleted in a next version of the algo, and to be replaced by a fixing of the Chopit()
  # estimation algorithm).
  individual.analysis.mean.boot <- try(IndividualAnalysisMean(chopit.boot))
  individual.analysis.var.boot <- try(IndividualAnalysisVar(chopit.boot))
  
  if (is.list(individual.analysis.mean.boot) & is.list(individual.analysis.var.boot)) {
    corr.mean.boot <- unlist(individual.analysis.mean.boot$corrlist)
    reg.coef.mean.boot <- unlist(lapply(individual.analysis.mean.boot$reglist, FUN = function(x) x$coefficients[, 1]))
    corr.var.boot <- unlist(individual.analysis.var.boot$corrlist)
    reg.coef.var.boot <- unlist(lapply(individual.analysis.var.boot$reglist, FUN = function(x) x$coefficients[, 1]))
    return(boot.iteration = c(corr.mean.boot, reg.coef.mean.boot, corr.var.boot, reg.coef.var.boot))
  }
  
  else {
  	kK  <- initial.estimation$constants$kK
  	return(boot.iteration  = rep("BootIteration: error on this iteration", times = 3 * kK + 3 * (kK - 1))) 
  }  
  

}



# # BootIndividualAnalysisBak <- function(data, formula, heterosk, initial.estimation = NULL, B = 200, parallel = "multicore", ncpus) {
  # # Old version that uses the boot library
  # # Bootstrap over BootIteration (see this function's description). Checked
  # # Args
  # #   initial.estimation: if the formula has been estimated on these data before, then initial.estimation of class Chopit can be provided. This allows
  # #                       to speed up the first step of boot(), where the model is re-estimated on the same data, by providing the appropriate initial 
  # #						  para.
  
  # # Library
  # library(boot)
  # # Sanity Checks
  # if (!is.null(initial.estimation)) {
  	# if (class(initial.estimation) != "Chopit") stop("BootIndividualAnalysis: initial.estimation is not of class 'Chopit'")
  	# if (initial.estimation$call.arg$heterosk != heterosk) stop("BootIndividualAnalysis: initial.estimation$call.arg$heterosk != heterosk")
  # }  
  
  # # Bootstrap
  # boot.results <- boot(data = data, statistic = BootIteration, formula = formula, heterosk = heterosk, initial.estimation = initial.estimation, R = B, 
                       # stype = "i", parallel = "multicore", ncpus = ncpus)
  
  # # Return
  # return(boot.results)
# }


