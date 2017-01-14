TableStructReduced <- function(ancreg.results.structural, ancreg.results.rf, names) {
  # TableStructReduced
  # Arg:
  #   ancreg.results.structural: ancreg results for the structural model
  #   ancreg.results.rf: ancreg results for the reduced-form model
  #   names: 2-col matrix ; first row are the original names that are to be replaced as in ancreg.object$table (without the "Group g" that preceeds 
  #          them) ; second column are the names one wants to use
  
  
  # PACKAGES
  library(xtable)


  # CONTROL VARIABLES
  missing.names <- missing(names)  # Logical, says whether `names' has been specified

  
  # SANITY CHECKS
  if (!missing.names) {
  	if (!is.matrix(names)) stop("names is not specified as a matrix")
  	if (ncol(names) != 2) stop("names does not have 2 columns")
  }
  
  
  # CONSTANTS
  
  
  # 
  
  
  # Return
  return()
}