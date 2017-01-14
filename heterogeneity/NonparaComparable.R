##############################################
# NON PARAMETRICALLY MADE COMPARABLE VARIABLE
##############################################
# Main function: NonparaComparable
# Sub-functions: Ci and Inconsistenciesi
# Description: computes the non-parametric variable as described in King & al (2004, APSR)


      
Ci <- function(xi, M) {
  # Takes x <- cbind(self.evaluation,vignette1, vignette2) as an argument. This function looks at the position of x[i,1] w.r.t. x[i,2], ..., x[i,M] and returns a vector of values of ci, with ci[j]=j-th value that c takes if any, NA else
  # Args:
  #   xi, vector of observation (self.eval,v1,...) for individual i
  #   M, the number of vignettes
  # Returns:
  #   ci, a numeric vector. ci is made of NA if there is any NA is i's self-evaluation of vignette evaluations. Else, entry j is equal to the label 
  #   of the j-th condition satisfied for this individual
 
  ci <- vector(mode = "numeric", length = 2 * M + 1) * NA  # initialization of ci
  if (any(is.na(xi) == TRUE)) {  # if any of these xi is na, then return a vector of NA for this individual
  	return(ci)	
  }
  j = 1  # initialization of j, the column for the next entry
  if (xi[1] < xi[2]) {
    ci[j] <- 1
  	j <- j + 1
  }
  if (xi[1] == xi[2]) {
   	ci[j] <- 2
   	j <- j + 1
  }
  if (M > 1) {
  	for (k in 1 : (M - 1)) {
      if (xi[k + 1] < xi[1] & xi[1] < xi[k + 2]) {
        ci[j] <- 2 * k + 1
        j <- j + 1
      }
      if (xi[1] == xi[k + 2]) {
        ci[j] <- 2 * k + 2
        j <- j + 1
      }
    }
  }
  if (xi[M + 1] < xi[1]) {
  	  ci[j] <- 2 * M + 1
  	  j <- j + 1
  }
  return(ci)
}



NonparaComparable <- function(self.eval, v1, ...) {
  # Args:
  #   self.eval is the respondent's self evaluation on the domain
  #   ... is a set of vignette evaluations for the same domain. These vignettes should be provided as separated vectors.
  #   Remark: 
  #     In the list of arguments, the vignettes should be ordered such that the first vignette corresponds to the level of domain closest to 
  #     label "1", and the last vignette to the level of domain closest to label "K"
  # Returns:
  #   c: matrix (c1,...,c(max)), where Cj is the j-th answer category the observation belongs to, and severalC says whether the observation belongs 
  #   to one category or more.
  #   multiple: logical vetor. Entry i is TRUE if there are two or more ("non-parametric") categories in which self-evaluation and vignettes 
  #   evaluations of individual i enters ; FALSE else.
  #   inconsistencies: logical vector, of which the i-th element is set to TRUE if individual i exhibits any inconsistency in its ranking, FALSE 
  #   else

  # Error handling
  #   If arguments are not numeric, error
  if (is.numeric(cbind(v1, ...)) == FALSE) stop("The arguments are not numeric")
  
  # Displaying informations
  cat("the vignettes should be ordered such that the first vignette corresponds to the level of domain closest to label `1', and the last vignette 
       to the level of domain closest to label `K'")
  
  # Basic computations
  self.eval <- cbind(self.eval)  # self.eval becomes a column vector
  v.mat <- cbind(v1, ...)  # all the vignette evaluation vectors are concatenated into a matrix
  M <- ncol(v.mat)  # Number of vignettes
  K <- max(cbind(self.eval, v.mat),na.rm = TRUE)  # Number of answers category ; maximum of the number of categories among self.eval,v1, ...
  
  
  # Computing c using the function Ci()
  x <- cbind(self.eval,v.mat)  # We concatenate the whole to perform a tapply on it 
  c <- apply(X = x, MARGIN = 1, FUN = Ci, M = M)
  c <- t(c)
   
  
  # Computing a vector for which entry i is TRUE if there are two or more ("non-parametric") categories in which self-evaluation and vignettes 
  # evaluations of individual i enters ; FALSE else.
  multiple <- !is.na(c[, 2])
  
  # Computing a logical vector, of which the i-th element is set to TRUE if individual i exhibits any inconsistency in its ranking, FALSE else.
  # If only 1 vignette (M = 1), all the elements are obviously set to FALSE
  # If M > 1, for all the individuals, we check whether UNION_{i=1}^{M-1} [UNION_{j=i+1}^{M} (v_i > v_j)]  ; if so, inconsistencies = TRUE.
  inconsistencies <- cbind(logical(nrow(v.mat))) #  Initialization as a column vector (all values = FALSE)
  if (M > 1) {
    for (i in (1 : (M - 1))) {
      for (j in ((i + 1) : M)) {
        inconsistencies <- inconsistencies | ((v.mat[, i] > v.mat[, j]))
      }
    }
  }
  inconsistencies[is.na(inconsistencies)] <- FALSE  # if NA, inconsistencies is set to FALSE to stick to the previous definition of inconsistencies
  
  #Return
  return(list(c = c, multiple = multiple, inconsistencies = inconsistencies))
}
