NonparaAnalysis <- function(c, v1, ...,  grouping.var, moment, jointly = FALSE, subset.index = NULL) {
  #  Clusters the population by grouping.var ; then checks whether a group moment in variable c varies with its moments in vignette evaluations.
  #  Performs a regression of the groups' moments of variable c on the groups' moments of variable v1 and, if any, v2, v3 etc. (jointly or not). 
  #  Plot this regression, and displays the t-stat.
  #
  # Args:
  #   c: the c variable
  #   v1: first vignettes
  #   ...: additional vignettes, if any
  #   grouping.var: grouping variables
  #   moment: name (in characters) of the moment function to be used (mean, var, etc.)
  #   jointly: should the link between self and vignette evaluations' moment be analyzed on the sample cbind(c(c_group, v1_group), c(c_group,v2)) 
  #   jointly, or on c(c_group, v1_group) and c(c_group, v2_group) separately ?
  #   subset.index: subset of observations on which the analysis should be performed. If NULL (default), all observations are used
  # Returns:
  #   plots the regression result(s)
  #   As a list:
  #     a matrix where the columns indicate the groups, first row gives for each group the moment of c, and the rest of the rows give for each 
  #     group 
  #     the moment of the each vignette variable
  #     as a sublist, the regression results
   
  
  # Error handling
  #   If c contains more than 1 variables, error
  if (!is.null(dim(c))) { 
  	if (((nrow(c) > 1 & ncol(c) > 1))) stop("More than one c variable \n") 
  }
  if (jointly == TRUE) stop("NonparaAnalysis: jointly = TRUE as not be implemented yet \n")
  #   moment should be of type character  
  if (!is.character(moment)) stop("NonparaAnalysis: moment should be of type character")
    
    
  # Vignettes-related objects
  vignettes <- data.frame(v1, ...)  # vignettes is a matrix, where each column is a vignette, and each row is an observation
  M <- ncol(vignettes)
  
  
  # Grouping variable related objects
  grouping.var <- interaction(grouping.var)
  card.groups <- length(levels(as.factor(grouping.var)))
  
  
  # Subsetting
  if (is.null(subset.index)) subset.index <- 1:length(c)  # if nothing specified as a subset.index, all the observations are used
  c <- c[subset.index]
  vignettes <- vignettes[subset.index, , drop = FALSE]
  grouping.var <- grouping.var[subset.index]
#  return(list(c = c, vignettes = vignettes, grouping.var = grouping.var))
  
  # CASE 1: analysis is be performed separately for each vignette statement variable
  if (jointly == FALSE) {  
  	# Computing the moments for each groups
  	groups.moment.c <- rbind(tapply(X = c, INDEX = grouping.var, FUN = get(moment), na.rm = TRUE, simplify = TRUE))
  	groups.moment.v <- numeric(0)
  	for (i in 1:M) {
  	  groups.moment.v <- rbind(groups.moment.v, tapply(X = vignettes[, i], INDEX = grouping.var, FUN = get(moment), na.rm = TRUE, simplify = TRUE))
  	}
#  	return(list(groups.moment.c, groups.moment.v))
  	# Performing the regressions separately, and plotting the results
  	if (M > 1) par(mfrow = c(M / 2 * (M %% 2 == 0) + ceiling(M / 2) * !(M %% 2 == 0), 2))
	if (M == 1) par(mfrow = c(1,1))
  	regression.results <- list()
  	for (i in 1:M) {
  	  reg <- lm(groups.moment.c[, , drop = TRUE] ~ groups.moment.v[i, , drop = TRUE])
  	  summary.reg <- summary(reg)
  	  regression.results[[i]] <- summary.reg
  	  pvalue <- 2 - 2 * (pt(abs(summary.reg$coefficients[2, "t value"]), df = length(summary.reg$residuals) - 2))  # pvalue for slope coef
 	  plot(groups.moment.v[i, , drop = TRUE], groups.moment.c[, , drop = TRUE], main = paste("Vignette", i),
   	       ylab = paste(moment, "C by group"),
   	       xlab = paste(moment, " V", i, " by group", sep = ""), 
  	       sub = paste("p-value:", round(pvalue, digits = 3), "; G =", length(reg$residuals)))
  	  abline(reg = reg)  # Adding the regression line
  	}
  }
  
    
  return(list(moments = rbind(groups.moment.c, groups.moment.v), regression.results = regression.results)) 
}
