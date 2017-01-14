SNonparaAnalysis <- function(self.eval, v1, ..., grouping.var, moment, jointly = FALSE) {
  # Simplified Non-Parametric Analysis. Takes, as an argument, the self-assessment and vignette variables, together with the grouping 
  # variables as a dataframe. Calculates by itself the "c" variable through NonparaComparable(), and then uses NonparaAnalysis to study the link
  # between moments in c and moments in the vignettes.
  # It should be remarked that we only perform the second step of the analysis with the "c"s that are not inconsistent. That is, we only keep the
  # observations for which there is only one c (no multiple c-category in which the observation could belong), et no inconsistency in the ranking
  # of the vignettes if there are multiple vignettes (the individual does not inverse the ranking of the vignettes in his/her evaluations).
  #
  # Args:
  #   self.eval: respondent self-evaluation
  #   v1: first  vignette
  #   ...: other vignettes
  #   Remark: In the list of arguments, the vignettes should be ordered such that the first vignette corresponds to the level of domain closest to 
  #   label "1", and the last vignette to the level of domain closest to label "J"
  #   grouping.var: data.frame containing the grouping variables
  #   moment: name (in characters) of the moment function to be used (mean, var, etc.)
  #
  # Returns:
  #   plots the regression result(s)
  #   As a list:
  #     a matrix where the columns indicate the groups, first row gives for each group the moment of c, and the rest of the rows give for each 
  #     group 
  #     the moment of the each vignette variable
  #     as a sublist, the regression results
  
  # Sanity Checks
  if (!is.data.frame(grouping.var)) stop("SNonparaAnalysis: grouping.var is not a data.frame \ n")
  if (!is.character(moment)) stop("SNonparaAnalysis: moment should be of type character")
  
  # Obtaining the c variable
  c.list <- NonparaComparable(self.eval = self.eval, v1 = v1, ... = ...)  # List of the c variables
  relevant.cases <- !is.na(c.list$c[, 1]) & c.list$multiple == 0 & c.list$inconsistencies == 0  # Index of the relevant observations for the 
                                                                                                # analysis

  # Studying the linear correlation between c and v. This evaluation is performed in the observations for which there is no inconsistency in the
  # c, as described above.
  assign(deparse(substitute), moment)
  non.para.analysis <- NonparaAnalysis(c = c.list$c[, 1], v1 = v1, ... = ..., grouping.var = grouping.var, moment = moment, jointly = jointly, 
                                       subset.index = relevant.cases)
  
  # Return
  return(non.para.analysis)
}