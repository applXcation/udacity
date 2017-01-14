summary.Chopit <- function(results) {
  # Creates summary objects for objects of class Chopit.
  #
  # Args
  #   results: object of class "Chopit", i.e. result of call to function Chopit()
  #
  # Returns
  #   a list of class summary.Chopit containing the following elements
  #     table: table of summary statistics for the coef
  #     contrib.number: table displaying the number of contrib to each sub-likelihood

  # Building the coef / stat table
  if (all(is.na(results$var.cov))) {
  	table <- cbind(Coef = results$optim.results$estimate, StdErr = NA, z.stat = NA, p.value = NA)
  }
  else {
    se <- sqrt(diag(results$var.cov))
    z.stat <- results$optim.results$estimate / se
    p.value <- 2 * pnorm(- abs(z.stat))
    table <- cbind(Coef = results$optim.results$estimate, StdErr = se, z.stat = z.stat, p.value = p.value)
  }
  
  # Building the table displaying the number of contrib to each sub-likelihood
 contrib.number <- cbind(unlist(results$contrib.number))
 rownames(contrib.number) <- c("Self-evaluation", paste("Vignette evaluation", 1:(length(contrib.number) - 1)))
  
 #Return
 summary <- list(table = table, contrib.number = contrib.number) ; class(summary) <- "summary.Chopit"
 return(summary) 
}


print.summary.Chopit <- function(summary) {
	cat("Coef table: \n")
	printCoefmat(summary$table, P.values = TRUE, has.Pvalue = TRUE)
	cat("\n Number of individual contributions to each sub-part of the likelihood \n")
	print(summary$contrib.number)
	cat("\n")
}