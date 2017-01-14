AncBootIteration <- function(data, index) {
  # Function that takes the data, and makes two estimations: one of the main model, and one of the alternative, "reduced-form" model.
  # The estimates are then returned, concacenated.
  indexed.data <- data[index, ]
  indiv.index <- indexed.data$pid
  group.index <- indexed.data$income_bi2
  
  formula.parents <- list(y = inv_ghql ~ lnUChhincome + age + agesq100 + factor(female) + nchild + factor(married) + factor(jbstat_short) + factor(educ) 
                            + factor(big_region) + factor(doiy4), 
						x.s = ~ death_parents)
  formula.parents.rf <- list(y = inv_ghql ~ death_parents + lnUChhincome + age + agesq100 + factor(female) + nchild + factor(married) +   
                                 factor(jbstat_short) + factor(educ) 
                                  + factor(big_region) + factor(doiy4), 
  						          x.s = ~ 1)
  ancreg.BHPS.parents.re <- AncReg(formula = formula.parents, data = indexed.data, indiv.index = indexed.data$pid, 
							group.index = indexed.data$income_bi2, H = 5)
  ancreg.BHPS.parents.re.rf <- AncReg(formula = formula.parents.rf, data = indexed.data, indiv.index = indexed.data$pid, 
							group.index = indexed.data$income_bi2, H = 5, test.rf = TRUE)
  return(c(ancreg.BHPS.parents.re$structural.para, ancreg.BHPS.parents.re.rf$structural.para.rf))

}
