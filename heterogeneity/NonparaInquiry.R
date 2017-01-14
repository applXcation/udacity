NonparaInquiry <- function(y.s, y.v) {
# See the "test non-para" sheets

  max.k <- max(y.s, y.v, na.rm = TRUE)
  table <- matrix(NA, nrow = 2, ncol = max.k)
  rownames(table) <- c("y.v > k", "y.v <= k") ; colnames(table) <- 1:max.k  # temporary colnames for table
  for (k in 1 : max.k) {
    indic.y.v.k <- (y.v <= k)
    p.y.s.0 <- sum(y.s[!indic.y.v.k] <= k, na.rm = TRUE) / sum(!is.na(y.s[!indic.y.v.k]))  # Probability that y.s <= k given that !(y.v <= k)
    p.y.s.1 <- sum(y.s[indic.y.v.k] <= k, na.rm = TRUE) / sum(!is.na(y.s[indic.y.v.k]))  # Probability that y.s <= k given that (y.v <= k)
    table[, k] <- c(p.y.s.0, p.y.s.1)
    colnames(table)[k] <- paste("P(y.s <= ", k, ")", sep = "")
  }
  return(table)
}
