TableAncThresh <- function(ancthresh) {
  # Displays tables like in my paper. Works only with two groups.
  
  # Getting back simple objects
  space <- ancthresh$space
  space.sd <- ancthresh$space.sd
  tstat.list <- ancthresh$tstat.list
  
  table <- numeric(0)
  for (i in 1:2) {
  table <- rbind(table,
  				 cbind(space[i, 2], space[i, 1], tstat.list[[i]][1,2], ((space[i, 1] - space[i, 2]) / space[i, 2]) * 100) ,
  				 cbind(space.sd[i,2], space.sd[i,1], NA, NA))
  }
  
  return(table)
}