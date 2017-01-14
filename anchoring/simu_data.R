# USEFUL FONCTION
XBar <- function(x, group.index) {
  # XBar: from a matrix x, creates the matrix x.bar.
  # Remark: x and group.index have to be ordered by the levels of group.index. This guarantees that the rows of x.bar correspond to the rows of x. 
  #         Otherwise, the program stops.
  
  # Sanity check
  if (any(group.index != group.index[order(group.index)])) stop("XBar: group.index is not sorted (meaning that x is not sorter by the levels of group.index either")
  if (nrow(x) != length(group.index)) stop("XBar: the number of row of x is not compatible with the length of group.index")
  
  # Constants
  nt.g <- as.vector(table(group.index))  # Matrix that contains the number of observations in each group
  kG <- length(nt.g)  # number of groups
  kL <- ncol(x) # Number of variables in x
  
  # Calculation
  for (g in 1:kG) {
  	if (g == 1) {
  	  x.bar <- cbind(x[1 : nt.g[1], ], matrix(0, nrow = nt.g[1], ncol = (kG - 1) * kL))
  	  row.indic <- nt.g[1]  # Indicates the latest computed row
  	}
  	if (g > 1 & g < kG) {
  	  x.bar <- rbind(x.bar,
  	  				 cbind(matrix(0, nrow = nt.g[g], ncol = (g - 1) * kL), x[(row.indic + 1) : (row.indic + nt.g[g]), ],
  	  				       matrix(0, nrow = nt.g[g], ncol = (kG - g) * kL)))  	  
  	  row.indic <- row.indic + nt.g[g]
  	}
  	if (g == kG) {
  	  x.bar <- rbind(x.bar,
  	  		   cbind(matrix(0, nrow = nt.g[kG], ncol = (kG - 1) * kL), x[(row.indic + 1) : nrow(x), ]))
  	}
  }
  
  # Return
  return(x.bar)
}

SigmaBar <- function(sigma, group.index) {
  # SigmaUBar: calculates a col-vector for wich each observation corresponds to an individual observation (i,t) of sigma.g, depending on the group
  # 		   the respondent belongs to.
  # Remark: group.index has to be sorted (and obviously, the ordering of the rows in x has to be conistent with this index). See XBar for
  #   		more details.
  
  # Sanity checks (desactivited because the two conditions are guaranteed in AncReg AND because it slows the algo down a lot)
#  if (any(group.index != group.index[order(group.index)])) stop("SigmaUBar: group.index is not sorted (meaning that x and index are not sorted by the levels of group.index either")
#  if (nrow(sigma) != length(levels(as.factor(group.index)))) stop("SigmaUBar: the number of levels in group.index is not equal to the length of sigma")

  # Constants
  nt.g <- as.vector(table(group.index))  # Matrix that contains the number of observations in each group
  kG <- length(nt.g)  # number of groups
  
  # Calculation
  sigma.bar <- numeric(0)
  for (g in 1 : kG) {
    sigma.bar <- c(sigma.bar, rep.int(sigma[g], nt.g[g]))
  }
  return(sigma.bar)
}



# SIMULATING DATA
# Parameters not meant to be changed
kG <- 3  # Number of groups
kK <- 3  # Number of statements
kLs <- 2  # Number of variables in x.s
kL <- 2  # Number of variables in x

# Number of periods and groups
kT <- 2
n.g <- c(300, 300, 300)   # Numbers of respondents in each group ; not to be confused with nt.g, number of observations in each group

# Parameters
sigma.eps <- matrix(c(1, 0.75, 0.75), ncol = 1)  # col-vector containing the sigmas for each group
sigma.u <- matrix(c(0.5, 0.75, 0.25), ncol = 1)
beta.s <- matrix(c(0.25, 0.5), ncol = 1)  # parameters for the special variables
beta.bar <- matrix(c(0.5, 0.75, 0.25, 0.5, 0.6, 0.7), ncol = 1)  # each column g is the vector of para beta.bar for group g
tau <- cbind(c(-Inf, -1.5, 1.5, +Inf), c(-Inf, -1.25, 1.25, +Inf), c(-Inf, -1.5, 1.5,+Inf))  # each column g is the vector of thresholds for group g

# Simulating the data
# Indexes
kN <- sum(n.g)  # number of respondents
kNT <- kN * kT # number of observations
index <- 1:kN %x% rep.int(1, kT)  # vector of individual index
group.index <- c(1 %x% rep.int(1, n.g[1] * kT), 2 %x% rep.int(1, n.g[2] * kT), 3 %x% rep.int(1, n.g[3] * kT))  # vector of group index
# Latent data
u <- rnorm(kN) %x% matrix(1, nrow = kT)
x.s <- cbind(rnorm(kN * kT, mean = 0.3, sd = 0.5), rnorm(kN * kT, sd = 0.5))
x <- cbind(rnorm(kN * kT, mean = 0.3, sd = 0.5), rnorm(kN * kT, sd = 0.5))
epsilon <- cbind(rnorm(kN * kT))

x.bar <- XBar(x, group.index)
sigma.u.bar <- SigmaBar(sigma.u, group.index)
sigma.eps.bar <- SigmaBar(sigma.eps, group.index)

y.star <- x.s %*% beta.s + x.bar %*% beta.bar + sigma.u.bar * u + sigma.eps.bar * epsilon

# Statements
y <- 1 * (tau[1,group.index] <= y.star & y.star < tau[2,group.index]) + 2 * (tau[2,group.index] <= y.star & y.star < tau[3,group.index]) + 3 * (tau[3,group.index] <= y.star & y.star < tau[4,group.index])

# Further computations
beta.bar.star <- cbind(c(beta.bar[1:2] / sigma.eps[1], beta.bar[3:4] / sigma.eps[2], beta.bar[5:6] / sigma.eps[3]))
tau.star <- cbind(tau[, 1] / sigma.eps[1], tau[, 2] / sigma.eps[2], tau[, 3] / sigma.eps[3])
tau.star.0 <- tau.star[c(-1, -(kK + 1)), ]
sigma.eps.0 <- sigma.eps[-1, , drop = FALSE]
kappa.eps.0 <- log(sigma.eps.0)
gamma.star.0 <- rbind(tau.star.0[1, ], log(tau.star.0[2, ] - tau.star.0[1, ]))
sigma.u.star <- sigma.u / sigma.eps

# Further computations - reduced-form model
b <- c(1, 2.5 / 3, 1)
b.0 <- b[-1]
Gamma.d.0 <- (tau[2, 1] - tau[cbind(2, c(2,3))] / b.0)
Gamma.bar.s <- c(beta.s, beta.s / b[2], beta.s / b[3])
Gamma.bar <- c(beta.bar[1:2] / b[1], beta.bar[3:4] / b[2], beta.bar[5:6] / b[3])
sigma.v <- sigma.u / b
sigma.eta <- sigma.eps / b

Gamma.star.d.0 <- Gamma.d.0 / sigma.eta[-1]
Gamma.bar.star.s <- c(Gamma.bar.s[1:2] / sigma.eta[1], Gamma.bar.s[3:4] / sigma.eta[2], Gamma.bar.s[5:6] / sigma.eta[3])
Gamma.bar.star <- c(Gamma.bar[1:2] / sigma.eta[1], Gamma.bar[3:4] / sigma.eta[2], Gamma.bar[5:6] / sigma.eta[3])
sigma.v.star <- sigma.v / sigma.eta
gamma.0 <- c(tau[2, 1], log(tau[3, 1] - tau[2, 1])) # gamma for the thresholds
kappa.eta <- log(sigma.eta)
kappa.v.star <- log(sigma.v.star)