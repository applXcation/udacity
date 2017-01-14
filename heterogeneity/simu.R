# Generating data

###########
# CONSTANTS
###########
# Scalars
kN <- 5e+2
heterosk <- TRUE
kK <- 4  # This parameter is not meant to be played around with
kV <- 2  # This parameter is not meant to be played around with

# Coefficients
beta <- rbind(0, 1.5, 2)
gamma <- cbind(c(-1.5, 0.5, 1),c(1, .05, .05),c(.6, .05, .05)) # gamma = (gamma.1, gamma.2, gamma.3)
theta <- cbind(3, 4)
kappa <- rbind(0, 0.5, -0.5)
sigma.v1 <- 1
sigma.v2 <- 1.3
sigma.v <- cbind(sigma.v1, sigma.v2)
sigma.tilde.v <- log(1 / sigma.v)

###########
# VARIABLES
###########
# Generating the xs
x <- cbind(replicate(kN, 1), replicate(2, rnorm(kN, mean = 1, sd = 1.5), 2))
names(x) <- c("const", "x1", "x2")

# Generating the sigma.s
sigma <- heterosk * exp(x %*% kappa) + (1-heterosk) * replicate(kN, 1)

# Generating the tau
tau.0 <- matrix(-Inf, nrow = kN) ; tau.4 <- matrix(+Inf, nrow = kN)
tau.1 <- x %*% gamma[, 1, drop = FALSE]
tau.2 <- tau.1 + exp(x %*% gamma[, 2, drop = FALSE])
tau.3 <- tau.2 + exp(x %*% gamma[, 3, drop = FALSE])
tau <- cbind(tau.0, tau.1, tau.2, tau.3, tau.4)

# Generating the self-assessment variable
epsilon.s <- rnorm(kN)
y.star.s <- x %*% beta + sigma * epsilon.s
y.s <- matrix(0, nrow = kN)
for (k in 1:kK) {
  y.s <- y.s + k * (tau[, k] <= y.star.s & y.star.s < tau[, k+1])
}

# Generating the vignette-asessment variables
epsilon.v1 <- rnorm(kN) ; epsilon.v2 <- rnorm(kN)
y.star.v1 <- theta[1] + sigma.v1 * epsilon.v1
y.star.v2 <- theta[2] + sigma.v2 * epsilon.v2
y.v1 <- matrix(0, nrow = kN)
y.v2 <- matrix(0, nrow = kN)
for (k in 1:kK) {
  y.v1 <- y.v1 + k * (tau[, k] <= y.star.v1 & y.star.v1 < tau[, k+1])
  y.v2 <- y.v2 + k * (tau[, k] <= y.star.v2 & y.star.v2 < tau[, k+1])
}

y.v <- cbind(y.v1, y.v2)

# Generating a grouping variable creating two groups in the population ; the first group has higher meann lower variance
group <- 1 * (x[, 2] < quantile(x[, 2], probs = .30) | x[, 2] > quantile(x[, 2], probs = .95))
group2 <- rbinom(kN,1,0.5)

# Validating the generated data
table(y.s)
table(y.v1)
table(y.v2)
tapply(y.star.s, group, mean)
tapply(y.star.s, group, var)
