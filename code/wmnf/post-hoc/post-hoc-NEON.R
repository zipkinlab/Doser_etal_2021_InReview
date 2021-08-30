# post-hoc-NEON.R: code to summarize NIMBLE results from model using only 
#                  NEON data and compute biodiversity metrics and perform
#                  post-hoc trend analysis.
# Author: Jeffrey W. Doser (doserjef@msu.edu)
# Citation: 

rm(list = ls())
library(coda)
library(dplyr)

# Read in data ------------------------------------------------------------
load("data/nimble-data.R")

# Read in results ---------------------------------------------------------
load("results/icom-NEON-results-450000-iterations-1-chain-2021-08-18.R")
samples.a <- samples
load("results/icom-NEON-results-450000-iterations-2-chain-2021-08-18.R")
samples.b <- samples
load("results/icom-NEON-results-450000-iterations-3-chain-2021-08-18.R")
samples.c <- samples
samples.list <- mcmc.list(samples.a, samples.b, samples.c)
param.names <- attr(samples.a, 'dimnames')[[2]]
# Indices for the latent occurrence values
z.neon.indx <- which(substr(param.names, 1, 6) %in% c('z.neon'))
psi.neon.indx <- which(substr(param.names, 1, 3) %in% c('psi'))
# Get the MCMC samples corresponding to the actual parameters (not latent states). 
# This makes computation easier. 
samples <- samples.list[, -c(z.neon.indx, psi.neon.indx)]
# Gelman-Rubin diagnostic for convergence check. 
r.hat.vals <- gelman.diag(samples)
# Subset for post-hoc analysis. 
z.neon.samples <- mcmc.list(mcmc(samples.a[5001:7500, z.neon.indx]))
psi.neon.samples <- mcmc.list(mcmc(samples.a[5001:7500, psi.neon.indx]))
sp <- c('AMRE', 'BAWW', 'BHVI', 'BLBW', 'BLPW', 'BTBW',
        'BTNW', 'CAWA', 'MAWA', 'NAWA', 'OVEN', 'REVI')
sp.names <- c('American Redstart', 'Black-and-white Warbler', 'Blue-headed Vireo', 
	      'Blackburnian Warbler', 'Blackpoll Warbler', 'Black-throated Blue Warbler', 
	      'Black-throated Green Warbler', 'Canada Warbler', 'Magnolia Warbler', 
	      'Nashville Warbler', 'Ovenbird', 'Red-eyed Vireo')

# Estimate species richness and Jaccard index -----------------------------
# NEON
I.neon <- n_distinct(v.1.df$Species)
J.neon <- n_distinct(v.1.df$Site)
T.neon <- n_distinct(v.1.df$Year)
z.neon.samples <- mcmc(do.call('rbind', z.neon.samples))
n.iter <- nrow(z.neon.samples)
z.neon.wide <- array(z.neon.samples, dim = c(n.iter, I.neon, J.neon, T.neon))
# Species richness --------------------
rich.neon.samples <- apply(z.neon.wide, c(1, 3, 4), sum)
# Jaccard Index -----------------------
jaccard.neon.samples <- array(NA, dim = dim(rich.neon.samples))
ref.site <- 1
for (a in 1:n.iter) {
  print(paste('Currently on iteration ', a, ' out of ', n.iter, sep = ''))
  for (j in 1:J.neon) {
    for (t in 1:T.neon) {
      jaccard.neon.samples[a, j, t] <- sum(z.neon.wide[a, , ref.site, t] * z.neon.wide[a, , j, t]) / 
	      (sum(z.neon.wide[a, , ref.site, t]) + sum(z.neon.wide[a, , j, t]) - 
	      sum(z.neon.wide[a, , ref.site, t] * z.neon.wide[a, , j, t]))
    } # t
  } # j
} # a 

# Post-hoc linear regression as a trend estimate --------------------------
psi.neon.samples <- mcmc(do.call('rbind', psi.neon.samples))
psi.all.wide <- array(psi.neon.samples, dim = c(n.iter, I.neon, J.neon, T.neon))
# Average annual species-specific occupancy probabilities
psi.avg.sp.year <- apply(psi.all.wide, c(1, 2, 4), mean, na.rm = TRUE)

# Gibbs Sampler -----------------------
rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p)))){stop("Dimension problem!")}
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
}

# Inverse Gamma random number generator
rigamma <- function(n, a, b){
  1/rgamma(n = n, shape = a, rate = b)
}


# Set up vectors and arrays -----------
n.beta <- 2
beta.samples <- array(0, c(I.neon, n.beta, n.iter))
sigma.sq.samples <- matrix(0, I.neon, n.iter)

# Priors -------------------------------

# beta: normal prior
mu.beta <- rep(0, n.beta)
Sigma.beta <- diag(1000, n.beta)

# sigma.sq: inverse gamma prior
a.sigma.sq <- .001
b.sigma.sq <- .001

# Need to change year values from the actual large values. 
# Change so the minimum year is 1
years.neon <- 2015:2018
year.un <- years.neon - min(years.neon, na.rm = TRUE) + 1

for (i in 1:I.neon) {
  print(i)
  # Initial Values --------------------
  beta <- coefficients(lm(apply(psi.avg.sp.year[, i, ], 2, mean) ~ year.un))
  sigma.sq <- 1
  X <- matrix(1, T.neon, n.beta)
  X[, 2] <- year.un
  XTX <- t(X) %*% X
  for (a in 1:n.iter) {
     # if (a %% 100 == 0) print(a)
     y <- psi.avg.sp.year[a, i, ]
     V <- t(X) %*% chol2inv(chol(diag(sigma.sq, T.neon))) %*% X + chol2inv(chol(Sigma.beta))
     v <- t(X) %*% chol2inv(chol(diag(sigma.sq, T.neon))) %*% y + 
	     chol2inv(chol(Sigma.beta)) %*% mu.beta
     V.inv <- chol2inv(chol(V))
     beta <- rmvn(1, V.inv %*% v, V.inv)

     # Sample sigma.sq --------------------
     a.sigma.post <- a.sigma.sq + T.neon / 2
     b.sigma.post <- b.sigma.sq + sum((y - X%*%beta)^2) / 2
     sigma.sq <- rigamma(1, a.sigma.post, b.sigma.post)

     # Save samples -----------------------
     beta.samples[i, ,a] <- beta
     sigma.sq.samples[i, a] <- sigma.sq
  }
}

save(samples, jaccard.neon.samples, rich.neon.samples, 
     psi.avg.sp.year, beta.samples, sigma.sq.samples, 
     file = 'results/icom-NEON-results.R')
