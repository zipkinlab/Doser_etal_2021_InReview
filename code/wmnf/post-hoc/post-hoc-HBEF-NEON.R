# post-hoc-HBEF-NEON.R: code to summarize NIMBLE results from model using 
#                       HBEF and NEON data and compute biodiversity metrics and perform
#                       post-hoc trend analysis.
# Author: Jeffrey W. Doser (doserjef@msu.edu)
# Citation: 
rm(list = ls())
library(coda)
library(dplyr)
# For binding arrays
library(abind)

# Read in data ------------------------------------------------------------
load("data/nimble-data.R")

# Read in results ---------------------------------------------------------
load("results/icom-HBEF-NEON-results-450000-iterations-1-chain-2021-08-22.R")
samples.a <- samples
load("results/icom-HBEF-NEON-results-450000-iterations-2-chain-2021-08-22.R")
samples.b <- samples
load("results/icom-HBEF-NEON-results-450000-iterations-3-chain-2021-08-22.R")
samples.c <- samples
samples.list <- mcmc.list(samples.a, samples.b, samples.c)
param.names <- attr(samples.a, 'dimnames')[[2]]
# Indices for latent occurrence values
z.hbef.indx <- which(substr(param.names, 1, 6) %in% c('z.hbef'))
z.neon.indx <- which(substr(param.names, 1, 6) %in% c('z.neon'))
psi.hbef.indx <- which(substr(param.names, 1, 8) %in% c('psi.hbef'))
psi.neon.indx <- which(substr(param.names, 1, 8) %in% c('psi.neon'))
# Get the MCMC samples corresponding to the actual parameters (not latent states). 
# This makes computation easier. 
samples <- samples.list[, -c(z.hbef.indx, z.neon.indx, 
			     psi.hbef.indx, psi.neon.indx)]
# Gelman-Rubin diagnostic for convergence check
r.hat.vals <- gelman.diag(samples)
# Subset for post-hoc analysis
z.hbef.samples <- mcmc.list(mcmc(samples.a[5001:7500, z.hbef.indx]))
z.neon.samples <- mcmc.list(mcmc(samples.a[5001:7500, z.neon.indx]))
psi.hbef.samples <- mcmc.list(mcmc(samples.a[5001:7500, psi.hbef.indx]))
psi.neon.samples <- mcmc.list(mcmc(samples.a[5001:7500, psi.neon.indx]))
sp <- c('AMRE', 'BAWW', 'BHVI', 'BLBW', 'BLPW', 'BTBW',
        'BTNW', 'CAWA', 'MAWA', 'NAWA', 'OVEN', 'REVI')
sp.names <- c('American Redstart', 'Black-and-white Warbler', 'Blue-headed Vireo', 
	      'Blackburnian Warbler', 'Blackpoll Warbler', 'Black-throated Blue Warbler', 
	      'Black-throated Green Warbler', 'Canada Warbler', 'Magnolia Warbler', 
	      'Nashville Warbler', 'Ovenbird', 'Red-eyed Vireo')


# Estimate species richness and Jaccard index -----------------------------
# Ordered by species, site, year
# HBEF
I.hbef <- n_distinct(y.df$Species)
J.hbef <- n_distinct(y.df$Site)
T.hbef <- n_distinct(y.df$Year)
z.hbef.samples <- mcmc(do.call('rbind', z.hbef.samples))
n.iter <- nrow(z.hbef.samples)
z.hbef.wide <- array(z.hbef.samples, dim = c(n.iter, I.hbef, J.hbef, T.hbef))
# Species richness
rich.hbef.samples <- apply(z.hbef.wide, c(1, 3, 4), sum)
# Compute Jaccard index
jaccard.hbef.samples <- array(NA, dim = dim(rich.hbef.samples))
ref.site <- 1
for (a in 1:n.iter) {
  print(paste('Currently on iteration ', a, ' out of ', n.iter, sep = ''))
  for (j in 1:J.hbef) {
    for (t in 1:T.hbef) {
      jaccard.hbef.samples[a, j, t] <- sum(z.hbef.wide[a, , ref.site, t] * z.hbef.wide[a, , j, t]) / 
	      (sum(z.hbef.wide[a, , ref.site, t]) + sum(z.hbef.wide[a, , j, t]) - 
	      sum(z.hbef.wide[a, , ref.site, t] * z.hbef.wide[a, , j, t]))
    } # t
  } # j
} # a 


# NEON
I.neon <- n_distinct(v.1.df$Species)
J.neon <- n_distinct(v.1.df$Site)
T.neon <- n_distinct(v.1.df$Year)
z.neon.samples <- mcmc(do.call('rbind', z.neon.samples))
n.iter <- nrow(z.neon.samples)
# Make this the same size as the HBEF data, just fill in the correct years
z.neon.wide <- array(z.neon.samples, dim = c(n.iter, I.neon, J.neon, T.neon))
# Species richness --------------------
rich.neon.samples <- apply(z.neon.wide, c(1, 3, 4), sum)
# Jaccard index -----------------------
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

# Temporal Trend Estimation 
neon.year.indx <- unique(v.1.df$Year.all)
neon.sp.indx <- unique(v.1.df$Species)
psi.hbef.samples <- mcmc(do.call('rbind', psi.hbef.samples))
psi.hbef.wide <- array(psi.hbef.samples, dim = c(n.iter, I.hbef, J.hbef, T.hbef))
psi.neon.samples <- mcmc(do.call('rbind', psi.neon.samples))
psi.neon.wide <- array(NA, dim = c(n.iter, I.hbef, J.neon, T.hbef))
psi.neon.wide[, neon.sp.indx, , neon.year.indx] <- array(psi.neon.samples, 
					                 dim = c(n.iter, I.neon, J.neon, T.neon))
psi.all.wide <- abind(psi.hbef.wide, psi.neon.wide, along = 3)
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
beta.samples <- array(0, c(I.hbef, n.beta, n.iter))
sigma.sq.samples <- matrix(0, I.hbef, n.iter)

# Priors -------------------------------

# beta: normal prior
mu.beta <- rep(0, n.beta)
Sigma.beta <- diag(1000, n.beta)

# sigma.sq: inverse gamma prior
a.sigma.sq <- .001
b.sigma.sq <- .001

# Need to change year values from the actual large values. 
# Change so the minimum year is 1
years.hbef <- 2010:2018
year.un <- years.hbef - min(years.hbef, na.rm = TRUE) + 1

for (i in 1:I.hbef) {
  print(i)
  # Initial Values --------------------
  beta <- coefficients(lm(apply(psi.avg.sp.year[, i, ], 2, mean) ~ year.un))
  sigma.sq <- 1
  X <- matrix(1, T.hbef, n.beta)
  X[, 2] <- year.un
  XTX <- t(X) %*% X
  for (a in 1:n.iter) {
     # if (a %% 100 == 0) print(a)
     y <- psi.avg.sp.year[a, i, ]
     V <- t(X) %*% chol2inv(chol(diag(sigma.sq, T.hbef))) %*% X + chol2inv(chol(Sigma.beta))
     v <- t(X) %*% chol2inv(chol(diag(sigma.sq, T.hbef))) %*% y + 
	     chol2inv(chol(Sigma.beta)) %*% mu.beta
     V.inv <- chol2inv(chol(V))
     beta <- rmvn(1, V.inv %*% v, V.inv)

     # Sample sigma.sq --------------------
     a.sigma.post <- a.sigma.sq + T.hbef / 2
     b.sigma.post <- b.sigma.sq + sum((y - X%*%beta)^2) / 2
     sigma.sq <- rigamma(1, a.sigma.post, b.sigma.post)

     # Save samples -----------------------
     beta.samples[i, ,a] <- beta
     sigma.sq.samples[i, a] <- sigma.sq
  }
}

save(samples, jaccard.neon.samples, rich.neon.samples, 
     jaccard.hbef.samples, rich.hbef.samples, psi.avg.sp.year, 
     beta.samples, sigma.sq.samples,
     file = 'results/icom-HBEF-NEON-results.R')
