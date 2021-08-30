# post-hoc-HBEF-BBS.R: code to summarize NIMBLE results from model using 
#                      NEON and BBS data and compute biodiversity metrics and perform
#                      post-hoc trend analysis.
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
load("results/icom-HBEF-BBS-results-450000-iterations-1-chain-2021-08-23.R")
samples.a <- samples
load("results/icom-HBEF-BBS-results-450000-iterations-2-chain-2021-08-23.R")
samples.b <- samples
load("results/icom-HBEF-BBS-results-450000-iterations-3-chain-2021-08-23.R")
samples.c <- samples
samples.list <- mcmc.list(samples.a, samples.b, samples.c)
param.names <- attr(samples.a, 'dimnames')[[2]]
# Indices for latent occurrence values
z.hbef.indx <- which(substr(param.names, 1, 6) %in% c('z.hbef'))
z.bbs.indx <- which(substr(param.names, 1, 5) %in% c('z.bbs'))
psi.hbef.indx <- which(substr(param.names, 1, 8) %in% c('psi.hbef'))
psi.bbs.indx <- which(substr(param.names, 1, 7) %in% c('psi.bbs'))
samples <- samples.list[, -c(z.hbef.indx, z.bbs.indx, 
			     psi.hbef.indx, psi.bbs.indx)]
# Gelman-Rubin diagnostic for convergence check
r.hat.vals <- gelman.diag(samples)
z.hbef.samples <- mcmc.list(mcmc(samples.a[5001:7500, z.hbef.indx]))
z.bbs.samples <- mcmc.list(mcmc(samples.a[5001:7500, z.bbs.indx]))
psi.hbef.samples <- mcmc.list(mcmc(samples.a[5001:7500, psi.hbef.indx]))
psi.bbs.samples <- mcmc.list(mcmc(samples.a[5001:7500, psi.bbs.indx]))
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
# Species richness --------------------
rich.hbef.samples <- apply(z.hbef.wide, c(1, 3, 4), sum)
# Jaccard index -----------------------
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


# BBS
I.bbs <- n_distinct(v.2.df$Species)
J.bbs <- n_distinct(v.2.df$Site)
T.bbs <- n_distinct(v.2.df$Year)
z.bbs.samples <- mcmc(do.call('rbind', z.bbs.samples))
n.iter <- nrow(z.bbs.samples)
z.bbs.wide <- array(z.bbs.samples, dim = c(n.iter, I.bbs, J.bbs, T.bbs))
# Species richness --------------------
rich.bbs.samples <- apply(z.bbs.wide, c(1, 3, 4), sum)
# Jaccard index -----------------------
jaccard.bbs.samples <- array(NA, dim = dim(rich.bbs.samples))
ref.site <- 1
for (a in 1:n.iter) {
  print(paste('Currently on iteration ', a, ' out of ', n.iter, sep = ''))
  for (t in 1:T.bbs) {
    for (j in 1:J.bbs) {
      jaccard.bbs.samples[a, j, t] <- sum(z.hbef.wide[a, , ref.site, t] * z.bbs.wide[a, , j, t]) / 
	      (sum(z.hbef.wide[a, , ref.site, t]) + sum(z.bbs.wide[a, , j, t]) - 
	      sum(z.hbef.wide[a, , ref.site, t] * z.bbs.wide[a, , j, t]))
    } # j
  } # t
} # a 

# Temporal Trend Estimation 
psi.bbs.samples <- mcmc(do.call('rbind', psi.bbs.samples))
psi.bbs.wide <- array(psi.bbs.samples, dim = c(n.iter, I.bbs, J.bbs, T.bbs))
psi.hbef.samples <- mcmc(do.call('rbind', psi.hbef.samples))
psi.hbef.wide <- array(psi.hbef.samples, dim = c(n.iter, I.hbef, J.hbef, T.hbef))
psi.all.wide <- abind(psi.hbef.wide, psi.bbs.wide, along = 3)
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
beta.samples <- array(0, c(I.bbs, n.beta, n.iter))
sigma.sq.samples <- matrix(0, I.bbs, n.iter)

# Priors -------------------------------

# beta: normal prior
mu.beta <- rep(0, n.beta)
Sigma.beta <- diag(1000, n.beta)

# sigma.sq: inverse gamma prior
a.sigma.sq <- .001
b.sigma.sq <- .001

# Need to change year values from the actual large values. 
# Change so the minimum year is 1
years.bbs <- 2010:2018
year.un <- years.bbs - min(years.bbs, na.rm = TRUE) + 1

for (i in 1:I.bbs) {
  print(i)
  # Initial Values --------------------
  beta <- coefficients(lm(apply(psi.avg.sp.year[, i, ], 2, mean) ~ year.un))
  sigma.sq <- 1
  X <- matrix(1, T.bbs, n.beta)
  X[, 2] <- year.un
  XTX <- t(X) %*% X
  for (a in 1:n.iter) {
     # if (a %% 100 == 0) print(a)
     y <- psi.avg.sp.year[a, i, ]
     V <- t(X) %*% chol2inv(chol(diag(sigma.sq, T.bbs))) %*% X + chol2inv(chol(Sigma.beta))
     v <- t(X) %*% chol2inv(chol(diag(sigma.sq, T.bbs))) %*% y + 
	     chol2inv(chol(Sigma.beta)) %*% mu.beta
     V.inv <- chol2inv(chol(V))
     beta <- rmvn(1, V.inv %*% v, V.inv)

     # Sample sigma.sq --------------------
     a.sigma.post <- a.sigma.sq + T.bbs / 2
     b.sigma.post <- b.sigma.sq + sum((y - X%*%beta)^2) / 2
     sigma.sq <- rigamma(1, a.sigma.post, b.sigma.post)

     # Save samples -----------------------
     beta.samples[i, ,a] <- beta
     sigma.sq.samples[i, a] <- sigma.sq
  }
}

save(samples, jaccard.bbs.samples, rich.bbs.samples, 
     jaccard.hbef.samples, rich.hbef.samples, psi.avg.sp.year, 
     beta.samples, sigma.sq.samples, 
     file = 'results/icom-HBEF-BBS-results.R')
