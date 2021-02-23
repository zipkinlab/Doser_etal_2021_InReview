rm(list = ls())
library(dplyr)
library(ggplot2)
library(jagsUI)
library(coda)
set.seed(1010)
load("../../data/bird-500m-3km-model-data.R")

# Data prep for analysis --------------------------------------------------
K <- dim(y.ebird)[2]
R <- nrow(pixels)
R.hb <- dim(hb.dat)[1]
J.hb <- dim(hb.dat)[2]
R.eb <- dim(y.ebird)[3]
J.eb <- dim(y.ebird)[1]
R.neon <- dim(y.neon)[1]
J.neon <- dim(y.neon)[2]
n.years <- dim(y.ebird)[4]
years <- as.numeric(attr(hb.day, 'dimnames')[[3]])
# Standardize covariates abundance covariates
years <- c(scale(years))
neon.years <- 6:9
pixels$forest <- c(scale(pixels$forest))
pixels$elevation <- c(scale(pixels$elevation))
ppt.vals <- scale(ppt.vals)
tmean.vals <- scale(tmean.vals)
pixels <- cbind(pixels, ppt.vals, tmean.vals)
names(pixels) <- c('xmin', 'ymin', 'xmax', 'ymax', 'pixel.num', 'cell.num', 
		   'elevation', 'forest', 'centroid.x', 'centroid.y', 
		   'ppt-1', 'ppt-2', 'ppt-3', 'ppt-4', 'ppt-5', 'ppt-6', 
		   'ppt-7', 'ppt-8', 'ppt-9', 'tmean-1', 'tmean-2', 'tmean-3', 
		   'tmean-4', 'tmean-5', 'tmean-6', 'tmean-7', 'tmean-8', 
		   'tmean-9')
# Standardize hubbard brook covariates
hb.day <- array(scale(hb.day), dim = c(R.hb, J.hb, n.years))
hb.tod <- array(scale(hb.tod), dim = c(R.hb, J.hb, n.years))
# Standardize eBird covariates
obsv <- array(scale(obsv), dim = c(J.eb, R.eb, n.years))
dist.trav <- array(scale(dist.trav), dim = dim(obsv))
day.eb <- array(scale(day.eb), dim = dim(obsv))
time.eb <- array(scale(time.eb), dim = dim(obsv))
length.eb <- array(scale(length.eb), dim = dim(obsv))
# Standardize neon data
n.years.neon <- dim(y.neon)[3]
years.neon <- years[6:9]
neon.t.ind <- which(years %in% years.neon)
neon.day <- array(scale(day.neon), dim = c(R.neon, n.years.neon))
neon.hour <- array(scale(hour.neon), dim = c(R.neon, n.years.neon))
# NOTE: this is the order the abundance is estimated in. 
pixels <- pixels %>% arrange(cell.num)
cells <- cells %>% arrange(cell.num)
pixel.hb <- sapply(hb.covs$pixel, FUN = function(a) {which(pixels$pixel.num == a)})
pixel.neon <- sapply(neon.covs$pixel, FUN = function(a) {which(pixels$pixel.num == a)})
# Get low and high vectors for the eBird data for summing across the Poisson process
low <- rep(NA, R.eb)
high <- rep(NA, R.eb)
indx <- 1
for (i in 1:R.eb) {
  curr.cell <- cells$cell.num[i]
  tmp <- sum(pixels$cell.num == curr.cell)
  low[i] <- indx
  high[i] <- indx + tmp - 1
  indx <- indx + tmp 
}
y.ebird <- aperm(y.ebird, c(1, 3, 2, 4))
y.neon <- aperm(y.neon, c(1, 2, 4, 3))
hb.dat <- aperm(hb.dat, c(1, 2, 4, 3))

# Convert all data to presence/absence
eta <- apply(y.neon, c(1, 3, 4), sum, na.rm = TRUE)
x <- ifelse(eta > 0, 1, eta)
C <- ifelse(hb.dat > 0, 1, hb.dat)

# Subset to smaller spatial region for testing
# This really doesn't even speed things up.
hb.indices <- unique(pixel.hb)
neon.indices <- unique(pixel.neon)
curr.cells <- unique(pixels$cell.num[c(hb.indices, neon.indices)])
unique.cells <- unique(pixels$cell.num)
cell.indices <- which(unique.cells %in% curr.cells)
pixels <- pixels %>% filter(cell.num %in% unique(curr.cells))
R <- nrow(pixels)
y.ebird <- y.ebird[, cell.indices, , ]
R.eb <- dim(y.ebird)[2]
obsv <- obsv[, cell.indices, ]
dist.trav <- dist.trav[, cell.indices, ]
day.eb <- day.eb[, cell.indices, ]
time.eb <- time.eb[, cell.indices, ]
length.eb <- length.eb[, cell.indices, ]
pixel.hb <- sapply(hb.covs$pixel, FUN = function(a) {which(pixels$pixel.num == a)})
pixel.neon <- sapply(neon.covs$pixel, FUN = function(a) {which(pixels$pixel.num == a)})
cells <- cells[cell.indices, ]
low <- rep(NA, R.eb)
high <- rep(NA, R.eb)
indx <- 1
for (i in 1:R.eb) {
  curr.cell <- cells$cell.num[i]
  tmp <- sum(pixels$cell.num == curr.cell)
  low[i] <- indx
  high[i] <- indx + tmp - 1
  indx <- indx + tmp
}


# Bundle data -------------------------------------------------------------
bugs.data <- list(K = K, n.years = n.years, n.years.neon = n.years.neon,
		  R = R, ELEV = pixels$elevation, R.hb = R.hb, R.eb = R.eb,
		  J.hb = J.hb, DAY.hb = hb.day, TOD.hb = hb.tod, J.eb = J.eb,
		  OBSV = obsv, DIST = dist.trav, DAY.eb = day.eb, 
		  TIME.eb = time.eb, LENGTH = length.eb, R.neon = R.neon, 
		  DAY.neon = neon.day, HOUR.neon = neon.hour, 
		  C = C, pixel.hb = pixel.hb, 
		  pixel.neon = pixel.neon, y = y.ebird, x = x, 
		  low = low, high = high, 
		  TMEAN = as.matrix(pixels[, 20:28]), 
		  PPT = as.matrix(pixels[, 11:19]), neon.years = neon.years)
		   
# Initial values ----------------------------------------------------------
z.init <- array(1, dim = c(R, K, n.years))
inits <- function() {
  list(
    z = z.init
  )
}

# Parameters monitored ----------------------------------------------------
params <- c('int.psi.mean', 'beta.psi.1.mean', 'beta.psi.2.mean', 
	    'int.phi.mean', 'beta.phi.1.mean', 'beta.phi.2.mean', 
	    'beta.phi.3.mean', 'beta.phi.4.mean', 'beta.phi.5.mean', 
	    'int.gamma.mean', 'beta.gamma.1.mean', 'beta.gamma.2.mean', 
	    'beta.gamma.3.mean', 'beta.gamma.4.mean', 'beta.gamma.5.mean',
	    'int.alpha.eb.mean', 'alpha.eb.1.mean', 'alpha.eb.2.mean', 
	    'alpha.eb.3.mean', 'alpha.eb.4.mean', 'alpha.eb.5.mean', 
	    'alpha.eb.6.mean', 'int.alpha.hb.mean', 'alpha.hb.1.mean', 
	    'alpha.hb.2.mean', 'alpha.hb.3.mean', 'int.alpha.neon.mean', 
	    'alpha.neon.1.mean', 'alpha.neon.2.mean', 'alpha.neon.3.mean', 
	    'int.psi', 'beta.psi.1', 'beta.psi.2', 'int.phi', 'beta.phi.1', 
	    'beta.phi.2', 'beta.phi.3', 'beta.phi.4', 'beta.phi.5', 
	    'int.gamma', 'beta.gamma.1', 'beta.gamma.2', 'beta.gamma.3', 
	    'beta.gamma.4', 'beta.gamma.5')
# MCMC settings -----------------------------------------------------------
n.iter <- 50000
n.thin <- 20
n.burn <- 30000
n.chain <- 3

out <- jags(bugs.data, inits, params, 'icm-birds-jags.txt',
            n.iter = n.iter, n.thin = n.thin, 
            n.burn = n.burn, n.chain = n.chain, parallel = TRUE)

date <- Sys.Date()
file.name <- paste('../../results/icm-3km-small-bird-results-', n.iter, '-iterations-', date, '.R', sep = '')
save(out, bugs.data, file = file.name)
