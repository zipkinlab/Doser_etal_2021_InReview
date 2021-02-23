rm(list = ls())
library(jagsUI)
library(coda)
library(dplyr)
library(ggplot2)

# Functions ---------------------------------------------------------------
# Logit transformation
logit <- function(theta, a = 0, b = 1) {
  log((theta-a)/(b-theta))
}

# Inverse logit transformation
logit.inv <- function(z, a = 0, b = 1) {
  b-(b-a)/(1+exp(z))
}


# Read in results ---------------------------------------------------------
# OVEN
load("results/isdm-small-OVEN-results-50000-iterations-2021-02-04.R")
oven.out <- out
# OVEN with entire 3 data sets and range
load("results/isdm-OVEN-results-50000-iterations-2021-02-06.R")
oven.full.out <- out
# OVEN with 3 data sets in small range
load("results/isdm-small-eb-OVEN-results-50000-iterations-2021-02-07.R")
oven.eb.out <- out
# BAWW
load("results/isdm-small-BAWW-results-50000-iterations-2021-02-13.R")
baww.out <- out
# REVI
load("results/isdm-small-REVI-results-50000-iterations-2021-02-13.R")
revi.out <- out
# BTBW
load("results/isdm-small-BTBW-results-50000-iterations-2021-02-04.R")
btbw.out <- out
# BTNW
load("results/isdm-small-BTNW-results-50000-iterations-2021-02-13.R")
btnw.out <- out
# AMRE
load("results/isdm-small-AMRE-results-50000-iterations-2021-02-04.R")
amre.out <- out
# NAWA
load("results/isdm-small-NAWA-results-50000-iterations-2021-02-04.R")
nawa.out <- out
# BHVI
load("results/isdm-small-BHVI-results-50000-iterations-2021-02-05.R")
bhvi.out <- out
# BLPW
load("results/isdm-small-BLPW-results-50000-iterations-2021-02-05.R")
blpw.out <- out
# BLBW
load("results/isdm-small-BLBW-results-50000-iterations-2021-02-05.R")
blbw.out <- out
# MAWA
load("results/isdm-small-MAWA-results-50000-iterations-2021-02-05.R")
mawa.out <- out
# CAWA
load("results/isdm-small-CAWA-results-50000-iterations-2021-02-13.R")
cawa.out <- out
# Species are: (1) American Redstart, (2) Black and white warbler, 
# (3) Blue-headed vireo, (4) Blackburnian warbler, (5) Blackpoll warbler, 
# (6) Black-throated blue, (7) Black throated green, (8) Canadian warbler, 
# (9) Magnolia warbler, (10) Nashville warbler, (11) Ovenbird, (12) red-eyed vireo

# Elevation Preference

# Low: BTBW, OVEN, REVI, BTNW, BAWW, CAWA, BLBW, AMRE, BHVI (although seems to range)
# High: BLPW, MAWA, NAWA

n.iter <- out$mcmc.info$n.samples
K <- 12
# Temperature effect on colonization
beta.gamma.1.samples <- matrix(NA, n.iter, K)
beta.gamma.1.samples[, 1] <- amre.out$sims.list$beta.gamma.1
beta.gamma.1.samples[, 2] <- baww.out$sims.list$beta.gamma.1
beta.gamma.1.samples[, 3] <- bhvi.out$sims.list$beta.gamma.1
beta.gamma.1.samples[, 4] <- blbw.out$sims.list$beta.gamma.1
beta.gamma.1.samples[, 5] <- blpw.out$sims.list$beta.gamma.1
beta.gamma.1.samples[, 6] <- btbw.out$sims.list$beta.gamma.1
beta.gamma.1.samples[, 7] <- btnw.out$sims.list$beta.gamma.1
beta.gamma.1.samples[, 8] <- cawa.out$sims.list$beta.gamma.1
beta.gamma.1.samples[, 9] <- mawa.out$sims.list$beta.gamma.1
beta.gamma.1.samples[, 10] <- nawa.out$sims.list$beta.gamma.1
beta.gamma.1.samples[, 11] <- oven.out$sims.list$beta.gamma.1
beta.gamma.1.samples[, 12] <- revi.out$sims.list$beta.gamma.1
sp <- c('AMRE', 'BAWW', 'BHVI', 'BLBW', 'BLPW', 'BTBW', 
	'BTNW', 'CAWA', 'MAWA', 'NAWA', 'OVEN', 'REVI')
colnames(beta.gamma.1.samples) <- sp
summary(mcmc(beta.gamma.1.samples))

# Precipitation effect on colonization
beta.gamma.2.samples <- matrix(NA, n.iter, K)
beta.gamma.2.samples[, 1] <- amre.out$sims.list$beta.gamma.2
beta.gamma.2.samples[, 2] <- baww.out$sims.list$beta.gamma.2
beta.gamma.2.samples[, 3] <- bhvi.out$sims.list$beta.gamma.2
beta.gamma.2.samples[, 4] <- blbw.out$sims.list$beta.gamma.2
beta.gamma.2.samples[, 5] <- blpw.out$sims.list$beta.gamma.2
beta.gamma.2.samples[, 6] <- btbw.out$sims.list$beta.gamma.2
beta.gamma.2.samples[, 7] <- btnw.out$sims.list$beta.gamma.2
beta.gamma.2.samples[, 8] <- cawa.out$sims.list$beta.gamma.2
beta.gamma.2.samples[, 9] <- mawa.out$sims.list$beta.gamma.2
beta.gamma.2.samples[, 10] <- nawa.out$sims.list$beta.gamma.2
beta.gamma.2.samples[, 11] <- oven.out$sims.list$beta.gamma.2
beta.gamma.2.samples[, 12] <- revi.out$sims.list$beta.gamma.2
sp <- c('AMRE', 'BAWW', 'BHVI', 'BLBW', 'BLPW', 'BTBW', 
	'BTNW', 'CAWA', 'MAWA', 'NAWA', 'OVEN', 'REVI')
colnames(beta.gamma.2.samples) <- sp
summary(mcmc(beta.gamma.2.samples))

# Temperature effect on persistence
beta.phi.1.samples <- matrix(NA, n.iter, K)
beta.phi.1.samples[, 1] <- amre.out$sims.list$beta.phi.1
beta.phi.1.samples[, 2] <- baww.out$sims.list$beta.phi.1
beta.phi.1.samples[, 3] <- bhvi.out$sims.list$beta.phi.1
beta.phi.1.samples[, 4] <- blbw.out$sims.list$beta.phi.1
beta.phi.1.samples[, 5] <- blpw.out$sims.list$beta.phi.1
beta.phi.1.samples[, 6] <- btbw.out$sims.list$beta.phi.1
beta.phi.1.samples[, 7] <- btnw.out$sims.list$beta.phi.1
beta.phi.1.samples[, 8] <- cawa.out$sims.list$beta.phi.1
beta.phi.1.samples[, 9] <- mawa.out$sims.list$beta.phi.1
beta.phi.1.samples[, 10] <- nawa.out$sims.list$beta.phi.1
beta.phi.1.samples[, 11] <- oven.out$sims.list$beta.phi.1
beta.phi.1.samples[, 12] <- revi.out$sims.list$beta.phi.1
sp <- c('AMRE', 'BAWW', 'BHVI', 'BLBW', 'BLPW', 'BTBW', 
	'BTNW', 'CAWA', 'MAWA', 'NAWA', 'OVEN', 'REVI')
colnames(beta.phi.1.samples) <- sp
summary(mcmc(beta.phi.1.samples))

# Precipitation effect on persistence
summary(mcmc(beta.phi.2.samples))
beta.phi.2.samples <- matrix(NA, n.iter, K)
beta.phi.2.samples[, 1] <- amre.out$sims.list$beta.phi.2
beta.phi.2.samples[, 2] <- baww.out$sims.list$beta.phi.2
beta.phi.2.samples[, 3] <- bhvi.out$sims.list$beta.phi.2
beta.phi.2.samples[, 4] <- blbw.out$sims.list$beta.phi.2
beta.phi.2.samples[, 5] <- blpw.out$sims.list$beta.phi.2
beta.phi.2.samples[, 6] <- btbw.out$sims.list$beta.phi.2
beta.phi.2.samples[, 7] <- btnw.out$sims.list$beta.phi.2
beta.phi.2.samples[, 8] <- cawa.out$sims.list$beta.phi.2
beta.phi.2.samples[, 9] <- mawa.out$sims.list$beta.phi.2
beta.phi.2.samples[, 10] <- nawa.out$sims.list$beta.phi.2
beta.phi.2.samples[, 11] <- oven.out$sims.list$beta.phi.2
beta.phi.2.samples[, 12] <- revi.out$sims.list$beta.phi.2
sp <- c('AMRE', 'BAWW', 'BHVI', 'BLBW', 'BLPW', 'BTBW', 
	'BTNW', 'CAWA', 'MAWA', 'NAWA', 'OVEN', 'REVI')
colnames(beta.phi.2.samples) <- sp
summary(mcmc(beta.phi.2.samples))


# Intercept on initial occupancy
summary(mcmc(int.psi.samples))
int.psi.samples <- matrix(NA, n.iter, K)
int.psi.samples[, 1] <- amre.out$sims.list$int.psi
int.psi.samples[, 2] <- baww.out$sims.list$int.psi
int.psi.samples[, 3] <- bhvi.out$sims.list$int.psi
int.psi.samples[, 4] <- blbw.out$sims.list$int.psi
int.psi.samples[, 5] <- blpw.out$sims.list$int.psi
int.psi.samples[, 6] <- btbw.out$sims.list$int.psi
int.psi.samples[, 7] <- btnw.out$sims.list$int.psi
int.psi.samples[, 8] <- cawa.out$sims.list$int.psi
int.psi.samples[, 9] <- mawa.out$sims.list$int.psi
int.psi.samples[, 10] <- nawa.out$sims.list$int.psi
int.psi.samples[, 11] <- oven.out$sims.list$int.psi
int.psi.samples[, 12] <- revi.out$sims.list$int.psi
sp <- c('AMRE', 'BAWW', 'BHVI', 'BLBW', 'BLPW', 'BTBW', 
	'BTNW', 'CAWA', 'MAWA', 'NAWA', 'OVEN', 'REVI')
colnames(int.psi.samples) <- sp
summary(mcmc(int.psi.samples))
