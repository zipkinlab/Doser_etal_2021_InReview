# cross-val.R: code to obtain predictive results from cross-validation
#              for the White Mountain National Forest case study using
#              single species integrated distribution models. 
# Author: Jeffrey W. Doser (doserjef@msu.edu)
# Citation: 

rm(list = ls())
library(coda)
library(tidyverse)

# Functions ---------------------------------------------------------------
logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}
# This function computes log pointwise-posterior predictive densities that 
# are used to assess predictive performance
source("code/cross-val-idm/lpd-sampler.R")

# Current hold-out set ----------------------------------------------------
# Takes value 1 or 2 depending on which set is the current hold out set
# Provided as user input when running script from terminal. 
curr.set <- as.numeric(commandArgs(trailingOnly = TRUE))
# curr.set <- 1
pred.set <- ifelse(curr.set == 1, 2, 1)

# Load in Data ------------------------------------------------------------
# AMRE 
samples.amre <- list()
load("results/cross-val/idm-results-250000-iterations-AMRE-species-1-holdout-2021-08-04.R")
samples.amre[[1]] <- samples
load("results/cross-val/idm-results-250000-iterations-AMRE-species-2-holdout-2021-08-04.R")
samples.amre[[2]] <- samples
# BAWW 
samples.baww <- list()
load("results/cross-val/idm-results-250000-iterations-BAWW-species-1-holdout-2021-08-04.R")
samples.baww[[1]] <- samples
load("results/cross-val/idm-results-250000-iterations-BAWW-species-2-holdout-2021-08-04.R")
samples.baww[[2]] <- samples
# BHVI 
samples.bhvi <- list()
load("results/cross-val/idm-results-250000-iterations-BHVI-species-1-holdout-2021-08-04.R")
samples.bhvi[[1]] <- samples
load("results/cross-val/idm-results-250000-iterations-BHVI-species-2-holdout-2021-08-04.R")
samples.bhvi[[2]] <- samples
# BLBW 
samples.blbw <- list()
load("results/cross-val/idm-results-250000-iterations-BLBW-species-1-holdout-2021-08-04.R")
samples.blbw[[1]] <- samples
load("results/cross-val/idm-results-250000-iterations-BLBW-species-2-holdout-2021-08-04.R")
samples.blbw[[2]] <- samples
# Not enough data for BLPW. 
# BTBW 
samples.btbw <- list()
load("results/cross-val/idm-results-250000-iterations-BTBW-species-1-holdout-2021-08-04.R")
samples.btbw[[1]] <- samples
load("results/cross-val/idm-results-250000-iterations-BTBW-species-2-holdout-2021-08-04.R")
samples.btbw[[2]] <- samples
# BTNW 
samples.btnw <- list()
load("results/cross-val/idm-results-250000-iterations-BTNW-species-1-holdout-2021-08-04.R")
samples.btnw[[1]] <- samples
load("results/cross-val/idm-results-250000-iterations-BTNW-species-2-holdout-2021-08-04.R")
samples.btnw[[2]] <- samples
# CAWA 
samples.cawa <- list()
load("results/cross-val/idm-results-250000-iterations-CAWA-species-1-holdout-2021-08-04.R")
samples.cawa[[1]] <- samples
load("results/cross-val/idm-results-250000-iterations-CAWA-species-2-holdout-2021-08-04.R")
samples.cawa[[2]] <- samples
# MAWA 
samples.mawa <- list()
load("results/cross-val/idm-results-250000-iterations-MAWA-species-1-holdout-2021-08-04.R")
samples.mawa[[1]] <- samples
load("results/cross-val/idm-results-250000-iterations-MAWA-species-2-holdout-2021-08-04.R")
samples.mawa[[2]] <- samples
# NAWA 
samples.nawa <- list()
load("results/cross-val/idm-results-250000-iterations-NAWA-species-1-holdout-2021-08-04.R")
samples.nawa[[1]] <- samples
load("results/cross-val/idm-results-250000-iterations-NAWA-species-2-holdout-2021-08-04.R")
samples.nawa[[2]] <- samples
# OVEN 
samples.oven <- list()
load("results/cross-val/idm-results-250000-iterations-OVEN-species-1-holdout-2021-08-04.R")
samples.oven[[1]] <- samples
load("results/cross-val/idm-results-250000-iterations-OVEN-species-2-holdout-2021-08-04.R")
samples.oven[[2]] <- samples
# REVI 
samples.revi <- list()
load("results/cross-val/idm-results-250000-iterations-REVI-species-1-holdout-2021-08-04.R")
samples.revi[[1]] <- samples
load("results/cross-val/idm-results-250000-iterations-REVI-species-2-holdout-2021-08-04.R")
samples.revi[[2]] <- samples

# Define constants --------------------------------------------------------
n.years <- 9
I <- 12
sp <- c('AMRE', 'BAWW', 'BHVI', 'BLBW', 'BLPW', 'BTBW',
        'BTNW', 'CAWA', 'MAWA', 'NAWA', 'OVEN', 'REVI')
sp.names <- c('American Redstart', 'Black-and-white Warbler', 'Blue-headed Vireo',
	      'Blackburnian Warbler', 'Blackpoll Warbler', 'Black-throated Blue Warbler',
	      'Black-throated Green Warbler', 'Canada Warbler', 'Magnolia Warbler',
	      'Nashville Warbler', 'Ovenbird', 'Red-eyed Vireo')
sp.indx <- 1:12
sp.neon.indx <- c(1:4, 6, 7, 8, 9, 11, 12)
I.neon <- length(sp.neon.indx)
n.years.indx <- 1:9
n.years.neon <- 4
n.years.neon.indx <- 6:9
elpd.hbef.vals <- rep(NA, I)
elpd.neon.vals <- rep(NA, I)
elpd.bbs.vals <- rep(NA, I)
# Load in raw data --------------------------------------------------------
load("data/final-bird-data.R")
load("data/cross-val-indices.R")
# Rename things for consistency
hbef.elev <- hb.elev
hbef.for <- hb.for

# The remainder of the code computes the lpd according to the description
# provided in Doser et al. (2021). A single log-pointwide predictive density
# measure is calculated for each species, for each data set location included
# in that model. For these single species models, we only compute
# the lpd for the full model for subsequent comparisons to the full model 
# for the ICOM. 

elev.fit <- c(hbef.elev[-hbef.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(hbef.elev[hbef.x.indices[[curr.set]]], 
	       bbs.elev[bbs.x.indices[[curr.set]]], 
	       neon.elev[neon.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(hbef.for[-hbef.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(hbef.for[hbef.x.indices[[curr.set]]], 
	       bbs.for[bbs.x.indices[[curr.set]]], 
	       neon.for[neon.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.years <- 9
n.iter <- 500
my.iter <- (nrow(samples.amre[[1]]) - n.iter + 1):nrow(samples.amre[[1]])
J.hbef <- length(hbef.x.indices[[curr.set]])
J.neon <- length(neon.x.indices[[curr.set]])
J.bbs <- length(bbs.x.indices[[curr.set]])
elev.hbef.pred <- elev.pred[1:J.hbef]
elev.neon.pred <- elev.pred[1:J.neon + J.hbef]
elev.bbs.pred <- elev.pred[1:J.bbs + J.neon + J.hbef]
for.hbef.pred <- for.pred[1:J.hbef]
for.neon.pred <- for.pred[1:J.neon + J.hbef]
for.bbs.pred <- for.pred[1:J.bbs + J.neon + J.hbef]

# AMRE --------------------------------------------------------------------
samples.fit <- samples.amre[[curr.set]]
samples.pred <- samples.amre[[pred.set]]

hbef.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.hbef.pred, 
		    for.pred = for.hbef.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    n.years = n.years, type = 'hbef')
elpd.hbef.vals[1] <- sum(log(hbef.elpd.vals.1))

neon.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		         samples.pred = samples.pred, 
		         elev.pred = elev.neon.pred, 
		         for.pred = for.neon.pred, 
		         my.iter = my.iter, 
		         n.iter = n.iter, 
		         n.years = n.years, type = 'neon')
elpd.neon.vals[1] <- sum(log(neon.elpd.vals.1))

bbs.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.bbs.pred, 
		    for.pred = for.bbs.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    n.years = n.years, type = 'bbs')

elpd.bbs.vals[1] <- sum(log(bbs.elpd.vals.1))
# BAWW --------------------------------------------------------------------
samples.fit <- samples.baww[[curr.set]]
samples.pred <- samples.baww[[pred.set]]

hbef.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.hbef.pred, 
		    for.pred = for.hbef.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    n.years = n.years, type = 'hbef')
elpd.hbef.vals[2] <- sum(log(hbef.elpd.vals.1))

neon.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		         samples.pred = samples.pred, 
		         elev.pred = elev.neon.pred, 
		         for.pred = for.neon.pred, 
		         my.iter = my.iter, 
		         n.iter = n.iter, 
		         n.years = n.years, type = 'neon')
elpd.neon.vals[2] <- sum(log(neon.elpd.vals.1))

bbs.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.bbs.pred, 
		    for.pred = for.bbs.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    n.years = n.years, type = 'bbs')

elpd.bbs.vals[2] <- sum(log(bbs.elpd.vals.1))
# BHVI --------------------------------------------------------------------
samples.fit <- samples.bhvi[[curr.set]]
samples.pred <- samples.bhvi[[pred.set]]

hbef.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.hbef.pred, 
		    for.pred = for.hbef.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    n.years = n.years, type = 'hbef')
elpd.hbef.vals[3] <- sum(log(hbef.elpd.vals.1))

neon.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		         samples.pred = samples.pred, 
		         elev.pred = elev.neon.pred, 
		         for.pred = for.neon.pred, 
		         my.iter = my.iter, 
		         n.iter = n.iter, 
		         n.years = n.years, type = 'neon')
elpd.neon.vals[3] <- sum(log(neon.elpd.vals.1))

bbs.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.bbs.pred, 
		    for.pred = for.bbs.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    n.years = n.years, type = 'bbs')

elpd.bbs.vals[3] <- sum(log(bbs.elpd.vals.1))
# BLBW --------------------------------------------------------------------
samples.fit <- samples.blbw[[curr.set]]
samples.pred <- samples.blbw[[pred.set]]

hbef.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.hbef.pred, 
		    for.pred = for.hbef.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    n.years = n.years, type = 'hbef')
elpd.hbef.vals[4] <- sum(log(hbef.elpd.vals.1))

neon.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		         samples.pred = samples.pred, 
		         elev.pred = elev.neon.pred, 
		         for.pred = for.neon.pred, 
		         my.iter = my.iter, 
		         n.iter = n.iter, 
		         n.years = n.years, type = 'neon')
elpd.neon.vals[4] <- sum(log(neon.elpd.vals.1))

bbs.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.bbs.pred, 
		    for.pred = for.bbs.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    n.years = n.years, type = 'bbs')

elpd.bbs.vals[4] <- sum(log(bbs.elpd.vals.1))
# BTBW --------------------------------------------------------------------
samples.fit <- samples.btbw[[curr.set]]
samples.pred <- samples.btbw[[pred.set]]

hbef.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.hbef.pred, 
		    for.pred = for.hbef.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    n.years = n.years, type = 'hbef')
elpd.hbef.vals[6] <- sum(log(hbef.elpd.vals.1))

neon.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		         samples.pred = samples.pred, 
		         elev.pred = elev.neon.pred, 
		         for.pred = for.neon.pred, 
		         my.iter = my.iter, 
		         n.iter = n.iter, 
		         n.years = n.years, type = 'neon')
elpd.neon.vals[6] <- sum(log(neon.elpd.vals.1))

bbs.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.bbs.pred, 
		    for.pred = for.bbs.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    n.years = n.years, type = 'bbs')

elpd.bbs.vals[6] <- sum(log(bbs.elpd.vals.1))
# BTNW --------------------------------------------------------------------
samples.fit <- samples.btnw[[curr.set]]
samples.pred <- samples.btnw[[pred.set]]

hbef.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.hbef.pred, 
		    for.pred = for.hbef.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    n.years = n.years, type = 'hbef')
elpd.hbef.vals[7] <- sum(log(hbef.elpd.vals.1))

neon.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		         samples.pred = samples.pred, 
		         elev.pred = elev.neon.pred, 
		         for.pred = for.neon.pred, 
		         my.iter = my.iter, 
		         n.iter = n.iter, 
		         n.years = n.years, type = 'neon')
elpd.neon.vals[7] <- sum(log(neon.elpd.vals.1))

bbs.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.bbs.pred, 
		    for.pred = for.bbs.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    n.years = n.years, type = 'bbs')

elpd.bbs.vals[7] <- sum(log(bbs.elpd.vals.1))
# CAWA --------------------------------------------------------------------
samples.fit <- samples.cawa[[curr.set]]
samples.pred <- samples.cawa[[pred.set]]

hbef.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.hbef.pred, 
		    for.pred = for.hbef.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    n.years = n.years, type = 'hbef')
elpd.hbef.vals[8] <- sum(log(hbef.elpd.vals.1))

neon.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		         samples.pred = samples.pred, 
		         elev.pred = elev.neon.pred, 
		         for.pred = for.neon.pred, 
		         my.iter = my.iter, 
		         n.iter = n.iter, 
		         n.years = n.years, type = 'neon')
elpd.neon.vals[8] <- sum(log(neon.elpd.vals.1))

bbs.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.bbs.pred, 
		    for.pred = for.bbs.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    n.years = n.years, type = 'bbs')

elpd.bbs.vals[8] <- sum(log(bbs.elpd.vals.1))
# MAWA --------------------------------------------------------------------
samples.fit <- samples.mawa[[curr.set]]
samples.pred <- samples.mawa[[pred.set]]

hbef.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.hbef.pred, 
		    for.pred = for.hbef.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    n.years = n.years, type = 'hbef')
elpd.hbef.vals[9] <- sum(log(hbef.elpd.vals.1))

neon.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		         samples.pred = samples.pred, 
		         elev.pred = elev.neon.pred, 
		         for.pred = for.neon.pred, 
		         my.iter = my.iter, 
		         n.iter = n.iter, 
		         n.years = n.years, type = 'neon')
elpd.neon.vals[9] <- sum(log(neon.elpd.vals.1))

bbs.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.bbs.pred, 
		    for.pred = for.bbs.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    n.years = n.years, type = 'bbs')

elpd.bbs.vals[9] <- sum(log(bbs.elpd.vals.1))
# NAWA --------------------------------------------------------------------
samples.fit <- samples.nawa[[curr.set]]
samples.pred <- samples.nawa[[pred.set]]

hbef.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.hbef.pred, 
		    for.pred = for.hbef.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    n.years = n.years, type = 'hbef')
elpd.hbef.vals[10] <- sum(log(hbef.elpd.vals.1))

bbs.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.bbs.pred, 
		    for.pred = for.bbs.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    n.years = n.years, type = 'bbs')

elpd.bbs.vals[10] <- sum(log(bbs.elpd.vals.1))
# OVEN --------------------------------------------------------------------
samples.fit <- samples.oven[[curr.set]]
samples.pred <- samples.oven[[pred.set]]

hbef.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.hbef.pred, 
		    for.pred = for.hbef.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    n.years = n.years, type = 'hbef')
elpd.hbef.vals[11] <- sum(log(hbef.elpd.vals.1))

neon.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		         samples.pred = samples.pred, 
		         elev.pred = elev.neon.pred, 
		         for.pred = for.neon.pred, 
		         my.iter = my.iter, 
		         n.iter = n.iter, 
		         n.years = n.years, type = 'neon')
elpd.neon.vals[11] <- sum(log(neon.elpd.vals.1))

bbs.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.bbs.pred, 
		    for.pred = for.bbs.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    n.years = n.years, type = 'bbs')

elpd.bbs.vals[11] <- sum(log(bbs.elpd.vals.1))
# REVI --------------------------------------------------------------------
samples.fit <- samples.revi[[curr.set]]
samples.pred <- samples.revi[[pred.set]]

hbef.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.hbef.pred, 
		    for.pred = for.hbef.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    n.years = n.years, type = 'hbef')
elpd.hbef.vals[12] <- sum(log(hbef.elpd.vals.1))

neon.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		         samples.pred = samples.pred, 
		         elev.pred = elev.neon.pred, 
		         for.pred = for.neon.pred, 
		         my.iter = my.iter, 
		         n.iter = n.iter, 
		         n.years = n.years, type = 'neon')
elpd.neon.vals[12] <- sum(log(neon.elpd.vals.1))

bbs.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.bbs.pred, 
		    for.pred = for.bbs.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    n.years = n.years, type = 'bbs')

elpd.bbs.vals[12] <- sum(log(bbs.elpd.vals.1))

save(elpd.hbef.vals, elpd.bbs.vals, elpd.neon.vals, 
     file = paste('results/cross-val/cross-val-idm-results-', curr.set, '.R', sep = ''))
