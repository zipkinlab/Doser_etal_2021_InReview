# cross-val.R: code to obtain predictive results from cross-validation
#              for the White Mountain National Forest case study. 
# Author: Jeffrey W. Doser (doserjef@msu.edu)
# Citation: 

rm(list = ls())
library(coda)
library(tidyverse)

# Functions ---------------------------------------------------------------
logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}
# These functions compute log pointwise-posterior predictive densities that
# are used to assess predictive performance. 
source("code/cross-val-icom/lpd-sampler.R")
source("code/cross-val-icom/lpd-fit-neon-sampler.R")
source("code/cross-val-icom/lpd-pred-neon-sampler.R")

# Current hold-out set ----------------------------------------------------
# Takes value 1 or 2 depending on which set is the current hold out set
# Provided as user input when running script from terminal. 
#curr.set <- as.numeric(commandArgs(trailingOnly = TRUE))
curr.set <- 2
pred.set <- ifelse(curr.set == 1, 2, 1)
# Read in all data sets ---------------------------------------------------
# HBEF --------------------------------
samples.hbef <- list()
load("results/cross-val/icom-HBEF-results-350000-iterations-1-holdout-2021-08-20.R")
samples.hbef[[1]] <- samples
load("results/cross-val/icom-HBEF-results-350000-iterations-2-holdout-2021-08-20.R")
samples.hbef[[2]] <- samples
# NEON --------------------------------
samples.neon <- list()
load("results/cross-val/icom-NEON-results-350000-iterations-1-holdout-2021-08-19.R")
samples.neon[[1]] <- samples
load("results/cross-val/icom-NEON-results-350000-iterations-2-holdout-2021-08-19.R")
samples.neon[[2]] <- samples
# BBS ---------------------------------
samples.bbs <- list()
load("results/cross-val/icom-bbs-results-350000-iterations-1-holdout-2021-08-19.R")
samples.bbs[[1]] <- samples
load("results/cross-val/icom-bbs-results-350000-iterations-2-holdout-2021-08-19.R")
samples.bbs[[2]] <- samples
# HBEF + NEON -------------------------
samples.hbef.neon <- list()
load("results/cross-val/icom-HBEF-NEON-results-350000-iterations-1-holdout-2021-08-20.R")
samples.hbef.neon[[1]] <- samples
load("results/cross-val/icom-HBEF-NEON-results-350000-iterations-2-holdout-2021-08-20.R")
samples.hbef.neon[[2]] <- samples
# HBEF + BBS --------------------------
samples.hbef.bbs <- list()
load("results/cross-val/icom-HBEF-BBS-results-350000-iterations-1-holdout-2021-08-20.R")
samples.hbef.bbs[[1]] <- samples
load("results/cross-val/icom-HBEF-BBS-results-350000-iterations-2-holdout-2021-08-20.R")
samples.hbef.bbs[[2]] <- samples
# NEON + BBS --------------------------
samples.neon.bbs <- list()
load("results/cross-val/icom-NEON-BBS-results-350000-iterations-1-holdout-2021-08-19.R")
samples.neon.bbs[[1]] <- samples
load("results/cross-val/icom-NEON-BBS-results-350000-iterations-2-holdout-2021-08-19.R")
samples.neon.bbs[[2]] <- samples
# HBEF + NEON + BBS -------------------
samples.hbef.neon.bbs <- list()
load("results/cross-val/icom-HBEF-NEON-BBS-results-350000-iterations-1-holdout-2021-08-20.R")
samples.hbef.neon.bbs[[1]] <- samples
load("results/cross-val/icom-HBEF-NEON-BBS-results-350000-iterations-2-holdout-2021-08-20.R")
samples.hbef.neon.bbs[[2]] <- samples
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
# Initiate arrays to store lpd values. 
elpd.hbef.vals <- array(NA, dim = c(I, 4, 7))
elpd.neon.vals <- array(NA, dim = c(I, 4, 7))
elpd.bbs.vals <- array(NA, dim = c(I, 4, 7))

# Load in raw data --------------------------------------------------------
load("data/final-bird-data.R")
load("data/cross-val-indices.R")
# Rename things for consistency
hbef.elev <- hb.elev
hbef.for <- hb.for

# The remainder of the code computes the lpd according to the description 
# provided in Doser et al. (2021). A single log-pointwide predictive density
# measure is calculated for each model, for each species, for each data set location included
# in that model, and for each model that could have potentially generated
# those data. This approach accounts for model uncertainty and enables 
# assessment of predictive performance for individual species, the entire 
# community, and at different parts of the entire region of interest. 

# HBEF as predicting model -------------------------------------------------
# HBEF ---------------------------------
samples.fit <- samples.hbef[[curr.set]]
samples.pred <- samples.hbef[[pred.set]]
elev.fit <- c(hbef.elev[-hbef.x.indices[[curr.set]]])
mean.elev.fit <- mean(elev.fit)
sd.elev.fit <- sd(elev.fit)
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(hbef.elev[hbef.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean.elev.fit) / sd.elev.fit 
for.fit <- c(hbef.for[-hbef.x.indices[[curr.set]]])
mean.for.fit <- mean(for.fit)
sd.for.fit <- sd(for.fit)
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(hbef.for[hbef.x.indices[[curr.set]]])
for.pred <- (for.pred - mean.for.fit) / sd.for.fit
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.hbef <- length(hbef.x.indices[[curr.set]])
elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.pred, 
		    for.pred = for.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'hbef')
elpd.hbef.vals[, 1, 1] <- apply(log(elpd.vals.1), 1, sum)

# NEON ---------------------------------
samples.fit <- samples.hbef[[curr.set]]
samples.pred <- samples.neon[[pred.set]]
elev.fit <- c(hbef.elev[-hbef.x.indices[[curr.set]]])
mean.elev.fit <- mean(elev.fit)
sd.elev.fit <- sd(elev.fit)
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(neon.elev[neon.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean.elev.fit) / sd.elev.fit 
for.fit <- c(hbef.for[-hbef.x.indices[[curr.set]]])
mean.for.fit <- mean(for.fit)
sd.for.fit <- sd(for.fit)
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(neon.for[neon.x.indices[[curr.set]]])
for.pred <- (for.pred - mean.for.fit) / sd.for.fit
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.neon <- length(neon.x.indices[[curr.set]])
I.neon <- 10
n.years.neon <- 4
elpd.vals.1 <- elpd.neon.pred(samples.fit = samples.fit, 
		              samples.pred = samples.pred, 
		              elev.pred = elev.pred, 
		              for.pred = for.pred, 
		              my.iter = my.iter, 
		              n.iter = n.iter, 
		              I = I.neon, I.full = I, n.years = n.years)
elpd.neon.vals[sp.neon.indx, 1, 1] <- apply((log(elpd.vals.1)), 1, sum)

# BBS ---------------------------------
samples.fit <- samples.hbef[[curr.set]]
samples.pred <- samples.bbs[[pred.set]]
elev.fit <- c(hbef.elev[-hbef.x.indices[[curr.set]]])
mean.elev.fit <- mean(elev.fit)
sd.elev.fit <- sd(elev.fit)
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(bbs.elev[bbs.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean.elev.fit) / sd.elev.fit 
for.fit <- c(hbef.for[-hbef.x.indices[[curr.set]]])
mean.for.fit <- mean(for.fit)
sd.for.fit <- sd(for.fit)
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(bbs.for[bbs.x.indices[[curr.set]]])
for.pred <- (for.pred - mean.for.fit) / sd.for.fit
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.bbs <- length(bbs.x.indices[[curr.set]])
elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.pred, 
		    for.pred = for.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'bbs')
elpd.bbs.vals[, 1, 1] <- apply(log(elpd.vals.1), 1, sum)

# HBEF + NEON -------------------------
samples.fit <- samples.hbef[[curr.set]]
samples.pred <- samples.hbef.neon[[pred.set]]
elev.fit <- c(hbef.elev[-hbef.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(hbef.elev[hbef.x.indices[[curr.set]]],
	       neon.elev[neon.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(hbef.for[-hbef.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(hbef.for[hbef.x.indices[[curr.set]]],
	       neon.for[neon.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.years <- 9
I <- 12
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.hbef <- length(hbef.x.indices[[curr.set]])
J.neon <- length(neon.x.indices[[curr.set]])
elev.hbef.pred <- elev.pred[1:J.hbef]
elev.neon.pred <- elev.pred[1:J.neon + J.hbef]
for.hbef.pred <- for.pred[1:J.hbef]
for.neon.pred <- for.pred[1:J.neon + J.hbef]

hbef.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		        samples.pred = samples.pred, 
		        elev.pred = elev.hbef.pred, 
		        for.pred = for.hbef.pred, 
		        my.iter = my.iter, 
		        n.iter = n.iter, 
		        I = I, n.years = n.years, type = 'hbef')
elpd.hbef.vals[, 2, 1] <- apply(log(hbef.elpd.vals.1), 1, sum)


neon.elpd.vals.1 <- elpd.neon.pred(samples.fit = samples.fit, 
		                   samples.pred = samples.pred, 
		                   elev.pred = elev.neon.pred, 
		                   for.pred = for.neon.pred, 
		                   my.iter = my.iter, 
		                   n.iter = n.iter, 
		                   I = I, I.full = I, n.years = n.years)

elpd.neon.vals[sp.neon.indx, 2, 1] <- apply(log(neon.elpd.vals.1), 1, sum)
# HBEF + BBS --------------------------
samples.fit <- samples.hbef[[curr.set]]
samples.pred <- samples.hbef.bbs[[pred.set]]
elev.fit <- c(hbef.elev[-hbef.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(hbef.elev[hbef.x.indices[[curr.set]]],
	       bbs.elev[bbs.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(hbef.for[-hbef.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(hbef.for[hbef.x.indices[[curr.set]]],
	       bbs.for[bbs.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.years <- 9
I <- 12
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.hbef <- length(hbef.x.indices[[curr.set]])
J.bbs <- length(bbs.x.indices[[curr.set]])
elev.hbef.pred <- elev.pred[1:J.hbef]
elev.bbs.pred <- elev.pred[1:J.bbs + J.hbef]
for.hbef.pred <- for.pred[1:J.hbef]
for.bbs.pred <- for.pred[1:J.bbs + J.hbef]

hbef.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.hbef.pred, 
		    for.pred = for.hbef.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'hbef')
elpd.hbef.vals[, 3, 1] <- apply(log(hbef.elpd.vals.1), 1, sum)


bbs.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.bbs.pred, 
		    for.pred = for.bbs.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'bbs')
elpd.bbs.vals[, 2, 1] <- apply(log(bbs.elpd.vals.1), 1, sum)
# NEON + BBS --------------------------
samples.fit <- samples.hbef[[curr.set]]
samples.pred <- samples.neon.bbs[[pred.set]]
elev.fit <- c(hbef.elev[-hbef.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(neon.elev[neon.x.indices[[curr.set]]],
	       bbs.elev[bbs.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(hbef.for[-hbef.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(neon.for[neon.x.indices[[curr.set]]],
	       bbs.for[bbs.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.years <- 9
I <- 12
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.neon <- length(neon.x.indices[[curr.set]])
J.bbs <- length(bbs.x.indices[[curr.set]])
elev.neon.pred <- elev.pred[1:J.neon]
elev.bbs.pred <- elev.pred[1:J.bbs + J.neon]
for.neon.pred <- for.pred[1:J.neon]
for.bbs.pred <- for.pred[1:J.bbs + J.neon]

neon.elpd.vals.1 <- elpd.neon.pred(samples.fit = samples.fit, 
		                   samples.pred = samples.pred, 
		                   elev.pred = elev.neon.pred, 
		                   for.pred = for.neon.pred, 
		                   my.iter = my.iter, 
		                   n.iter = n.iter, 
		                   I = I, I.full = I, n.years = n.years)
elpd.neon.vals[sp.neon.indx, 3, 1] <- apply(log(neon.elpd.vals.1), 1, sum)

bbs.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.bbs.pred, 
		    for.pred = for.bbs.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'bbs')

elpd.bbs.vals[, 3, 1] <- apply(log(bbs.elpd.vals.1), 1, sum)

# HBEF + NEON + BBS -------------------
samples.fit <- samples.hbef[[curr.set]]
samples.pred <- samples.hbef.neon.bbs[[pred.set]]
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
I <- 12
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.hbef <- length(hbef.x.indices[[curr.set]])
J.neon <- length(neon.x.indices[[curr.set]])
J.bbs <- length(bbs.x.indices[[curr.set]])
elev.hbef.pred <- elev.pred[1:J.hbef]
elev.neon.pred <- elev.pred[1:J.neon + J.hbef]
elev.bbs.pred <- elev.pred[1:J.bbs + J.neon + J.hbef]
for.hbef.pred <- for.pred[1:J.hbef]
for.neon.pred <- for.pred[1:J.neon + J.hbef]
for.bbs.pred <- for.pred[1:J.bbs + J.neon + J.hbef]

hbef.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.hbef.pred, 
		    for.pred = for.hbef.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'hbef')
elpd.hbef.vals[, 4, 1] <- apply(log(hbef.elpd.vals.1), 1, sum)

neon.elpd.vals.1 <- elpd.neon.pred(samples.fit = samples.fit, 
		                   samples.pred = samples.pred, 
		                   elev.pred = elev.neon.pred, 
		                   for.pred = for.neon.pred, 
		                   my.iter = my.iter, 
		                   n.iter = n.iter, 
		                   I = I, I.full = I, n.years = n.years)
elpd.neon.vals[sp.neon.indx, 4, 1] <- apply(log(neon.elpd.vals.1), 1, sum)

bbs.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.bbs.pred, 
		    for.pred = for.bbs.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'bbs')

elpd.bbs.vals[, 4, 1] <- apply(log(bbs.elpd.vals.1), 1, sum)

# Only assess the NEON only model to predict the NEON data set itself, 
# since it is from a smaller time frame. 
# NEON as predicting model -------------------------------------------------
# NEON ---------------------------------
samples.fit <- samples.neon[[curr.set]]
samples.pred <- samples.neon[[pred.set]]
elev.fit <- c(neon.elev[-neon.x.indices[[curr.set]]])
mean.elev.fit <- mean(elev.fit)
sd.elev.fit <- sd(elev.fit)
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(neon.elev[neon.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean.elev.fit) / sd.elev.fit 
for.fit <- c(neon.for[-neon.x.indices[[curr.set]]])
mean.for.fit <- mean(for.fit)
sd.for.fit <- sd(for.fit)
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(neon.for[neon.x.indices[[curr.set]]])
for.pred <- (for.pred - mean.for.fit) / sd.for.fit
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.neon <- length(neon.x.indices[[curr.set]])
elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.pred, 
		    for.pred = for.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I.neon, n.years = n.years.neon, type = 'neon')
elpd.neon.vals[sp.neon.indx, 1, 2] <- apply(log(elpd.vals.1), 1, sum)

# HBEF + NEON -------------------------
samples.fit <- samples.neon[[curr.set]]
samples.pred <- samples.hbef.neon[[pred.set]]
elev.fit <- c(neon.elev[-neon.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(hbef.elev[hbef.x.indices[[curr.set]]],
	       neon.elev[neon.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(neon.for[-neon.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(hbef.for[hbef.x.indices[[curr.set]]],
	       neon.for[neon.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.years <- 9
I <- 12
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.neon <- length(neon.x.indices[[curr.set]])
elev.neon.pred <- elev.pred[1:J.neon + J.hbef]
for.neon.pred <- for.pred[1:J.neon + J.hbef]

neon.elpd.vals.1 <- elpd.neon.fit.pred(samples.fit = samples.fit,
		                       samples.pred = samples.pred,
		                       elev.pred = elev.neon.pred,
		                       for.pred = for.neon.pred,
		                       my.iter = my.iter,
		                       n.iter = n.iter,
		                       I = I, n.years = n.years)

elpd.neon.vals[sp.neon.indx, 2, 2] <- apply(log(neon.elpd.vals.1), 1, sum)

# NEON + BBS --------------------------
samples.fit <- samples.neon[[curr.set]]
samples.pred <- samples.neon.bbs[[pred.set]]
elev.fit <- c(neon.elev[-neon.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(neon.elev[neon.x.indices[[curr.set]]],
	       bbs.elev[bbs.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(neon.for[-neon.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(neon.for[neon.x.indices[[curr.set]]],
	       bbs.for[bbs.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.years <- 9
I <- 12
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.neon <- length(neon.x.indices[[curr.set]])
elev.neon.pred <- elev.pred[1:J.neon]
for.neon.pred <- for.pred[1:J.neon]

neon.elpd.vals.1 <- elpd.neon.fit.pred(samples.fit = samples.fit, 
		                       samples.pred = samples.pred, 
		                       elev.pred = elev.neon.pred, 
		                       for.pred = for.neon.pred, 
		                       my.iter = my.iter, 
		                       n.iter = n.iter, 
		                       I = I, n.years = n.years)
elpd.neon.vals[sp.neon.indx, 3, 2] <- apply(log(neon.elpd.vals.1), 1, sum)

# HBEF + NEON + BBS -------------------
samples.fit <- samples.neon[[curr.set]]
samples.pred <- samples.hbef.neon.bbs[[curr.set]]
elev.fit <- c(neon.elev[-neon.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(hbef.elev[hbef.x.indices[[curr.set]]], 
	       bbs.elev[bbs.x.indices[[curr.set]]], 
	       neon.elev[neon.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(neon.for[-neon.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(hbef.for[hbef.x.indices[[curr.set]]], 
	       bbs.for[bbs.x.indices[[curr.set]]], 
	       neon.for[neon.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.years <- 9
I <- 12
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.neon <- length(neon.x.indices[[curr.set]])
elev.neon.pred <- elev.pred[1:J.neon + J.hbef]
for.neon.pred <- for.pred[1:J.neon + J.hbef]

neon.elpd.vals.1 <- elpd.neon.fit.pred(samples.fit = samples.fit, 
		                       samples.pred = samples.pred, 
		                       elev.pred = elev.neon.pred, 
		                       for.pred = for.neon.pred, 
		                       my.iter = my.iter, 
		                       n.iter = n.iter, 
		                       I = I, n.years = n.years)
elpd.neon.vals[sp.neon.indx, 4, 2] <- apply(log(neon.elpd.vals.1), 1, sum)



# BBS as predicting model -------------------------------------------------
# HBEF ---------------------------------
samples.fit <- samples.bbs[[curr.set]]
samples.pred <- samples.hbef[[curr.set]]
elev.fit <- c(bbs.elev[-bbs.x.indices[[curr.set]]])
mean.elev.fit <- mean(elev.fit)
sd.elev.fit <- sd(elev.fit)
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(hbef.elev[hbef.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean.elev.fit) / sd.elev.fit
for.fit <- c(bbs.for[-bbs.x.indices[[curr.set]]])
mean.for.fit <- mean(for.fit)
sd.for.fit <- sd(for.fit)
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(hbef.for[hbef.x.indices[[curr.set]]])
for.pred <- (for.pred - mean.for.fit) / sd.for.fit
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.hbef <- length(hbef.x.indices[[curr.set]])
elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.pred, 
		    for.pred = for.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'hbef')
elpd.hbef.vals[, 1, 3] <- apply(log(elpd.vals.1), 1, sum)

# NEON ---------------------------------
samples.fit <- samples.bbs[[curr.set]]
samples.pred <- samples.neon[[pred.set]]
elev.fit <- c(bbs.elev[-bbs.x.indices[[curr.set]]])
mean.elev.fit <- mean(elev.fit)
sd.elev.fit <- sd(elev.fit)
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(neon.elev[neon.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean.elev.fit) / sd.elev.fit
for.fit <- c(bbs.for[-bbs.x.indices[[curr.set]]])
mean.for.fit <- mean(for.fit)
sd.for.fit <- sd(for.fit)
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(neon.for[neon.x.indices[[curr.set]]])
for.pred <- (for.pred - mean.for.fit) / sd.for.fit
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.neon <- length(neon.x.indices[[curr.set]])
I.neon <- 10
n.years.neon <- 4
elpd.vals.1 <- elpd.neon.pred(samples.fit = samples.fit, 
		              samples.pred = samples.pred, 
		              elev.pred = elev.pred, 
		              for.pred = for.pred, 
		              my.iter = my.iter, 
		              n.iter = n.iter, 
		              I = I.neon, I.full = I, n.years = n.years)
elpd.neon.vals[sp.neon.indx, 1, 3] <- apply((log(elpd.vals.1)), 1, sum)

# BBS ---------------------------------
samples.fit <- samples.bbs[[curr.set]]
samples.pred <- samples.bbs[[pred.set]]
elev.fit <- c(bbs.elev[-bbs.x.indices[[curr.set]]])
mean.elev.fit <- mean(elev.fit)
sd.elev.fit <- sd(elev.fit)
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(bbs.elev[bbs.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean.elev.fit) / sd.elev.fit
for.fit <- c(bbs.for[-bbs.x.indices[[curr.set]]])
mean.for.fit <- mean(for.fit)
sd.for.fit <- sd(for.fit)
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(bbs.for[bbs.x.indices[[curr.set]]])
for.pred <- (for.pred - mean.for.fit) / sd.for.fit
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.bbs <- length(bbs.x.indices[[curr.set]])
elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.pred, 
		    for.pred = for.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'bbs')
elpd.bbs.vals[, 1, 3] <- apply(log(elpd.vals.1), 1, sum)

# HBEF + NEON -------------------------
samples.fit <- samples.bbs[[curr.set]]
samples.pred <- samples.hbef.neon[[pred.set]]
elev.fit <- c(bbs.elev[-bbs.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(hbef.elev[hbef.x.indices[[curr.set]]],
	       neon.elev[neon.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(bbs.for[-bbs.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(hbef.for[hbef.x.indices[[curr.set]]],
	       neon.for[neon.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.years <- 9
I <- 12
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.hbef <- length(hbef.x.indices[[curr.set]])
J.neon <- length(neon.x.indices[[curr.set]])
elev.hbef.pred <- elev.pred[1:J.hbef]
elev.neon.pred <- elev.pred[1:J.neon + J.hbef]
for.hbef.pred <- for.pred[1:J.hbef]
for.neon.pred <- for.pred[1:J.neon + J.hbef]

hbef.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		        samples.pred = samples.pred, 
		        elev.pred = elev.hbef.pred, 
		        for.pred = for.hbef.pred, 
		        my.iter = my.iter, 
		        n.iter = n.iter, 
		        I = I, n.years = n.years, type = 'hbef')
elpd.hbef.vals[, 2, 3] <- apply(log(hbef.elpd.vals.1), 1, sum)


neon.elpd.vals.1 <- elpd.neon.pred(samples.fit = samples.fit, 
		                   samples.pred = samples.pred, 
		                   elev.pred = elev.neon.pred, 
		                   for.pred = for.neon.pred, 
		                   my.iter = my.iter, 
		                   n.iter = n.iter, 
		                   I = I, I.full = I, n.years = n.years)

elpd.neon.vals[sp.neon.indx, 2, 3] <- apply(log(neon.elpd.vals.1), 1, sum)
# HBEF + BBS --------------------------
samples.fit <- samples.bbs[[curr.set]]
samples.pred <- samples.hbef.bbs[[pred.set]]
elev.fit <- c(bbs.elev[-bbs.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(hbef.elev[hbef.x.indices[[curr.set]]],
	       bbs.elev[bbs.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(bbs.for[-bbs.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(hbef.for[hbef.x.indices[[curr.set]]],
	       bbs.for[bbs.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.years <- 9
I <- 12
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)

J.hbef <- length(hbef.x.indices[[curr.set]])
J.bbs <- length(bbs.x.indices[[curr.set]])
elev.hbef.pred <- elev.pred[1:J.hbef]
elev.bbs.pred <- elev.pred[1:J.bbs + J.hbef]
for.hbef.pred <- for.pred[1:J.hbef]
for.bbs.pred <- for.pred[1:J.bbs + J.hbef]

hbef.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.hbef.pred, 
		    for.pred = for.hbef.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'hbef')
elpd.hbef.vals[, 3, 3] <- apply(log(hbef.elpd.vals.1), 1, sum)


bbs.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.bbs.pred, 
		    for.pred = for.bbs.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'bbs')
elpd.bbs.vals[, 2, 3] <- apply(log(bbs.elpd.vals.1), 1, sum)
# NEON + BBS --------------------------
samples.fit <- samples.bbs[[curr.set]]
samples.pred <- samples.neon.bbs[[pred.set]]
elev.fit <- c(bbs.elev[-bbs.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(neon.elev[neon.x.indices[[curr.set]]],
	       bbs.elev[bbs.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(bbs.for[-bbs.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(neon.for[neon.x.indices[[curr.set]]],
	       bbs.for[bbs.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.years <- 9
I <- 12
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.neon <- length(neon.x.indices[[curr.set]])
J.bbs <- length(bbs.x.indices[[curr.set]])
elev.neon.pred <- elev.pred[1:J.neon]
elev.bbs.pred <- elev.pred[1:J.bbs + J.neon]
for.neon.pred <- for.pred[1:J.neon]
for.bbs.pred <- for.pred[1:J.bbs + J.neon]

neon.elpd.vals.1 <- elpd.neon.pred(samples.fit = samples.fit, 
		                   samples.pred = samples.pred, 
		                   elev.pred = elev.neon.pred, 
		                   for.pred = for.neon.pred, 
		                   my.iter = my.iter, 
		                   n.iter = n.iter, 
		                   I = I, I.full = I, n.years = n.years)
elpd.neon.vals[sp.neon.indx, 3, 3] <- apply(log(neon.elpd.vals.1), 1, sum)

bbs.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.bbs.pred, 
		    for.pred = for.bbs.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'bbs')

elpd.bbs.vals[, 3, 3] <- apply(log(bbs.elpd.vals.1), 1, sum)
# HBEF + NEON + BBS -------------------
samples.fit <- samples.bbs[[curr.set]]
samples.pred <- samples.hbef.neon.bbs[[pred.set]]
elev.fit <- c(bbs.elev[-bbs.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(hbef.elev[hbef.x.indices[[curr.set]]],
	       bbs.elev[bbs.x.indices[[curr.set]]],
	       neon.elev[neon.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(bbs.for[-bbs.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(hbef.for[hbef.x.indices[[curr.set]]],
	       bbs.for[bbs.x.indices[[curr.set]]],
	       neon.for[neon.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.years <- 9
I <- 12
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.hbef <- length(hbef.x.indices[[curr.set]])
J.neon <- length(neon.x.indices[[curr.set]])
J.bbs <- length(bbs.x.indices[[curr.set]])
elev.hbef.pred <- elev.pred[1:J.hbef]
elev.neon.pred <- elev.pred[1:J.neon + J.hbef]
elev.bbs.pred <- elev.pred[1:J.bbs + J.neon + J.hbef]
for.hbef.pred <- for.pred[1:J.hbef]
for.neon.pred <- for.pred[1:J.neon + J.hbef]
for.bbs.pred <- for.pred[1:J.bbs + J.neon + J.hbef]

hbef.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.hbef.pred, 
		    for.pred = for.hbef.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'hbef')
elpd.hbef.vals[, 4, 3] <- apply(log(hbef.elpd.vals.1), 1, sum)

neon.elpd.vals.1 <- elpd.neon.pred(samples.fit = samples.fit, 
		                   samples.pred = samples.pred, 
		                   elev.pred = elev.neon.pred, 
		                   for.pred = for.neon.pred, 
		                   my.iter = my.iter, 
		                   n.iter = n.iter, 
		                   I = I, I.full = I, n.years = n.years)
elpd.neon.vals[sp.neon.indx, 4, 3] <- apply(log(neon.elpd.vals.1), 1, sum)

bbs.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.bbs.pred, 
		    for.pred = for.bbs.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'bbs')

elpd.bbs.vals[, 4, 3] <- apply(log(bbs.elpd.vals.1), 1, sum)


# HBEF + NEON as predicting model -----------------------------------
# HBEF ---------------------------------
samples.fit <- samples.hbef.neon[[curr.set]]
samples.pred <- samples.hbef[[pred.set]]
elev.fit <- c(hbef.elev[-hbef.x.indices[[curr.set]]], 
	      neon.elev[-neon.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(hbef.elev[hbef.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(hbef.for[-hbef.x.indices[[curr.set]]], 
	      neon.for[-neon.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(hbef.for[hbef.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.hbef <- length(hbef.x.indices[[curr.set]])
elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.pred, 
		    for.pred = for.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'hbef')
elpd.hbef.vals[, 1, 4] <- apply(log(elpd.vals.1), 1, sum)
# NEON ---------------------------------
samples.fit <- samples.hbef.neon[[curr.set]]
samples.pred <- samples.neon[[pred.set]]
elev.fit <- c(hbef.elev[-hbef.x.indices[[curr.set]]], 
	      neon.elev[-neon.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(neon.elev[neon.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(hbef.for[-hbef.x.indices[[curr.set]]], 
	      neon.for[-neon.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(neon.for[neon.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.neon <- length(neon.x.indices[[curr.set]])
I.neon <- 10
n.years.neon <- 4
elpd.vals.1 <- elpd.neon.pred(samples.fit = samples.fit, 
		              samples.pred = samples.pred, 
		              elev.pred = elev.pred, 
		              for.pred = for.pred, 
		              my.iter = my.iter, 
		              n.iter = n.iter, 
		              I = I.neon, I.full = I, n.years = n.years)
elpd.neon.vals[sp.neon.indx, 1, 4] <- apply((log(elpd.vals.1)), 1, sum)
# BBS ---------------------------------
samples.fit <- samples.hbef.neon[[curr.set]]
samples.pred <- samples.bbs[[pred.set]]
elev.fit <- c(hbef.elev[-hbef.x.indices[[curr.set]]], 
	      neon.elev[-neon.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(bbs.elev[bbs.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(hbef.for[-hbef.x.indices[[curr.set]]], 
	      neon.for[-neon.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(bbs.for[bbs.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.bbs <- length(bbs.x.indices[[curr.set]])
elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.pred, 
		    for.pred = for.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'bbs')
elpd.bbs.vals[, 1, 4] <- apply(log(elpd.vals.1), 1, sum)
# HBEF + NEON -------------------------
samples.fit <- samples.hbef.neon[[curr.set]]
samples.pred <- samples.hbef.neon[[pred.set]]
elev.fit <- c(hbef.elev[-hbef.x.indices[[curr.set]]], 
	      neon.elev[-neon.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(hbef.elev[hbef.x.indices[[curr.set]]], 
	       neon.elev[neon.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(hbef.for[-hbef.x.indices[[curr.set]]], 
	      neon.for[-neon.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(hbef.for[hbef.x.indices[[curr.set]]], 
	       neon.for[neon.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.years <- 9
I <- 12
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.hbef <- length(hbef.x.indices[[curr.set]])
J.neon <- length(neon.x.indices[[curr.set]])
elev.hbef.pred <- elev.pred[1:J.hbef]
elev.neon.pred <- elev.pred[1:J.neon + J.hbef]
for.hbef.pred <- for.pred[1:J.hbef]
for.neon.pred <- for.pred[1:J.neon + J.hbef]

hbef.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		        samples.pred = samples.pred, 
		        elev.pred = elev.hbef.pred, 
		        for.pred = for.hbef.pred, 
		        my.iter = my.iter, 
		        n.iter = n.iter, 
		        I = I, n.years = n.years, type = 'hbef')
elpd.hbef.vals[, 2, 4] <- apply(log(hbef.elpd.vals.1), 1, sum)


neon.elpd.vals.1 <- elpd.neon.pred(samples.fit = samples.fit, 
		                   samples.pred = samples.pred, 
		                   elev.pred = elev.neon.pred, 
		                   for.pred = for.neon.pred, 
		                   my.iter = my.iter, 
		                   n.iter = n.iter, 
		                   I = I, I.full = I, n.years = n.years)

elpd.neon.vals[sp.neon.indx, 2, 4] <- apply(log(neon.elpd.vals.1), 1, sum)
# HBEF + BBS --------------------------
samples.fit <- samples.hbef.neon[[curr.set]]
samples.pred <- samples.hbef.bbs[[pred.set]]
elev.fit <- c(hbef.elev[-hbef.x.indices[[curr.set]]], 
	      neon.elev[-neon.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(hbef.elev[hbef.x.indices[[curr.set]]], 
	       bbs.elev[bbs.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(hbef.for[-hbef.x.indices[[curr.set]]], 
	      neon.for[-neon.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(hbef.for[hbef.x.indices[[curr.set]]], 
	       bbs.for[bbs.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.years <- 9
I <- 12
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.hbef <- length(hbef.x.indices[[curr.set]])
J.bbs <- length(bbs.x.indices[[curr.set]])
elev.hbef.pred <- elev.pred[1:J.hbef]
elev.bbs.pred <- elev.pred[1:J.bbs + J.hbef]
for.hbef.pred <- for.pred[1:J.hbef]
for.bbs.pred <- for.pred[1:J.bbs + J.hbef]

hbef.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.hbef.pred, 
		    for.pred = for.hbef.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'hbef')
elpd.hbef.vals[, 3, 4] <- apply(log(hbef.elpd.vals.1), 1, sum)


bbs.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.bbs.pred, 
		    for.pred = for.bbs.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'bbs')
elpd.bbs.vals[, 2, 4] <- apply(log(bbs.elpd.vals.1), 1, sum)
# NEON + BBS --------------------------
samples.fit <- samples.hbef.neon[[curr.set]]
samples.pred <- samples.neon.bbs[[pred.set]]
elev.fit <- c(hbef.elev[-hbef.x.indices[[curr.set]]], 
	      neon.elev[-neon.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(neon.elev[neon.x.indices[[curr.set]]], 
	       bbs.elev[bbs.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(hbef.for[-hbef.x.indices[[curr.set]]], 
	      neon.for[-neon.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(neon.for[neon.x.indices[[curr.set]]], 
	       bbs.for[bbs.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.years <- 9
I <- 12
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.neon <- length(neon.x.indices[[curr.set]])
J.bbs <- length(bbs.x.indices[[curr.set]])
elev.neon.pred <- elev.pred[1:J.neon]
elev.bbs.pred <- elev.pred[1:J.bbs + J.neon]
for.neon.pred <- for.pred[1:J.neon]
for.bbs.pred <- for.pred[1:J.bbs + J.neon]

neon.elpd.vals.1 <- elpd.neon.pred(samples.fit = samples.fit, 
		                   samples.pred = samples.pred, 
		                   elev.pred = elev.neon.pred, 
		                   for.pred = for.neon.pred, 
		                   my.iter = my.iter, 
		                   n.iter = n.iter, 
		                   I = I, I.full = I, n.years = n.years)
elpd.neon.vals[sp.neon.indx, 3, 4] <- apply(log(neon.elpd.vals.1), 1, sum)

bbs.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.bbs.pred, 
		    for.pred = for.bbs.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'bbs')

elpd.bbs.vals[, 3, 4] <- apply(log(bbs.elpd.vals.1), 1, sum)
# HBEF + NEON + BBS -------------------
samples.fit <- samples.hbef.neon[[curr.set]]
samples.pred <- samples.hbef.neon.bbs[[pred.set]]
elev.fit <- c(hbef.elev[-hbef.x.indices[[curr.set]]], 
	      neon.elev[-neon.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(hbef.elev[hbef.x.indices[[curr.set]]], 
	       bbs.elev[bbs.x.indices[[curr.set]]], 
	       neon.elev[neon.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(hbef.for[-hbef.x.indices[[curr.set]]], 
	      neon.for[-neon.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(hbef.for[hbef.x.indices[[curr.set]]], 
	       bbs.for[bbs.x.indices[[curr.set]]], 
	       neon.for[neon.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.years <- 9
I <- 12
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.hbef <- length(hbef.x.indices[[curr.set]])
J.neon <- length(neon.x.indices[[curr.set]])
J.bbs <- length(bbs.x.indices[[curr.set]])
elev.hbef.pred <- elev.pred[1:J.hbef]
elev.neon.pred <- elev.pred[1:J.neon + J.hbef]
elev.bbs.pred <- elev.pred[1:J.bbs + J.neon + J.hbef]
for.hbef.pred <- for.pred[1:J.hbef]
for.neon.pred <- for.pred[1:J.neon + J.hbef]
for.bbs.pred <- for.pred[1:J.bbs + J.neon + J.hbef]

hbef.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.hbef.pred, 
		    for.pred = for.hbef.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'hbef')
elpd.hbef.vals[, 4, 4] <- apply(log(hbef.elpd.vals.1), 1, sum)

neon.elpd.vals.1 <- elpd.neon.pred(samples.fit = samples.fit, 
		                   samples.pred = samples.pred, 
		                   elev.pred = elev.neon.pred, 
		                   for.pred = for.neon.pred, 
		                   my.iter = my.iter, 
		                   n.iter = n.iter, 
		                   I = I, I.full = I, n.years = n.years)
elpd.neon.vals[sp.neon.indx, 4, 4] <- apply(log(neon.elpd.vals.1), 1, sum)

bbs.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.bbs.pred, 
		    for.pred = for.bbs.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'bbs')

elpd.bbs.vals[, 4, 4] <- apply(log(bbs.elpd.vals.1), 1, sum)

# HBEF + BBS as predicting model -----------------------------------
# HBEF ---------------------------------
samples.fit <- samples.hbef.bbs[[curr.set]]
samples.pred <- samples.hbef[[pred.set]]
elev.fit <- c(hbef.elev[-hbef.x.indices[[curr.set]]], 
	      bbs.elev[-bbs.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(hbef.elev[hbef.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(hbef.for[-hbef.x.indices[[curr.set]]], 
	      bbs.for[-bbs.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(hbef.for[hbef.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.hbef <- length(hbef.x.indices[[curr.set]])
elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.pred, 
		    for.pred = for.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'hbef')
elpd.hbef.vals[, 1, 5] <- apply(log(elpd.vals.1), 1, sum)
# NEON ---------------------------------
samples.fit <- samples.hbef.bbs[[curr.set]]
samples.pred <- samples.neon[[pred.set]]
elev.fit <- c(hbef.elev[-hbef.x.indices[[curr.set]]], 
	      bbs.elev[-bbs.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(neon.elev[neon.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(hbef.for[-hbef.x.indices[[curr.set]]], 
	      bbs.for[-bbs.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(neon.for[neon.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.neon <- length(neon.x.indices[[curr.set]])
I.neon <- 10
n.years.neon <- 4
elpd.vals.1 <- elpd.neon.pred(samples.fit = samples.fit, 
		              samples.pred = samples.pred, 
		              elev.pred = elev.pred, 
		              for.pred = for.pred, 
		              my.iter = my.iter, 
		              n.iter = n.iter, 
		              I = I.neon, I.full = I, n.years = n.years)
elpd.neon.vals[sp.neon.indx, 1, 5] <- apply((log(elpd.vals.1)), 1, sum)
# BBS ---------------------------------
samples.fit <- samples.hbef.bbs[[curr.set]]
samples.pred <- samples.bbs[[pred.set]]
elev.fit <- c(hbef.elev[-hbef.x.indices[[curr.set]]], 
	      bbs.elev[-bbs.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(bbs.elev[bbs.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(hbef.for[-hbef.x.indices[[curr.set]]], 
	      bbs.for[-bbs.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(bbs.for[bbs.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.bbs <- length(bbs.x.indices[[curr.set]])
elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.pred, 
		    for.pred = for.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'bbs')
elpd.bbs.vals[, 1, 5] <- apply(log(elpd.vals.1), 1, sum)
# HBEF + NEON -------------------------
samples.fit <- samples.hbef.bbs[[curr.set]]
samples.pred <- samples.hbef.neon[[pred.set]]
elev.fit <- c(hbef.elev[-hbef.x.indices[[curr.set]]], 
	      bbs.elev[-bbs.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(hbef.elev[hbef.x.indices[[curr.set]]], 
	       neon.elev[neon.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(hbef.for[-hbef.x.indices[[curr.set]]], 
	      bbs.for[-bbs.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(hbef.for[hbef.x.indices[[curr.set]]], 
	       neon.for[neon.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.years <- 9
I <- 12
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)

J.hbef <- length(hbef.x.indices[[curr.set]])
J.neon <- length(neon.x.indices[[curr.set]])
elev.hbef.pred <- elev.pred[1:J.hbef]
elev.neon.pred <- elev.pred[1:J.neon + J.hbef]
for.hbef.pred <- for.pred[1:J.hbef]
for.neon.pred <- for.pred[1:J.neon + J.hbef]

hbef.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		        samples.pred = samples.pred, 
		        elev.pred = elev.hbef.pred, 
		        for.pred = for.hbef.pred, 
		        my.iter = my.iter, 
		        n.iter = n.iter, 
		        I = I, n.years = n.years, type = 'hbef')
elpd.hbef.vals[, 2, 5] <- apply(log(hbef.elpd.vals.1), 1, sum)


neon.elpd.vals.1 <- elpd.neon.pred(samples.fit = samples.fit, 
		                   samples.pred = samples.pred, 
		                   elev.pred = elev.neon.pred, 
		                   for.pred = for.neon.pred, 
		                   my.iter = my.iter, 
		                   n.iter = n.iter, 
		                   I = I, I.full = I, n.years = n.years)

elpd.neon.vals[sp.neon.indx, 2, 5] <- apply(log(neon.elpd.vals.1), 1, sum)
# HBEF + BBS --------------------------
samples.fit <- samples.hbef.bbs[[curr.set]]
samples.pred <- samples.hbef.bbs[[pred.set]]
elev.fit <- c(hbef.elev[-hbef.x.indices[[curr.set]]], 
	      bbs.elev[-bbs.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(hbef.elev[hbef.x.indices[[curr.set]]], 
	       bbs.elev[bbs.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(hbef.for[-hbef.x.indices[[curr.set]]], 
	      bbs.for[-bbs.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(hbef.for[hbef.x.indices[[curr.set]]], 
	       bbs.for[bbs.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.years <- 9
I <- 12
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.hbef <- length(hbef.x.indices[[curr.set]])
J.bbs <- length(bbs.x.indices[[curr.set]])
elev.hbef.pred <- elev.pred[1:J.hbef]
elev.bbs.pred <- elev.pred[1:J.bbs + J.hbef]
for.hbef.pred <- for.pred[1:J.hbef]
for.bbs.pred <- for.pred[1:J.bbs + J.hbef]

hbef.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.hbef.pred, 
		    for.pred = for.hbef.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'hbef')
elpd.hbef.vals[, 3, 5] <- apply(log(hbef.elpd.vals.1), 1, sum)


bbs.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.bbs.pred, 
		    for.pred = for.bbs.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'bbs')
elpd.bbs.vals[, 2, 5] <- apply(log(bbs.elpd.vals.1), 1, sum)
# NEON + BBS --------------------------
samples.fit <- samples.hbef.bbs[[curr.set]]
samples.pred <- samples.neon.bbs[[pred.set]]
elev.fit <- c(hbef.elev[-hbef.x.indices[[curr.set]]], 
	      bbs.elev[-bbs.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(neon.elev[neon.x.indices[[curr.set]]], 
	       bbs.elev[bbs.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(hbef.for[-hbef.x.indices[[curr.set]]], 
	      bbs.for[-bbs.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(neon.for[neon.x.indices[[curr.set]]], 
	       bbs.for[bbs.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.years <- 9
I <- 12
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.neon <- length(neon.x.indices[[curr.set]])
J.bbs <- length(bbs.x.indices[[curr.set]])
elev.neon.pred <- elev.pred[1:J.neon]
elev.bbs.pred <- elev.pred[1:J.bbs + J.neon]
for.neon.pred <- for.pred[1:J.neon]
for.bbs.pred <- for.pred[1:J.bbs + J.neon]

neon.elpd.vals.1 <- elpd.neon.pred(samples.fit = samples.fit, 
		                   samples.pred = samples.pred, 
		                   elev.pred = elev.neon.pred, 
		                   for.pred = for.neon.pred, 
		                   my.iter = my.iter, 
		                   n.iter = n.iter, 
		                   I = I, I.full = I, n.years = n.years)
elpd.neon.vals[sp.neon.indx, 3, 5] <- apply(log(neon.elpd.vals.1), 1, sum)

bbs.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.bbs.pred, 
		    for.pred = for.bbs.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'bbs')

elpd.bbs.vals[, 3, 5] <- apply(log(bbs.elpd.vals.1), 1, sum)
# HBEF + NEON + BBS -------------------
samples.fit <- samples.hbef.bbs[[curr.set]]
samples.pred <- samples.hbef.neon.bbs[[pred.set]]
elev.fit <- c(hbef.elev[-hbef.x.indices[[curr.set]]], 
	      bbs.elev[-bbs.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(hbef.elev[hbef.x.indices[[curr.set]]], 
	       bbs.elev[bbs.x.indices[[curr.set]]], 
	       neon.elev[neon.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(hbef.for[-hbef.x.indices[[curr.set]]], 
	      bbs.for[-bbs.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(hbef.for[hbef.x.indices[[curr.set]]], 
	       bbs.for[bbs.x.indices[[curr.set]]], 
	       neon.for[neon.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.years <- 9
I <- 12
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.hbef <- length(hbef.x.indices[[curr.set]])
J.neon <- length(neon.x.indices[[curr.set]])
J.bbs <- length(bbs.x.indices[[curr.set]])
elev.hbef.pred <- elev.pred[1:J.hbef]
elev.neon.pred <- elev.pred[1:J.neon + J.hbef]
elev.bbs.pred <- elev.pred[1:J.bbs + J.neon + J.hbef]
for.hbef.pred <- for.pred[1:J.hbef]
for.neon.pred <- for.pred[1:J.neon + J.hbef]
for.bbs.pred <- for.pred[1:J.bbs + J.neon + J.hbef]

hbef.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.hbef.pred, 
		    for.pred = for.hbef.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'hbef')
elpd.hbef.vals[, 4, 5] <- apply(log(hbef.elpd.vals.1), 1, sum)

neon.elpd.vals.1 <- elpd.neon.pred(samples.fit = samples.fit, 
		                   samples.pred = samples.pred, 
		                   elev.pred = elev.neon.pred, 
		                   for.pred = for.neon.pred, 
		                   my.iter = my.iter, 
		                   n.iter = n.iter, 
		                   I = I, I.full = I, n.years = n.years)
elpd.neon.vals[sp.neon.indx, 4, 5] <- apply(log(neon.elpd.vals.1), 1, sum)

bbs.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.bbs.pred, 
		    for.pred = for.bbs.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'bbs')

elpd.bbs.vals[, 4, 5] <- apply(log(bbs.elpd.vals.1), 1, sum)

# NEON + BBS as predicting model -----------------------------------
# HBEF ---------------------------------
samples.fit <- samples.neon.bbs[[curr.set]]
samples.pred <- samples.hbef[[pred.set]]
elev.fit <- c(bbs.elev[-bbs.x.indices[[curr.set]]],
	      neon.elev[-neon.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(hbef.elev[hbef.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(bbs.for[-bbs.x.indices[[curr.set]]],
	     neon.for[-neon.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(hbef.for[hbef.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.hbef <- length(hbef.x.indices[[curr.set]])
elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.pred, 
		    for.pred = for.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'hbef')
elpd.hbef.vals[, 1, 6] <- apply(log(elpd.vals.1), 1, sum)
# NEON ---------------------------------
samples.fit <- samples.neon.bbs[[curr.set]]
samples.pred <- samples.neon[[pred.set]]
elev.fit <- c(bbs.elev[-bbs.x.indices[[curr.set]]],
	      neon.elev[-neon.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(neon.elev[neon.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(bbs.for[-bbs.x.indices[[curr.set]]],
	      neon.for[-neon.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(neon.for[neon.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.neon <- length(neon.x.indices[[curr.set]])
I.neon <- 10
n.years.neon <- 4
elpd.vals.1 <- elpd.neon.pred(samples.fit = samples.fit, 
		              samples.pred = samples.pred, 
		              elev.pred = elev.pred, 
		              for.pred = for.pred, 
		              my.iter = my.iter, 
		              n.iter = n.iter, 
		              I = I.neon, I.full = I, n.years = n.years)
elpd.neon.vals[sp.neon.indx, 1, 6] <- apply((log(elpd.vals.1)), 1, sum)
# BBS ---------------------------------
samples.fit <- samples.neon.bbs[[curr.set]]
samples.pred <- samples.bbs[[pred.set]]
elev.fit <- c(bbs.elev[-bbs.x.indices[[curr.set]]],
	      neon.elev[-neon.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(bbs.elev[bbs.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(bbs.for[-bbs.x.indices[[curr.set]]],
	      neon.for[-neon.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(bbs.for[bbs.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.bbs <- length(bbs.x.indices[[curr.set]])
elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.pred, 
		    for.pred = for.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'bbs')
elpd.bbs.vals[, 1, 6] <- apply(log(elpd.vals.1), 1, sum)
# HBEF + NEON -------------------------
samples.fit <- samples.neon.bbs[[curr.set]]
samples.pred <- samples.hbef.neon[[pred.set]]
elev.fit <- c(bbs.elev[-bbs.x.indices[[curr.set]]],
	      neon.elev[-neon.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(hbef.elev[hbef.x.indices[[curr.set]]],
	       neon.elev[neon.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(bbs.for[-bbs.x.indices[[curr.set]]],
	      neon.for[-neon.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(hbef.for[hbef.x.indices[[curr.set]]],
	       neon.for[neon.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.years <- 9
I <- 12
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)

J.hbef <- length(hbef.x.indices[[curr.set]])
J.neon <- length(neon.x.indices[[curr.set]])
elev.hbef.pred <- elev.pred[1:J.hbef]
elev.neon.pred <- elev.pred[1:J.neon + J.hbef]
for.hbef.pred <- for.pred[1:J.hbef]
for.neon.pred <- for.pred[1:J.neon + J.hbef]

hbef.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		        samples.pred = samples.pred, 
		        elev.pred = elev.hbef.pred, 
		        for.pred = for.hbef.pred, 
		        my.iter = my.iter, 
		        n.iter = n.iter, 
		        I = I, n.years = n.years, type = 'hbef')
elpd.hbef.vals[, 2, 6] <- apply(log(hbef.elpd.vals.1), 1, sum)


neon.elpd.vals.1 <- elpd.neon.pred(samples.fit = samples.fit, 
		                   samples.pred = samples.pred, 
		                   elev.pred = elev.neon.pred, 
		                   for.pred = for.neon.pred, 
		                   my.iter = my.iter, 
		                   n.iter = n.iter, 
		                   I = I, I.full = I, n.years = n.years)

elpd.neon.vals[sp.neon.indx, 2, 6] <- apply(log(neon.elpd.vals.1), 1, sum)
# HBEF + BBS --------------------------
samples.fit <- samples.neon.bbs[[curr.set]]
samples.pred <- samples.hbef.bbs[[pred.set]]
elev.fit <- c(bbs.elev[-bbs.x.indices[[curr.set]]],
	      neon.elev[-neon.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(hbef.elev[hbef.x.indices[[curr.set]]],
	       bbs.elev[bbs.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(bbs.for[-bbs.x.indices[[curr.set]]],
	      neon.for[-neon.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(hbef.for[hbef.x.indices[[curr.set]]],
	       bbs.for[bbs.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.years <- 9
I <- 12
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.hbef <- length(hbef.x.indices[[curr.set]])
J.bbs <- length(bbs.x.indices[[curr.set]])
elev.hbef.pred <- elev.pred[1:J.hbef]
elev.bbs.pred <- elev.pred[1:J.bbs + J.hbef]
for.hbef.pred <- for.pred[1:J.hbef]
for.bbs.pred <- for.pred[1:J.bbs + J.hbef]

hbef.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.hbef.pred, 
		    for.pred = for.hbef.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'hbef')
elpd.hbef.vals[, 3, 6] <- apply(log(hbef.elpd.vals.1), 1, sum)


bbs.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.bbs.pred, 
		    for.pred = for.bbs.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'bbs')
elpd.bbs.vals[, 2, 6] <- apply(log(bbs.elpd.vals.1), 1, sum)
# NEON + BBS --------------------------
samples.fit <- samples.neon.bbs[[curr.set]]
samples.pred <- samples.neon.bbs[[pred.set]]
elev.fit <- c(bbs.elev[-bbs.x.indices[[curr.set]]],
	      neon.elev[-neon.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(neon.elev[neon.x.indices[[curr.set]]],
	       bbs.elev[bbs.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(bbs.for[-bbs.x.indices[[curr.set]]],
	      neon.for[-neon.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(neon.for[neon.x.indices[[curr.set]]],
	       bbs.for[bbs.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.years <- 9
I <- 12
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.neon <- length(neon.x.indices[[curr.set]])
J.bbs <- length(bbs.x.indices[[curr.set]])
elev.neon.pred <- elev.pred[1:J.neon]
elev.bbs.pred <- elev.pred[1:J.bbs + J.neon]
for.neon.pred <- for.pred[1:J.neon]
for.bbs.pred <- for.pred[1:J.bbs + J.neon]

neon.elpd.vals.1 <- elpd.neon.pred(samples.fit = samples.fit, 
		                   samples.pred = samples.pred, 
		                   elev.pred = elev.neon.pred, 
		                   for.pred = for.neon.pred, 
		                   my.iter = my.iter, 
		                   n.iter = n.iter, 
		                   I = I, I.full = I, n.years = n.years)
elpd.neon.vals[sp.neon.indx, 3, 6] <- apply(log(neon.elpd.vals.1), 1, sum)

bbs.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.bbs.pred, 
		    for.pred = for.bbs.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'bbs')

elpd.bbs.vals[, 3, 6] <- apply(log(bbs.elpd.vals.1), 1, sum)
# HBEF + NEON + BBS -------------------
samples.fit <- samples.neon.bbs[[curr.set]]
samples.pred <- samples.hbef.neon.bbs[[pred.set]]
elev.fit <- c(bbs.elev[-bbs.x.indices[[curr.set]]],
	      neon.elev[-neon.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(hbef.elev[hbef.x.indices[[curr.set]]],
	       bbs.elev[bbs.x.indices[[curr.set]]],
	       neon.elev[neon.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(bbs.for[-bbs.x.indices[[curr.set]]],
	      neon.for[-neon.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(hbef.for[hbef.x.indices[[curr.set]]],
	       bbs.for[bbs.x.indices[[curr.set]]],
	       neon.for[neon.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.years <- 9
I <- 12
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.hbef <- length(hbef.x.indices[[curr.set]])
J.neon <- length(neon.x.indices[[curr.set]])
J.bbs <- length(bbs.x.indices[[curr.set]])
elev.hbef.pred <- elev.pred[1:J.hbef]
elev.neon.pred <- elev.pred[1:J.neon + J.hbef]
elev.bbs.pred <- elev.pred[1:J.bbs + J.neon + J.hbef]
for.hbef.pred <- for.pred[1:J.hbef]
for.neon.pred <- for.pred[1:J.neon + J.hbef]
for.bbs.pred <- for.pred[1:J.bbs + J.neon + J.hbef]

hbef.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.hbef.pred, 
		    for.pred = for.hbef.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'hbef')
elpd.hbef.vals[, 4, 6] <- apply(log(hbef.elpd.vals.1), 1, sum)

neon.elpd.vals.1 <- elpd.neon.pred(samples.fit = samples.fit, 
		                   samples.pred = samples.pred, 
		                   elev.pred = elev.neon.pred, 
		                   for.pred = for.neon.pred, 
		                   my.iter = my.iter, 
		                   n.iter = n.iter, 
		                   I = I, I.full = I, n.years = n.years)
elpd.neon.vals[sp.neon.indx, 4, 6] <- apply(log(neon.elpd.vals.1), 1, sum)

bbs.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.bbs.pred, 
		    for.pred = for.bbs.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'bbs')

elpd.bbs.vals[, 4, 6] <- apply(log(bbs.elpd.vals.1), 1, sum)

# HBEF + NEON + BBS as predicting model -----------------------------------
# HBEF ---------------------------------
samples.fit <- samples.hbef.neon.bbs[[curr.set]]
samples.pred <- samples.hbef[[pred.set]]
elev.fit <- c(hbef.elev[-hbef.x.indices[[curr.set]]], 
	      bbs.elev[-bbs.x.indices[[curr.set]]], 
	      neon.elev[-neon.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(hbef.elev[hbef.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(hbef.for[-hbef.x.indices[[curr.set]]], 
	      bbs.for[-bbs.x.indices[[curr.set]]], 
	      neon.for[-neon.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(hbef.for[hbef.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.hbef <- length(hbef.x.indices[[curr.set]])
elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.pred, 
		    for.pred = for.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'hbef')
elpd.hbef.vals[, 1, 7] <- apply(log(elpd.vals.1), 1, sum)
# NEON ---------------------------------
samples.fit <- samples.hbef.neon.bbs[[curr.set]]
samples.pred <- samples.neon[[pred.set]]
elev.fit <- c(hbef.elev[-hbef.x.indices[[curr.set]]], 
	      bbs.elev[-bbs.x.indices[[curr.set]]], 
	      neon.elev[-neon.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(neon.elev[neon.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(hbef.for[-hbef.x.indices[[curr.set]]], 
	      bbs.for[-bbs.x.indices[[curr.set]]], 
	      neon.for[-neon.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(neon.for[neon.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.neon <- length(neon.x.indices[[curr.set]])
I.neon <- 10
n.years.neon <- 4
elpd.vals.1 <- elpd.neon.pred(samples.fit = samples.fit, 
		              samples.pred = samples.pred, 
		              elev.pred = elev.pred, 
		              for.pred = for.pred, 
		              my.iter = my.iter, 
		              n.iter = n.iter, 
		              I = I.neon, I.full = I, n.years = n.years)
elpd.neon.vals[sp.neon.indx, 1, 7] <- apply((log(elpd.vals.1)), 1, sum)
# BBS ---------------------------------
samples.fit <- samples.hbef.neon.bbs[[curr.set]]
samples.pred <- samples.bbs[[pred.set]]
elev.fit <- c(hbef.elev[-hbef.x.indices[[curr.set]]], 
	      bbs.elev[-bbs.x.indices[[curr.set]]], 
	      neon.elev[-neon.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(bbs.elev[bbs.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(hbef.for[-hbef.x.indices[[curr.set]]], 
	      bbs.for[-bbs.x.indices[[curr.set]]], 
	      neon.for[-neon.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(bbs.for[bbs.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.bbs <- length(bbs.x.indices[[curr.set]])
elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.pred, 
		    for.pred = for.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'bbs')
elpd.bbs.vals[, 1, 7] <- apply(log(elpd.vals.1), 1, sum)
# HBEF + NEON -------------------------
samples.fit <- samples.hbef.neon.bbs[[curr.set]]
samples.pred <- samples.hbef.neon[[pred.set]]
elev.fit <- c(hbef.elev[-hbef.x.indices[[curr.set]]], 
	      bbs.elev[-bbs.x.indices[[curr.set]]], 
	      neon.elev[-neon.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(hbef.elev[hbef.x.indices[[curr.set]]], 
	       neon.elev[neon.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(hbef.for[-hbef.x.indices[[curr.set]]], 
	      bbs.for[-bbs.x.indices[[curr.set]]], 
	      neon.for[-neon.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(hbef.for[hbef.x.indices[[curr.set]]], 
	       neon.for[neon.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.years <- 9
I <- 12
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.hbef <- length(hbef.x.indices[[curr.set]])
J.neon <- length(neon.x.indices[[curr.set]])
elev.hbef.pred <- elev.pred[1:J.hbef]
elev.neon.pred <- elev.pred[1:J.neon + J.hbef]
for.hbef.pred <- for.pred[1:J.hbef]
for.neon.pred <- for.pred[1:J.neon + J.hbef]

hbef.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		        samples.pred = samples.pred, 
		        elev.pred = elev.hbef.pred, 
		        for.pred = for.hbef.pred, 
		        my.iter = my.iter, 
		        n.iter = n.iter, 
		        I = I, n.years = n.years, type = 'hbef')
elpd.hbef.vals[, 2, 7] <- apply(log(hbef.elpd.vals.1), 1, sum)


neon.elpd.vals.1 <- elpd.neon.pred(samples.fit = samples.fit, 
		                   samples.pred = samples.pred, 
		                   elev.pred = elev.neon.pred, 
		                   for.pred = for.neon.pred, 
		                   my.iter = my.iter, 
		                   n.iter = n.iter, 
		                   I = I, I.full = I, n.years = n.years)

elpd.neon.vals[sp.neon.indx, 2, 7] <- apply(log(neon.elpd.vals.1), 1, sum)
# HBEF + BBS --------------------------
samples.fit <- samples.hbef.neon.bbs[[curr.set]]
samples.pred <- samples.hbef.bbs[[pred.set]]
elev.fit <- c(hbef.elev[-hbef.x.indices[[curr.set]]], 
	      bbs.elev[-bbs.x.indices[[curr.set]]], 
	      neon.elev[-neon.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(hbef.elev[hbef.x.indices[[curr.set]]], 
	       bbs.elev[bbs.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(hbef.for[-hbef.x.indices[[curr.set]]], 
	      bbs.for[-bbs.x.indices[[curr.set]]], 
	      neon.for[-neon.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(hbef.for[hbef.x.indices[[curr.set]]], 
	       bbs.for[bbs.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.years <- 9
I <- 12
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.hbef <- length(hbef.x.indices[[curr.set]])
J.bbs <- length(bbs.x.indices[[curr.set]])
elev.hbef.pred <- elev.pred[1:J.hbef]
elev.bbs.pred <- elev.pred[1:J.bbs + J.hbef]
for.hbef.pred <- for.pred[1:J.hbef]
for.bbs.pred <- for.pred[1:J.bbs + J.hbef]

hbef.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.hbef.pred, 
		    for.pred = for.hbef.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'hbef')
elpd.hbef.vals[, 3, 7] <- apply(log(hbef.elpd.vals.1), 1, sum)


bbs.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.bbs.pred, 
		    for.pred = for.bbs.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'bbs')
elpd.bbs.vals[, 2, 7] <- apply(log(bbs.elpd.vals.1), 1, sum)
# NEON + BBS --------------------------
samples.fit <- samples.hbef.neon.bbs[[curr.set]]
samples.pred <- samples.neon.bbs[[pred.set]]
elev.fit <- c(hbef.elev[-hbef.x.indices[[curr.set]]], 
	      bbs.elev[-bbs.x.indices[[curr.set]]], 
	      neon.elev[-neon.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(neon.elev[neon.x.indices[[curr.set]]], 
	       bbs.elev[bbs.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(hbef.for[-hbef.x.indices[[curr.set]]], 
	      bbs.for[-bbs.x.indices[[curr.set]]], 
	      neon.for[-neon.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(neon.for[neon.x.indices[[curr.set]]], 
	       bbs.for[bbs.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.years <- 9
I <- 12
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.neon <- length(neon.x.indices[[curr.set]])
J.bbs <- length(bbs.x.indices[[curr.set]])
elev.neon.pred <- elev.pred[1:J.neon]
elev.bbs.pred <- elev.pred[1:J.bbs + J.neon]
for.neon.pred <- for.pred[1:J.neon]
for.bbs.pred <- for.pred[1:J.bbs + J.neon]

neon.elpd.vals.1 <- elpd.neon.pred(samples.fit = samples.fit, 
		                   samples.pred = samples.pred, 
		                   elev.pred = elev.neon.pred, 
		                   for.pred = for.neon.pred, 
		                   my.iter = my.iter, 
		                   n.iter = n.iter, 
		                   I = I, I.full = I, n.years = n.years)
elpd.neon.vals[sp.neon.indx, 3, 7] <- apply(log(neon.elpd.vals.1), 1, sum)

bbs.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.bbs.pred, 
		    for.pred = for.bbs.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'bbs')

elpd.bbs.vals[, 3, 7] <- apply(log(bbs.elpd.vals.1), 1, sum)
# HBEF + NEON + BBS -------------------
samples.fit <- samples.hbef.neon.bbs[[curr.set]]
samples.pred <- samples.hbef.neon.bbs[[pred.set]]
elev.fit <- c(hbef.elev[-hbef.x.indices[[curr.set]]], 
	      bbs.elev[-bbs.x.indices[[curr.set]]], 
	      neon.elev[-neon.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(hbef.elev[hbef.x.indices[[curr.set]]], 
	       bbs.elev[bbs.x.indices[[curr.set]]], 
	       neon.elev[neon.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(hbef.for[-hbef.x.indices[[curr.set]]], 
	      bbs.for[-bbs.x.indices[[curr.set]]], 
	      neon.for[-neon.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(hbef.for[hbef.x.indices[[curr.set]]], 
	       bbs.for[bbs.x.indices[[curr.set]]], 
	       neon.for[neon.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)
n.years <- 9
I <- 12
n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.hbef <- length(hbef.x.indices[[curr.set]])
J.neon <- length(neon.x.indices[[curr.set]])
J.bbs <- length(bbs.x.indices[[curr.set]])
elev.hbef.pred <- elev.pred[1:J.hbef]
elev.neon.pred <- elev.pred[1:J.neon + J.hbef]
elev.bbs.pred <- elev.pred[1:J.bbs + J.neon + J.hbef]
for.hbef.pred <- for.pred[1:J.hbef]
for.neon.pred <- for.pred[1:J.neon + J.hbef]
for.bbs.pred <- for.pred[1:J.bbs + J.neon + J.hbef]

hbef.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.hbef.pred, 
		    for.pred = for.hbef.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'hbef')
elpd.hbef.vals[, 4, 7] <- apply(log(hbef.elpd.vals.1), 1, sum)

neon.elpd.vals.1 <- elpd.neon.pred(samples.fit = samples.fit, 
		                   samples.pred = samples.pred, 
		                   elev.pred = elev.neon.pred, 
		                   for.pred = for.neon.pred, 
		                   my.iter = my.iter, 
		                   n.iter = n.iter, 
		                   I = I, I.full = I, n.years = n.years)
elpd.neon.vals[sp.neon.indx, 4, 7] <- apply(log(neon.elpd.vals.1), 1, sum)

bbs.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    elev.pred = elev.bbs.pred, 
		    for.pred = for.bbs.pred, 
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, n.years = n.years, type = 'bbs')

elpd.bbs.vals[, 4, 7] <- apply(log(bbs.elpd.vals.1), 1, sum)

save(elpd.hbef.vals, elpd.neon.vals, 
     elpd.bbs.vals, file = paste('results/cross-val/cross-val-results-', 
				 curr.set, '.R', sep = ''))

