# summary.R: file to summarize cross-validation results for White Mountain 
#            National Forest case study. 
# Author: Jeffrey W. Doser (doserjef@msu.edu)
# Citation: 

rm(list = ls())
library(coda)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(ggpubr)

# Read in results ---------------------------------------------------------
load("results/cross-val/cross-val-results-1.R")
elpd.bbs.vals.1 <- elpd.bbs.vals
elpd.hbef.vals.1 <- elpd.hbef.vals
elpd.neon.vals.1 <- elpd.neon.vals
load("results/cross-val/cross-val-results-2.R")
elpd.bbs.vals.2 <- elpd.bbs.vals
elpd.hbef.vals.2 <- elpd.hbef.vals
elpd.neon.vals.2 <- elpd.neon.vals

# Take average of two-fold cross-validation. 
elpd.bbs.vals.icom <- (elpd.bbs.vals.1 + elpd.bbs.vals.2) / 2
elpd.neon.vals.icom <- (elpd.neon.vals.1 + elpd.neon.vals.2) / 2
elpd.hbef.vals.icom <- (elpd.hbef.vals.1 + elpd.hbef.vals.2) / 2

# Species by model results
bbs.sp.model <- apply(elpd.bbs.vals.icom, c(1, 3), mean, na.rm = TRUE)
hbef.sp.model <- apply(elpd.hbef.vals.icom, c(1, 3), mean, na.rm = TRUE)
neon.sp.model <- apply(elpd.neon.vals.icom, c(1, 3), mean, na.rm = TRUE)

# Order of maximum ELPD values by species
# Making them negative to make orders in reverse
bbs.orders <- apply(-bbs.sp.model, 1, order)
hbef.orders <- apply(-hbef.sp.model, 1, order)
neon.orders <- apply(-neon.sp.model, 1, order)
# Compute average rank for a given model and data set
# HBEF
mean(apply(hbef.orders[-7, ], 2, function(a) which(a == 1)))
mean(apply(hbef.orders[-7, ], 2, function(a) which(a == 2)))
mean(apply(hbef.orders[-7, ], 2, function(a) which(a == 3)))
mean(apply(hbef.orders[-7, ], 2, function(a) which(a == 4)))
mean(apply(hbef.orders[-7, ], 2, function(a) which(a == 5)))
mean(apply(hbef.orders[-7, ], 2, function(a) which(a == 6)))
mean(apply(hbef.orders[-7, ], 2, function(a) which(a == 7)))

# BBS
mean(apply(bbs.orders[-7, ], 2, function(a) which(a == 1)))
mean(apply(bbs.orders[-7, ], 2, function(a) which(a == 2)))
mean(apply(bbs.orders[-7, ], 2, function(a) which(a == 3)))
mean(apply(bbs.orders[-7, ], 2, function(a) which(a == 4)))
mean(apply(bbs.orders[-7, ], 2, function(a) which(a == 5)))
mean(apply(bbs.orders[-7, ], 2, function(a) which(a == 6)))
mean(apply(bbs.orders[-7, ], 2, function(a) which(a == 7)))

# NEON
mean(apply(neon.orders, 2, function(a) which(a == 1)))
mean(apply(neon.orders, 2, function(a) which(a == 2)))
mean(apply(neon.orders, 2, function(a) which(a == 3)))
mean(apply(neon.orders, 2, function(a) which(a == 4)))
mean(apply(neon.orders, 2, function(a) which(a == 5)))
mean(apply(neon.orders, 2, function(a) which(a == 6)))
mean(apply(neon.orders, 2, function(a) which(a == 7)))

# All
all.sp.model <- bbs.sp.model + hbef.sp.model
neon.sp.indx <- c(1, 2, 3, 4, 6, 7, 8, 9, 11, 12)
all.sp.model[neon.sp.indx, ] <- all.sp.model[neon.sp.indx, ] + neon.sp.model[neon.sp.indx, ]
all.orders <- apply(-all.sp.model, 1, order)
mean(apply(all.orders[-7, ], 2, function(a) which(a == 1)))
mean(apply(all.orders[-7, ], 2, function(a) which(a == 2)))
mean(apply(all.orders[-7, ], 2, function(a) which(a == 3)))
mean(apply(all.orders[-7, ], 2, function(a) which(a == 4)))
mean(apply(all.orders[-7, ], 2, function(a) which(a == 5)))
mean(apply(all.orders[-7, ], 2, function(a) which(a == 6)))
mean(apply(all.orders[-7, ], 2, function(a) which(a == 7)))

# All species model results
bbs.model <- apply(bbs.sp.model, 2, sum, na.rm = TRUE)
hbef.model <- apply(hbef.sp.model, 2, sum, na.rm = TRUE)
neon.model <- apply(neon.sp.model, 2, sum, na.rm = TRUE)

# All species all data sets
bbs.model + hbef.model + neon.model
order(bbs.model + hbef.model + neon.model)

# Compare with ISDM results -----------------------------------------------
load("results/cross-val/cross-val-idm-results-1.R")
elpd.bbs.idm.vals.1 <- elpd.bbs.vals
elpd.hbef.idm.vals.1 <- elpd.hbef.vals
elpd.neon.idm.vals.1 <- elpd.neon.vals
load("results/cross-val/cross-val-idm-results-2.R")
elpd.bbs.idm.vals.2 <- elpd.bbs.vals
elpd.hbef.idm.vals.2 <- elpd.hbef.vals
elpd.neon.idm.vals.2 <- elpd.neon.vals

# Take average of two-fold cross validation
elpd.bbs.idm.vals <- (elpd.bbs.idm.vals.1 + elpd.bbs.idm.vals.2) / 2
elpd.neon.idm.vals <- (elpd.neon.idm.vals.1 + elpd.neon.idm.vals.2) / 2
elpd.hbef.idm.vals <- (elpd.hbef.idm.vals.1 + elpd.hbef.idm.vals.2) / 2

# HBEF ---------------------------------
elpd.hbef.vals.icom[, 4, 7]
elpd.hbef.idm.vals
elpd.hbef.vals.icom[, 4, 7] > elpd.hbef.idm.vals

# BBS ---------------------------------
elpd.bbs.vals.icom[, 4, 7]
elpd.bbs.idm.vals
elpd.bbs.vals.icom[, 4, 7] > elpd.bbs.idm.vals

# NEON --------------------------------
elpd.neon.vals.icom[, 4, 7]
elpd.neon.idm.vals
elpd.neon.vals.icom[, 4, 7] > elpd.neon.idm.vals

# Overall
# HBEF
sum(elpd.hbef.vals.icom[-5, 4, 7])
sum(elpd.hbef.idm.vals, na.rm = TRUE)
# BBS
sum(elpd.bbs.vals.icom[-5, 4, 7])
sum(elpd.bbs.idm.vals, na.rm = TRUE)
# NEON
sum(elpd.neon.vals.icom[, 4, 7], na.rm = TRUE)
sum(elpd.neon.idm.vals, na.rm = TRUE)
