# main-NEON-BBS.R: code for running integrated community model with
#                  NEON and BBS data. 
# Author: Jeffrey W. Doser (doserjef@msu.edu)
# Citation: 

rm(list = ls())

library(coda)
library(tidyverse)
library(nimble)
library(coda)
# This index is used to run these chains across multiple cores. By 
# providing the number of the chain, you will create different output 
# files that you can then combine into a single mcmc.list
# For example: to run the model using chain 1, run the following on the 
#              command line from the home directory of this project
#   Rscript code/wmnf/main-NEON-BBS.R 1 & 
chain <- as.numeric(commandArgs(trailingOnly = TRUE))
if(length(chain) == 0) base::stop('Need to tell NIMBLE the chain number')

# Load in all the data. 
load("data/final-bird-data.R")
# The loaded data file consists of the following objects: 
#    bbs.elev: elevation for BBS sites
#    bbs.for: percent forest for BBS sites
#    day.bbs: day of survey for BBS data
#    day.neon: day of survey for NEON data
#    hb.covs: covariate information for Hubbard Brook data
#    hb.dat: bird observation data for Hubbard Brook
#    hb.day: day of survey for HBEF data
#    hb.elev: elevation for HBEF sites
#    hb.tod: time of day of HBEF surveys
#    hour.neon: hour survey began for NEON surveys
#    neon.elev: elevation for NEON sites
#    neon.for: forest cover for NEON sites
#    obsv.bbs: observer identity for BBS data
#    v.1: NEON observations
#    v.2: BBS observations
source("code/wmnf/nimble-code/icom-NEON-BBS.R")

# Data prep for analysis --------------------------------------------------
# Number of species
I <- dim(v.2)[1]
# Number of BBS sites
J.bbs <- dim(v.2)[2]
# Number of NEON sites
J.neon <- dim(v.1)[1]
# Total number of sites
J <- J.neon + J.bbs
# Number of years
n.years <- dim(v.2)[3]
# Specific years
years <- as.numeric(attr(hb.day, 'dimnames')[[3]])
# Standardize year covariate
years <- c(scale(years))
# Sites are organized by NEON, then BBS
# Elevation
elev <- c(neon.elev, bbs.elev)
elev.real <- elev
elev <- c(scale(elev))
# Percent forest
percent.for <- c(neon.for, bbs.for)
percent.for.real <- percent.for
percent.for <- c(scale(percent.for))
# Standardize neon data
n.years.neon <- dim(v.1)[3]
years.neon <- years[6:9]
neon.years <- which(years %in% years.neon)
neon.day <- array(scale(day.neon), dim = c(J.neon, n.years.neon))
neon.hour <- array(scale(hour.neon), dim = c(J.neon, n.years.neon))
# Standardize BBS Data
day.bbs <- array(scale(day.bbs), dim = c(J.bbs, n.years))
n.obsv.bbs <- length(unique(c(obsv.bbs))) 
obsv.bbs <- array(as.numeric(factor(obsv.bbs)), dim = c(J.bbs, n.years))
# Reorganize arrays
v.1 <- aperm(v.1, c(1, 4, 2, 3))
v.2 <- aperm(v.2, c(2, 1, 3))

# Convert all count data to detection/non-detection data. 
v.1 <- ifelse(v.1 > 0, 1, v.1)
v.2 <- ifelse(v.2 > 0, 1, v.2)

# Set all subsequent detections of a species to 0 if detected in 
# previous time bin. 
for (j in 1:J.neon) {
  for (i in 1:I) {
    for (t in 1:n.years.neon) {
      if (!is.na(v.1[j, i, 1, t])) {
        if (v.1[j, i, 1, t] == 1) {
          v.1[j, i, 2:3, t] <- NA
        } else if (v.1[j, i, 2, t] == 1) {
          v.1[j, i, 3, t] <- NA
        }
      }
    } # t
  } # i
} # j

# Get All Data in Long Format ---------------------------------------------
# Modeling the data in this long format drastically speeds up run times
# NEON --------------------------------
v.1.df <- as.data.frame.table(v.1)
names(v.1.df) <- c('Site', 'Species', 'Visit', 'Year', 'Count')
v.1.df <- v.1.df %>%
  mutate(Site = as.numeric(Site),
	 Species = as.numeric(Species),
	 Visit = as.numeric(Visit),
	 Year = as.numeric(Year))

# Add DAY to data frame
neon.day.df <- as.data.frame.table(neon.day)
names(neon.day.df) <- c('Site', 'Year', 'Day')
neon.day.df <- neon.day.df %>%
  mutate(Site = as.numeric(Site),
	 Year = as.numeric(Year))
# Join to full data
v.1.df <- left_join(v.1.df, neon.day.df, by = c('Site', 'Year'))
# Add TOD to data frame
neon.hour.df <- as.data.frame.table(neon.hour)
names(neon.hour.df) <- c('Site', 'Year', 'TOD')
neon.hour.df <- neon.hour.df %>%
  mutate(Site = as.numeric(Site),
	 Year = as.numeric(Year))
# Join to full data
v.1.df <- left_join(v.1.df, neon.hour.df, by = c('Site', 'Year'))
# Remove rows without data
v.1.df <- v.1.df %>%
  filter(!is.na(Count))
n.vals.neon <- nrow(v.1.df)
v.1.df$Year.all <- neon.years[v.1.df$Year]
# BBS ---------------------------------
v.2.df <- as.data.frame.table(v.2)
names(v.2.df) <- c('Site', 'Species', 'Year', 'Count')
v.2.df <- v.2.df %>%
  mutate(Site = as.numeric(Site),
	 Species = as.numeric(Species),
	 Year = as.numeric(Year))

# Add DAY to data frame
day.bbs.df <- as.data.frame.table(day.bbs)
names(day.bbs.df) <- c('Site', 'Year', 'Day')
day.bbs.df <- day.bbs.df %>%
  mutate(Site = as.numeric(Site),
	 Year = as.numeric(Year))
# Join to full data
v.2.df <- left_join(v.2.df, day.bbs.df, by = c('Site', 'Year'))
# Add observer identity to data frame
obsv.bbs.df <- as.data.frame.table(obsv.bbs)
names(obsv.bbs.df) <- c('Site', 'Year', 'obsv')
obsv.bbs.df <- obsv.bbs.df %>%
  mutate(Site = as.numeric(Site),
	 Year = as.numeric(Year))
# Join to full data
v.2.df <- left_join(v.2.df, obsv.bbs.df, by = c('Site', 'Year'))
# Remove rows without data
v.2.df <- v.2.df %>%
  filter(!is.na(Count))

# Only estimate variables for species that were observed. 
v.1.sp.obs <- v.1.df %>%
  group_by(Species) %>%
  summarize(counts = sum(Count)) %>%
  filter(counts > 0) %>%
  pull(Species) 
v.1.df <- v.1.df %>%
  filter(Species %in% v.1.sp.obs)
n.sp.neon <- length(v.1.sp.obs)
v.1.df$p.sp <- as.numeric(factor(v.1.df$Species))
n.vals.neon <- nrow(v.1.df)
v.2.sp.obs <- v.2.df %>%
  group_by(Species) %>%
  summarize(counts = sum(Count)) %>%
  filter(counts > 0) %>%
  pull(Species)
v.2.df <- v.2.df %>%
  filter(Species %in% v.2.sp.obs)
n.sp.bbs <- length(v.2.sp.obs)
v.2.df$p.sp <- as.numeric(factor(v.2.df$Species))
n.obsv.bbs <- n_distinct(v.2.df$obsv)
n.vals.bbs <- nrow(v.2.df)

# Get indices necessary for Bayesian p-value
# low.bp: lower index for each species and year
# high.bp: higher index for each species and year
v.2.df <- v.2.df %>% arrange(Species, Year)
low.bp.v.2 <- matrix(NA, nrow = I, ncol = n.years)
high.bp.v.2 <- matrix(NA, nrow = I, ncol = n.years)
for (i in 1:I) {
  for (t in 1:n.years) {
    tmp <- which(v.2.df$p.sp == i & v.2.df$Year == t)
    low.bp.v.2[i, t] <- tmp[1]
    high.bp.v.2[i, t] <- tmp[length(tmp)]
  } # t
} # i

# NEON is observed over less years and for less species, so need to 
# make this one a bit special. 
v.1.df <- v.1.df %>% arrange(Species, Year)
low.bp.v.1 <- matrix(NA, nrow = n.sp.neon, ncol = n.years.neon)
high.bp.v.1 <- matrix(NA, nrow = n.sp.neon, ncol = n.years.neon)
for (i in 1:n.sp.neon) {
  for (t in 1:n.years.neon) {
    tmp <- which(v.1.df$p.sp == i & v.1.df$Year == t)
    low.bp.v.1[i, t] <- tmp[1]
    high.bp.v.1[i, t] <- tmp[length(tmp)]
  } # t
} # i

elev.neon <- elev[1:J.neon]
elev.bbs <- elev[1:J.bbs + J.neon]
for.neon <- percent.for[1:J.neon]
for.bbs <- percent.for[1:J.bbs + J.neon]
years.neon <- 6:9
sp.neon <- unique(v.1.df$Species)

# Constants -------------------------------------------------------------
icom.consts <- list(n.years = n.years, n.years.neon = n.years.neon, I = I, 
		   n.vals.neon = n.vals.neon, n.vals.bbs = n.vals.bbs,
		   years = years, 
		   sp.indx.bbs = v.2.df$Species, 
		   year.indx.bbs = v.2.df$Year, obsv.bbs = v.2.df$obsv, 
		   site.bbs = v.2.df$Site, sp.indx.neon = v.1.df$Species, 
		   year.indx.neon = v.1.df$Year, site.neon = v.1.df$Site, 
		   year.indx.neon.all = v.1.df$Year.all, n.obsv.bbs = n.obsv.bbs, 
		   I.bbs = n.sp.bbs, I.neon = n.sp.neon,
		   sp.indx.neon.p = v.1.df$p.sp, sp.indx.bbs.p = v.2.df$p.sp,
		   J.neon = J.neon, J.bbs = J.bbs, years.neon = years.neon, 
		   sp.neon = sp.neon, 
		   low.bp.v.1 = low.bp.v.1, high.bp.v.1 = high.bp.v.1, 
		   low.bp.v.2 = low.bp.v.2, high.bp.v.2 = high.bp.v.2, 
		   e = 0.0001)

# Data ------------------------------------------------------------------
icom.data <- list(ELEV.bbs = elev.bbs, FOREST.bbs = for.bbs, 
		 ELEV.neon = elev.neon, FOREST.neon = for.neon, 
		 v.1 = v.1.df$Count,
		 DAY.neon = v.1.df$Day, HOUR.neon = v.1.df$TOD,
      	         DAY.bbs = v.2.df$Day, v.2 = v.2.df$Count)
# Initial values --------------------------------------------------------
z.neon.init <- array(1, dim = c(n.sp.neon, J.neon, n.years.neon))
z.bbs.init <- array(1, dim = c(I, J.bbs, n.years))
icom.inits <- list(z.bbs = z.bbs.init, 
		  z.neon = z.neon.init, 
		  int.beta.mean = runif(n.years, 0.3, 0.9), 
		  beta.1.mean = rnorm(1), beta.2.mean = rnorm(1), 
		  beta.3.mean = rnorm(1), 
		  phi.mean = rnorm(1),
		  int.gamma.1.mean = runif(1, 0.3, 0.8), 
      	          int.gamma.2.mean = runif(1, 0.2, 0.8), gamma.1.1.mean = rnorm(1), 
      	          gamma.1.2.mean = rnorm(1), tau.gamma.2.3 = runif(1, 0.1, 2), 
		  gamma.2.1.mean = rnorm(1), gamma.2.2.mean = rnorm(1), 
		  gamma.1.3.mean = rnorm(1), 
		  tau.beta.0 = runif(n.years, 0.1, 2), tau.beta.1 = runif(1, 0.1, 2),
		  tau.beta.2 = runif(1, 0.1, 2), tau.beta.3 = runif(1, 0.1, 2),
		  tau.phi = runif(1, 0.1, 2),
		  tau.gamma.1.0 = runif(1, 0.1, 2), 
		  tau.gamma.1.1 = runif(1, 0.1, 2), tau.gamma.1.2 = runif(1, 0.1, 2), 
		  tau.gamma.1.3 = runif(1, 0.1, 2), 
		  tau.gamma.2.0 = runif(1, 0.1, 2), tau.gamma.2.1 = runif(1, 0.1, 2), 
		  tau.gamma.2.2 = runif(1, 0.1, 2))
# Create the model ------------------------------------------------------
icom.model <- nimbleModel(code = icom.code, name = 'icom', constants = icom.consts,
		          data = icom.data, inits = icom.inits)

# Configure MCMC --------------------------------------------------------
icom.conf <- configureMCMC(icom.model, monitors = c('int.beta', 'beta.1', 
      					          'beta.2', 'beta.3', 'phi', 
						  'int.beta.mean', 'beta.1.mean', 
					          'beta.2.mean', 'beta.3.mean', 
						  'chi.2.v.1', 'chi.2.rep.v.1', 
						  'chi.2.bbs', 'chi.2.rep.bbs', 
						  'phi.mean','z.neon', 'psi.neon', 
						  'z.bbs', 'psi.bbs'))
# Create an MCMC function -------------------------------------------
icom.mcmc <- buildMCMC(icom.conf)
# Compile model
icom.c.model <- compileNimble(icom.model)
icom.c.mcmc <- compileNimble(icom.mcmc, project = icom.model)
# Number of iterations --------------------------------------------------
n.iter <- 350000
n.burn <- 100000
n.thin <- 20
n.chain <- 1
samples <- runMCMC(icom.c.mcmc, niter = n.iter, nburnin = n.burn,
	           thin = n.thin, nchains = n.chain, samplesAsCodaMCMC = TRUE)
date <- Sys.Date()

file.name <- paste('results/icom-NEON-BBS-results-', n.iter, '-iterations-', chain, '-chain-', 
		   date, '.R', sep = '')

save(samples, file = file.name)
