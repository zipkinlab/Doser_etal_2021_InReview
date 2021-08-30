# main-HBEF.R: code to run community model using only Hubbard Brook data for
#              cross validation assessment. See code/wmnf/main-HBEF.R for more
#              detailed commenting regarding different objects. 
# Author: Jeffrey W. Doser (doserjef@msu.edu)
# Citation: 

rm(list = ls())

library(coda)
library(tidyverse)
library(nimble)
library(coda)
# This index is used to run different hold out sets across multiple cores. By 
# providing the number of the chain, you will create different output 
# files that you can then combine into a single mcmc.list
curr.set <- as.numeric(commandArgs(trailingOnly = TRUE))
if(length(curr.set) == 0) base::stop('Need to tell NIMBLE the hold out number')
# Alternatively, can specify manually. 
# curr.set <- 1

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
source("code/cross-val-icom/nimble-code/icom-HBEF.R")

# Data prep for analysis --------------------------------------------------
# Number of species
I <- dim(hb.dat)[4]
# Number of HBEF sites
J.hbef <- dim(hb.dat)[1]
# Repeat visits at HBEF
K.hbef <- dim(hb.dat)[2]
# Total number of sites
J <- J.hbef
# Number of years
n.years <- dim(hb.dat)[3]
# Specific years
years.real <- as.numeric(attr(hb.day, 'dimnames')[[3]])
# Standardize year covariate
years <- c(scale(years.real))
# Elevation
elev <- c(hb.elev)
elev.real <- elev
elev <- c(scale(elev))
# Percent forest
percent.for <- c(hb.for)
percent.for.real <- percent.for
percent.for <- c(scale(percent.for))
# Standardize hubbard brook covariates
hb.day <- array(scale(hb.day), dim = c(J.hbef, K.hbef, n.years))
hb.tod <- array(scale(hb.tod), dim = c(J.hbef, K.hbef, n.years))
# Reorganize data array
hb.dat <- aperm(hb.dat, c(1, 2, 4, 3))

# Convert all count data to detection/non-detection data. 
y <- ifelse(hb.dat > 0, 1, hb.dat)

# Create Cross Validation Data Sets ---------------------------------------
# Load the indices used to separate data into two different sets. 
# These are produced in code/cross-val-icom/main-HBEF-NEON-BBS.R
load("data/cross-val-indices.R")

y.fit <- y[-hbef.x.indices[[curr.set]], , , ]
y.pred <- y[hbef.x.indices[[curr.set]], , , ]

# Get occupancy covariates ------------
# Elevation
elev.fit <- c(hb.elev[-hbef.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(hb.elev[hbef.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(hb.for[-hbef.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(hb.for[hbef.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)

# Get detection covariates ------------
# HBEF
hbef.day.fit <- hb.day[-hbef.x.indices[[curr.set]], , ]
hbef.tod.fit <- hb.tod[-hbef.x.indices[[curr.set]], , ]
hbef.day.pred <- hb.day[hbef.x.indices[[curr.set]], , ]
hbef.tod.pred <- hb.tod[hbef.x.indices[[curr.set]], , ]

# Get All Data in Long Format ---------------------------------------------
# Modeling the data in this long format drastically speeds up run times
# HBEF ----------------------------------
# Fit 
y.fit.df <- as.data.frame.table(y.fit)
names(y.fit.df) <- c('Site', 'Visit', 'Species', 'Year', 'Count')
y.fit.df <- y.fit.df %>%
  mutate(Site = as.numeric(Site),
	 Visit = as.numeric(Visit),
	 Species = as.numeric(Species),
	 Year = as.numeric(Year))

# Add DAY to data frame
hbef.day.fit.df <- as.data.frame.table(hbef.day.fit)
names(hbef.day.fit.df) <- c('Site', 'Visit', 'Year', 'Day')
hbef.day.fit.df <- hbef.day.fit.df %>%
  mutate(Site = as.numeric(Site),
	 Visit = as.numeric(Visit),
	 Year = as.numeric(Year))
# Join to full data
y.fit.df <- left_join(y.fit.df, hbef.day.fit.df, by = c('Site', 'Visit', 'Year'))
# Add TOD to data frame
hbef.tod.fit.df <- as.data.frame.table(hbef.tod.fit)
names(hbef.tod.fit.df) <- c('Site', 'Visit', 'Year', 'TOD')
hbef.tod.fit.df <- hbef.tod.fit.df %>%
  mutate(Site = as.numeric(Site),
	 Visit = as.numeric(Visit),
	 Year = as.numeric(Year))
# Join to full data
y.fit.df <- left_join(y.fit.df, hbef.tod.fit.df, by = c('Site', 'Visit', 'Year'))
# Remove rows without data
y.fit.df <- y.fit.df %>%
  filter(!is.na(Count))
y.fit.df$pixel <- y.fit.df$Site
y.fit.df <- y.fit.df %>% arrange(Species, Year)
# Determine observed species 
y.fit.sp.obs <- y.fit.df %>%
  group_by(Species) %>%
  summarize(counts = sum(Count)) %>%
  filter(counts > 0) %>%
  pull(Species)
y.fit.df <- y.fit.df %>%
  filter(Species %in% y.fit.sp.obs)
n.sp.hbef <- length(y.fit.sp.obs)
y.fit.df$p.sp <- as.numeric(factor(y.fit.df$Species))
n.vals.hbef.fit <- nrow(y.fit.df)
# Pred
y.pred.df <- as.data.frame.table(y.pred)
names(y.pred.df) <- c('Site', 'Visit', 'Species', 'Year', 'Count')
y.pred.df <- y.pred.df %>%
  mutate(Site = as.numeric(Site),
	 Visit = as.numeric(Visit),
	 Species = as.numeric(Species),
	 Year = as.numeric(Year))

# Add DAY to data frame
hbef.day.pred.df <- as.data.frame.table(hbef.day.pred)
names(hbef.day.pred.df) <- c('Site', 'Visit', 'Year', 'Day')
hbef.day.pred.df <- hbef.day.pred.df %>%
  mutate(Site = as.numeric(Site),
	 Visit = as.numeric(Visit),
	 Year = as.numeric(Year))
# Join to full data
y.pred.df <- left_join(y.pred.df, hbef.day.pred.df, by = c('Site', 'Visit', 'Year'))
# Add TOD to data frame
hbef.tod.pred.df <- as.data.frame.table(hbef.tod.pred)
names(hbef.tod.pred.df) <- c('Site', 'Visit', 'Year', 'TOD')
hbef.tod.pred.df <- hbef.tod.pred.df %>%
  mutate(Site = as.numeric(Site),
	 Visit = as.numeric(Visit),
	 Year = as.numeric(Year))
# Join to full data
y.pred.df <- left_join(y.pred.df, hbef.tod.pred.df, by = c('Site', 'Visit', 'Year'))
# Remove rows without data
y.pred.df <- y.pred.df %>%
  filter(!is.na(Count), !is.na(Day), !is.na(TOD))
n.vals.hbef.pred <- nrow(y.pred.df)
y.pred.df$pixel <- y.pred.df$Site
y.pred.df <- y.pred.df %>% arrange(Species, Year)
y.pred.df <- y.pred.df %>%
  filter(Species %in% y.fit.sp.obs)
y.pred.df$p.sp <- as.numeric(factor(y.pred.df$Species))
n.vals.hbef.pred <- nrow(y.pred.df)

# Redfine all constants ---------------------------------------------------
J.fit <- sum(n_distinct(y.fit.df$Site))
J.pred <- sum(n_distinct(y.pred.df$Site))

# Constants -------------------------------------------------------------
icom.consts <- list(n.years = n.years, I = I, 
		   J = J.fit, n.vals.hbef = n.vals.hbef.fit, 
		   sp.indx.hbef = y.fit.df$Species, year.indx.hbef = y.fit.df$Year, 
		   site.hbef = y.fit.df$Site, years = years, 
		   I.hbef = n.sp.hbef, 
                   sp.indx.hbef.p = y.fit.df$p.sp)

# Data ------------------------------------------------------------------
icom.data <- list(ELEV = elev.fit, FOREST = for.fit,
		 y = y.fit.df$Count, 
		 DAY.hbef = y.fit.df$Day, TOD.hbef = y.fit.df$TOD)
# Initial values --------------------------------------------------------
z.init <- array(1, dim = c(I, J.fit, n.years))
icom.inits <- list(z.hbef = z.init, int.beta.mean = runif(n.years, 0.3, 0.9), 
		  beta.1.mean = rnorm(1), beta.2.mean = rnorm(1), 
		  beta.3.mean = rnorm(1), 
		  phi.mean = rnorm(1),
		  int.alpha.mean = runif(1, 0.3, 0.9), alpha.1.mean = rnorm(1), 
		  alpha.2.mean = rnorm(1), 
		  alpha.3.mean = rnorm(1), 
		  tau.beta.0 = runif(n.years, 0.1, 2), tau.beta.1 = runif(1, 0.1, 2),
		  tau.beta.2 = runif(1, 0.1, 2), tau.beta.3 = runif(1, 0.1, 2),
		  tau.phi = runif(1, 0.1, 2),
		  tau.alpha.0 = runif(1, 0.1, 2), 
		  tau.alpha.1 = runif(1, 0.1, 2), 
		  tau.alpha.2 = runif(1, 0.1, 2), tau.alpha.3 = runif(1, 0.1, 2))
# Create the model ------------------------------------------------------
icom.model <- nimbleModel(code = icom.code, name = 'icom', constants = icom.consts,
		          data = icom.data, inits = icom.inits)

# Configure MCMC --------------------------------------------------------
icom.conf <- configureMCMC(icom.model, monitors = c('int.beta', 'beta.1', 
      					          'beta.2', 'beta.3', 'phi', 
						  'int.beta.mean', 'beta.1.mean', 
					          'beta.2.mean', 'beta.3.mean', 
						  'phi.mean',  
						  'z.hbef'))
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

file.name <- paste('results/icom-HBEF-results-', n.iter, '-iterations-', curr.set, '-holdout-', 
		   date, '.R', sep = '')

save(samples, file = file.name)
