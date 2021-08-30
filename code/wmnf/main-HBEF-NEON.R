# main-HBEF-NEON.R: code to run integrated community model using
#                   HBEF and NEON data. 
# Author: Jeffrey W. Doser
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
#   Rscript code/wmnf/main-HBEF-NEON.R 1 & 
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
source("code/wmnf/nimble-code/icom-HBEF-NEON.R")

# Data prep for analysis --------------------------------------------------
# Number of species
I <- dim(hb.dat)[4]
# Number of DET sites
J.hbef <- dim(hb.dat)[1]
# Repeat visits at DET
K.hbef <- dim(hb.dat)[2]
# Number of NEON sites
J.neon <- dim(v.1)[1]
# Total number of sites
J <- J.hbef + J.neon
# Number of years
n.years <- dim(hb.dat)[3]
# Specific years
years <- as.numeric(attr(hb.day, 'dimnames')[[3]])
# Standardize year covariate
years <- c(scale(years))
# Sites are organized by DET, then NEON
# Elevation
elev <- c(hb.elev, neon.elev)
elev.real <- elev
elev <- c(scale(elev))
# Percent Forest
percent.for <- c(hb.for, neon.for)
percent.for.real <- percent.for
percent.for <- c(scale(percent.for))
# Standardize hubbard brook covariates
hb.day <- array(scale(hb.day), dim = c(J.hbef, K.hbef, n.years))
hb.tod <- array(scale(hb.tod), dim = c(J.hbef, K.hbef, n.years))
# Standardize neon data
n.years.neon <- dim(v.1)[3]
years.neon <- years[6:9]
neon.years <- which(years %in% years.neon)
neon.day <- array(scale(day.neon), dim = c(J.neon, n.years.neon))
neon.hour <- array(scale(hour.neon), dim = c(J.neon, n.years.neon))
# Reorganize arrays
hb.dat <- aperm(hb.dat, c(1, 2, 4, 3))
v.1 <- aperm(v.1, c(1, 4, 2, 3))

# Convert all count data to hbefection/non-hbefection data. 
v.1 <- ifelse(v.1 > 0, 1, v.1)
y <- ifelse(hb.dat > 0, 1, hb.dat)

# Set all subsequent hbefections of a species to 0 if hbefected in 
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
# DET ----------------------------------
y.df <- as.data.frame.table(y)
names(y.df) <- c('Site', 'Visit', 'Species', 'Year', 'Count')
y.df <- y.df %>%
  mutate(Site = as.numeric(Site),
	 Visit = as.numeric(Visit),
	 Species = as.numeric(Species),
	 Year = as.numeric(Year))

# Add DAY to data frame
hb.day.df <- as.data.frame.table(hb.day)
names(hb.day.df) <- c('Site', 'Visit', 'Year', 'Day')
hb.day.df <- hb.day.df %>%
  mutate(Site = as.numeric(Site),
	 Visit = as.numeric(Visit),
	 Year = as.numeric(Year))
# Join to full data
y.df <- left_join(y.df, hb.day.df, by = c('Site', 'Visit', 'Year'))
# Add TOD to data frame
hb.tod.df <- as.data.frame.table(hb.tod)
names(hb.tod.df) <- c('Site', 'Visit', 'Year', 'TOD')
hb.tod.df <- hb.tod.df %>%
  mutate(Site = as.numeric(Site),
	 Visit = as.numeric(Visit),
	 Year = as.numeric(Year))
# Join to full data
y.df <- left_join(y.df, hb.tod.df, by = c('Site', 'Visit', 'Year'))
# Remove rows without data
y.df <- y.df %>%
  filter(!is.na(Count))
n.vals.hbef <- nrow(y.df)
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

# Only estimate variables for species that were observed. 
y.sp.obs <- y.df %>%
  group_by(Species) %>%
  summarize(counts = sum(Count)) %>%
  filter(counts > 0) %>%
  pull(Species)
y.df <- y.df %>%
  filter(Species %in% y.sp.obs)
n.sp.hbef <- length(y.sp.obs)
y.df$p.sp <- as.numeric(factor(y.df$Species))
n.vals.hbef <- nrow(y.df)
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
# Include site identity for each data set to match up in full data set
y.df$pixel <- y.df$Site
v.1.df$pixel <- v.1.df$Site + J.hbef

# Get indices necessary for Bayesian p-value
# low.bp: lower index for each species and year
# high.bp: higher index for each species and year
y.df <- y.df %>% arrange(Species, Year)
low.bp.y <- matrix(NA, nrow = I, ncol = n.years)
high.bp.y <- matrix(NA, nrow = I, ncol = n.years)
for (i in 1:I) {
  for (t in 1:n.years) {
    tmp <- which(y.df$Species == i & y.df$Year == t)
    low.bp.y[i, t] <- tmp[1]
    high.bp.y[i, t] <- tmp[length(tmp)]
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

J.hbef <- n_distinct(y.df$Site)
J.neon <- n_distinct(v.1.df$Site)
J <- sum(n_distinct(y.df$Site),
	 n_distinct(v.1.df$Site))
elev.hbef <- elev[1:J.hbef]
elev.neon <- elev[1:J.neon + J.hbef]
for.hbef <- percent.for[1:J.hbef]
for.neon <- percent.for[1:J.neon + J.hbef]
years.neon <- 6:9
sp.neon <- unique(v.1.df$Species)

# Constants -------------------------------------------------------------
icom.consts <- list(n.years = n.years, n.years.neon = n.years.neon, I = I, 
		   n.vals.hbef = n.vals.hbef,  
		   sp.indx.hbef = y.df$Species, year.indx.hbef = y.df$Year, 
		   site.hbef = y.df$Site, years = years, 
		   sp.indx.neon = v.1.df$Species, 
		   year.indx.neon = v.1.df$Year, site.neon = v.1.df$Site, 
		   year.indx.neon.all = v.1.df$Year.all, 
		   I.hbef = n.sp.hbef, I.neon = n.sp.neon, 
		   sp.indx.hbef.p = y.df$p.sp, sp.indx.neon.p = v.1.df$p.sp, 
		   J.hbef = J.hbef, J.neon = J.neon, years.neon = years.neon, 
		   sp.neon = sp.neon, 
		   low.bp.y = low.bp.y, high.bp.y = high.bp.y, 
		   low.bp.v.1 = low.bp.v.1, high.bp.v.1 = high.bp.v.1, 
		   e = 0.0001)
# Data ------------------------------------------------------------------
icom.data <- list(ELEV.hbef = elev.hbef, FOREST.hbef = for.hbef, 
		 ELEV.neon = elev.neon, FOREST.neon = for.neon,
		 y = y.df$Count, v.1 = v.1.df$Count,
		 DAY.neon = v.1.df$Day, HOUR.neon = v.1.df$TOD,
		 DAY.hbef = y.df$Day, TOD.hbef = y.df$TOD)
# Initial values --------------------------------------------------------
z.hbef.init <- array(1, dim = c(I, J.hbef, n.years))
z.neon.init <- array(1, dim = c(n.sp.neon, J.neon, n.years.neon))
icom.inits <- list(z.hbef = z.hbef.init, 
		  z.neon = z.neon.init, 
		  int.beta.mean = runif(n.years, 0.3, 0.9), 
		  beta.1.mean = rnorm(1), beta.2.mean = rnorm(1), 
		  beta.3.mean = rnorm(1), 
		  phi.mean = rnorm(1),
		  int.alpha.mean = runif(1, 0.3, 0.9), alpha.1.mean = rnorm(1), 
		  alpha.2.mean = rnorm(1), 
		  alpha.3.mean = rnorm(1), int.gamma.1.mean = runif(1, 0.3, 0.8), 
      	          gamma.1.1.mean = rnorm(1), 
      	          gamma.1.2.mean = rnorm(1), 
		  gamma.1.3.mean = rnorm(1), 
		  tau.beta.0 = runif(n.years, 0.1, 2), tau.beta.1 = runif(1, 0.1, 2),
		  tau.beta.2 = runif(1, 0.1, 2), tau.beta.3 = runif(1, 0.1, 2),
		  tau.phi = runif(1, 0.1, 2),
		  tau.alpha.0 = runif(1, 0.1, 2), 
		  tau.alpha.1 = runif(1, 0.1, 2), 
		  tau.alpha.2 = runif(1, 0.1, 2), tau.alpha.3 = runif(1, 0.1, 2), 
		  tau.gamma.1.0 = runif(1, 0.1, 2), 
		  tau.gamma.1.1 = runif(1, 0.1, 2), tau.gamma.1.2 = runif(1, 0.1, 2), 
		  tau.gamma.1.3 = runif(1, 0.1, 2))
# Create the model ------------------------------------------------------
icom.model <- nimbleModel(code = icom.code, name = 'icom', constants = icom.consts,
		          data = icom.data, inits = icom.inits)

# Configure MCMC --------------------------------------------------------
icom.conf <- configureMCMC(icom.model, monitors = c('int.beta', 'beta.1', 
      					          'beta.2', 'beta.3', 'phi', 
						  'int.beta.mean', 'beta.1.mean', 
					          'beta.2.mean', 'beta.3.mean', 
						  'chi.2.v.1', 'chi.2.rep.v.1', 
						  'chi.2.y', 'chi.2.rep.y', 
						  'phi.mean', 'z.hbef', 'psi.hbef', 
						  'z.neon', 'psi.neon'))
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

file.name <- paste('results/icom-HBEF-NEON-results-', n.iter, '-iterations-', chain, '-chain-', 
		   date, '.R', sep = '')

save(samples, file = file.name)

