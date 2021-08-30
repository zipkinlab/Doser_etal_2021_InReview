# main-BBS.R: code to run community model using only BBS data
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
#   Rscript code/wmnf/main-BBS.R 1 & 
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
source("code/wmnf/nimble-code/icom-BBS.R")

# Data prep for analysis --------------------------------------------------
# Number of species
I <- dim(v.2)[1]
# Number of BBS sites
J.bbs <- dim(v.2)[2]
# Total number of sites
J <- J.bbs
# Number of years
n.years <- dim(v.2)[3]
# Specific years
years <- as.numeric(attr(hb.day, 'dimnames')[[3]])
# Standardize year covariate
years <- c(scale(years))
# Elevation
elev <- c(bbs.elev)
elev.real <- elev
elev <- c(scale(elev))
# Percent forest
percent.for <- c(bbs.for)
percent.for.real <- percent.for
percent.for <- c(scale(percent.for))
# Standardize BBS Data
day.bbs <- array(scale(day.bbs), dim = c(J.bbs, n.years))
n.obsv.bbs <- length(unique(c(obsv.bbs))) 
obsv.bbs <- array(as.numeric(factor(obsv.bbs)), dim = c(J.bbs, n.years))

# Convert all count data to detection/non-detection data. 
v.2 <- ifelse(v.2 > 0, 1, v.2)
# Reorder v.2
v.2 <- aperm(v.2, c(2, 1, 3))

# Get All Data in Long Format ---------------------------------------------
# Modeling the data in this long format drastically speeds up run times
# BBS ---------------------------------
v.2.bbs.df <- as.data.frame.table(v.2)
names(v.2.bbs.df) <- c('Site', 'Species', 'Year', 'Count')
v.2.bbs.df <- v.2.bbs.df %>%
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
v.2.bbs.df <- left_join(v.2.bbs.df, day.bbs.df, by = c('Site', 'Year'))
# Add observer identity to data frame
obsv.bbs.df <- as.data.frame.table(obsv.bbs)
names(obsv.bbs.df) <- c('Site', 'Year', 'obsv')
obsv.bbs.df <- obsv.bbs.df %>%
  mutate(Site = as.numeric(Site),
	 Year = as.numeric(Year))
# Join to full data
v.2.bbs.df <- left_join(v.2.bbs.df, obsv.bbs.df, by = c('Site', 'Year'))
# Remove rows without data
v.2.bbs.df <- v.2.bbs.df %>%
  filter(!is.na(Count))

# Only estimate variables for species that were observed. 
v.2.sp.obs <- v.2.bbs.df %>%
  group_by(Species) %>%
  summarize(counts = sum(Count)) %>%
  filter(counts > 0) %>%
  pull(Species)
v.2.bbs.df <- v.2.bbs.df %>%
  filter(Species %in% v.2.sp.obs)
n.sp.bbs <- length(v.2.sp.obs)
v.2.bbs.df$p.sp <- as.numeric(factor(v.2.bbs.df$Species))
n.obsv.bbs <- n_distinct(v.2.bbs.df$obsv)
n.vals.bbs <- nrow(v.2.bbs.df)
# Include site identity for each data set to match up in full data set
v.2.bbs.df$pixel <- v.2.bbs.df$Site

# Get indices necessary for Bayesian p-value
# low.bp: lower index for each species and year
# high.bp: higher index for each species and year
v.2.bbs.df <- v.2.bbs.df %>% arrange(Species, Year)
low.bp.v.2 <- matrix(NA, nrow = I, ncol = n.years)
high.bp.v.2 <- matrix(NA, nrow = I, ncol = n.years)
for (i in 1:I) {
  for (t in 1:n.years) {
    tmp <- which(v.2.bbs.df$p.sp == i & v.2.bbs.df$Year == t)
    low.bp.v.2[i, t] <- tmp[1]
    high.bp.v.2[i, t] <- tmp[length(tmp)]
  } # t
} # i

# Constants -------------------------------------------------------------
icom.consts <- list(n.years = n.years, I = I, 
		   J = J, 
		   n.vals.bbs = n.vals.bbs,
		   years = years, 
		   sp.indx.bbs = v.2.bbs.df$p.sp, 
		   year.indx.bbs = v.2.bbs.df$Year, obsv.bbs = v.2.bbs.df$obsv, 
		   site.bbs = v.2.bbs.df$Site, n.obsv.bbs = n.obsv.bbs, 
		   low.bp.v.2 = low.bp.v.2, high.bp.v.2 = high.bp.v.2, 
		   e = 0.0001)
# Data ------------------------------------------------------------------
icom.data <- list(ELEV = elev, FOREST = percent.for, 
      	         DAY.bbs = v.2.bbs.df$Day, v.2 = v.2.bbs.df$Count)
# Initial values --------------------------------------------------------
z.bbs.init <- array(1, dim = c(I, J, n.years))
icom.inits <- list(z.bbs = z.bbs.init, int.beta.mean = runif(n.years, 0.3, 0.9), 
		  beta.1.mean = rnorm(1), beta.2.mean = rnorm(1), 
		  beta.3.mean = rnorm(1), 
		  phi.mean = rnorm(1),
      	          int.gamma.2.mean = runif(1, 0.2, 0.8), tau.gamma.2.3 = runif(1, 0.1, 2), 
		  gamma.2.1.mean = rnorm(1), gamma.2.2.mean = rnorm(1), 
		  tau.beta.0 = runif(n.years, 0.1, 2), tau.beta.1 = runif(1, 0.1, 2),
		  tau.beta.2 = runif(1, 0.1, 2), tau.beta.3 = runif(1, 0.1, 2),
		  tau.phi = runif(1, 0.1, 2),
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
						  'phi.mean', 'z.bbs', 'psi', 
                                                  'chi.2.v.2', 'chi.2.rep.v.2'))
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

file.name <- paste('results/icom-BBS-results-', n.iter, '-iterations-', chain, '-chain-', 
		   date, '.R', sep = '')

save(samples, file = file.name)

