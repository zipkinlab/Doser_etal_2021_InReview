# main-NEON.R: code to run community model using only NEON data for cross 
#              validation. See code/wmnf/main-NEON.R for more detailed
#              commenting regarding different objects. 
# Author: Jeffrey W. Doser (doserjef@msu.edu)
# Citation: 

rm(list = ls())

library(coda)

library(tidyverse)
library(nimble)
library(coda)
# This index is used to run multiple hold out datasets across multiple cores. By 
# providing the number of the chain, you will create different output 
# files that you can then combine into a single mcmc.list
curr.set <- as.numeric(commandArgs(trailingOnly = TRUE))
if(length(curr.set) == 0) base::stop('Need to tell NIMBLE the hold out number')

# curr.set <- 1

# Load in all the data. 
load("data/final-bird-data.R")
source("code/cross-val-icom/nimble-code/icom-NEON.R")

# Data prep for analysis --------------------------------------------------
# Number of species
I <- dim(hb.dat)[4]
# Number of NEON sites
J.neon <- dim(v.1)[1]
# Total number of sites
J <- J.neon
# Number of years
n.years <- dim(v.1)[3]
# Specific years
years <- 2015:2018
n.years <- length(years)
# Standardize year covariate
years <- c(scale(years))
# Get elevation for all sites
# Sites are organized by HBEF, NEON, and BBS
elev <- c(neon.elev)
elev.real <- elev
elev <- c(scale(elev))
percent.for <- c(neon.for)
percent.for.real <- percent.for
percent.for <- c(scale(percent.for))
# The indices of the years for which we have NEON data
neon.years <- 6:9
# Standardize covariates abundance covariates
years <- c(scale(years))
# Standardize neon data
n.years.neon <- dim(v.1)[3]
years.neon <- years[6:9]
neon.day <- array(scale(day.neon), dim = c(J.neon, n.years.neon))
neon.hour <- array(scale(hour.neon), dim = c(J.neon, n.years.neon))
# Reorganize arrays
v.1 <- aperm(v.1, c(1, 4, 2, 3))

# Convert all count data to hbefection/non-hbefection data. 
v.1 <- ifelse(v.1 > 0, 1, v.1)
# Set all subsequent hbefections of a species to 0 if hbefected in 
# previous time bin. 
for (j in 1:J.neon) {
  for (i in 1:I) {
    for (t in 1:n.years) {
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

# Create Cross Validation Data Sets ---------------------------------------
load("data/cross-val-indices.R")
v.1.fit <- v.1[-neon.x.indices[[curr.set]], , , ]
v.1.pred <- v.1[neon.x.indices[[curr.set]], , , ]

# Get occupancy covariates ------------
# Elevation
elev.fit <- c(neon.elev[-neon.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(neon.elev[neon.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(neon.for[-neon.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(neon.for[neon.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)

# Get hbefection covariates ------------
# NEON
neon.day.fit <- neon.day[-neon.x.indices[[curr.set]], ]
neon.hour.fit <- neon.hour[-neon.x.indices[[curr.set]], ]
neon.day.pred <- neon.day[neon.x.indices[[curr.set]], ]
neon.hour.pred <- neon.hour[neon.x.indices[[curr.set]], ]

# Get All Data in Long Format ---------------------------------------------
# Modeling the data in this long format drastically speeds up run times
# NEON --------------------------------
# Fit
v.1.fit.df <- as.data.frame.table(v.1.fit)
names(v.1.fit.df) <- c('Site', 'Species', 'Visit', 'Year', 'Count')
v.1.fit.df <- v.1.fit.df %>%
  mutate(Site = as.numeric(Site),
	 Species = as.numeric(Species),
	 Visit = as.numeric(Visit),
	 Year = as.numeric(Year))

# Add DAY to data frame
neon.day.fit.df <- as.data.frame.table(neon.day.fit)
names(neon.day.fit.df) <- c('Site', 'Year', 'Day')
neon.day.fit.df <- neon.day.fit.df %>%
  mutate(Site = as.numeric(Site),
	 Year = as.numeric(Year))
# Join to full data
v.1.fit.df <- left_join(v.1.fit.df, neon.day.fit.df, by = c('Site', 'Year'))
# Add TOD to data frame
neon.hour.fit.df <- as.data.frame.table(neon.hour.fit)
names(neon.hour.fit.df) <- c('Site', 'Year', 'TOD')
neon.hour.fit.df <- neon.hour.fit.df %>%
  mutate(Site = as.numeric(Site),
	 Year = as.numeric(Year))
# Join to full data
v.1.fit.df <- left_join(v.1.fit.df, neon.hour.fit.df, by = c('Site', 'Year'))
# Remove rows without data
v.1.fit.df <- v.1.fit.df %>%
  filter(!is.na(Count))
n.vals.neon <- nrow(v.1.fit.df)
n.vals.neon.fit <- nrow(v.1.fit.df)
v.1.fit.df$pixel <- v.1.fit.df$Site
v.1.fit.df <- v.1.fit.df %>% arrange(Species, Year)
# Determine observed species 
v.1.fit.sp.obs <- v.1.fit.df %>%
  group_by(Species) %>%
  summarize(counts = sum(Count)) %>%
  filter(counts > 0) %>%
  pull(Species)
v.1.fit.df <- v.1.fit.df %>%
  filter(Species %in% v.1.fit.sp.obs)
n.sp.neon <- length(v.1.fit.sp.obs)
I <- n.sp.neon
v.1.fit.df$p.sp <- as.numeric(factor(v.1.fit.df$Species))
n.vals.neon.fit <- nrow(v.1.fit.df)
# Pred 
v.1.pred.df <- as.data.frame.table(v.1.pred)
names(v.1.pred.df) <- c('Site', 'Species', 'Visit', 'Year', 'Count')
v.1.pred.df <- v.1.pred.df %>%
  mutate(Site = as.numeric(Site),
	 Species = as.numeric(Species),
         Visit = as.numeric(Visit), 
	 Year = as.numeric(Year))

# Add DAY to data frame
neon.day.pred.df <- as.data.frame.table(neon.day.pred)
names(neon.day.pred.df) <- c('Site', 'Year', 'Day')
neon.day.pred.df <- neon.day.pred.df %>%
  mutate(Site = as.numeric(Site),
	 Year = as.numeric(Year))
# Join to full data
v.1.pred.df <- left_join(v.1.pred.df, neon.day.pred.df, by = c('Site', 'Year'))
# Add TOD to data frame
neon.hour.pred.df <- as.data.frame.table(neon.hour.pred)
names(neon.hour.pred.df) <- c('Site', 'Year', 'TOD')
neon.hour.pred.df <- neon.hour.pred.df %>%
  mutate(Site = as.numeric(Site),
	 Year = as.numeric(Year))
# Join to full data
v.1.pred.df <- left_join(v.1.pred.df, neon.hour.pred.df, by = c('Site', 'Year'))
# Remove rows without data
v.1.pred.df <- v.1.pred.df %>%
  filter(!is.na(Count), !is.na(Day), !is.na(TOD))
n.vals.neon <- nrow(v.1.pred.df)
v.1.pred.df$pixel <- v.1.pred.df$Site
v.1.pred.df <- v.1.pred.df %>% arrange(Species, Year)
v.1.pred.df <- v.1.pred.df %>%
  filter(Species %in% v.1.fit.sp.obs)
v.1.pred.df$p.sp <- as.numeric(factor(v.1.pred.df$Species))
n.vals.neon.pred <- nrow(v.1.pred.df)

# Redfine all constants ---------------------------------------------------
J.fit <- sum(n_distinct(v.1.fit.df$Site))
J.pred <- sum(n_distinct(v.1.pred.df$Site))

# Constants -------------------------------------------------------------
icom.consts <- list(n.years = n.years, I = I, 
		   J = J.fit, 
		   n.vals.neon = n.vals.neon.fit, years = years, 
		   sp.indx.neon = v.1.fit.df$p.sp, 
		   year.indx.neon = v.1.fit.df$Year, site.neon = v.1.fit.df$Site, 
		   year.indx.neon.all = v.1.fit.df$Year, I.neon = n.sp.neon, 
                   sp.indx.neon.p = v.1.fit.df$p.sp)

# Data ------------------------------------------------------------------
icom.data <- list(ELEV = elev.fit, FOREST = for.fit,
		 v.1 = v.1.fit.df$Count,
		 DAY.neon = v.1.fit.df$Day, HOUR.neon = v.1.fit.df$TOD)
# Initial values --------------------------------------------------------
z.init <- array(1, dim = c(I, J.fit, n.years))
icom.inits <- list(z.neon = z.init, int.beta.mean = runif(n.years, 0.3, 0.9), 
		  beta.1.mean = rnorm(1), beta.2.mean = rnorm(1), 
		  beta.3.mean = rnorm(1), 
		  phi.mean = rnorm(1),
		  int.gamma.1.mean = runif(1, 0.3, 0.8), 
      	          gamma.1.1.mean = rnorm(1), 
      	          gamma.1.2.mean = rnorm(1), 
		  gamma.1.3.mean = rnorm(1), 
		  tau.beta.0 = runif(n.years, 0.1, 2), tau.beta.1 = runif(1, 0.1, 2),
		  tau.beta.2 = runif(1, 0.1, 2), tau.beta.3 = runif(1, 0.1, 2),
		  tau.phi = runif(1, 0.1, 2),
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
						  'phi.mean',
						  'z.neon'))
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


file.name <- paste('results/icom-NEON-results-', n.iter, '-iterations-', curr.set, '-holdout-', 
		   date, '.R', sep = '')

save(samples, file = file.name)

