# main-BBS.R: code to run community model using only BBS data for cross
#             validation. See code/wmnf/main-BBS.R for more detailed
#             commenting regarding different objects. 
# Author: Jeffrey W. Doser (doserjef@msu.edu)
# Citation: 

rm(list = ls())

library(coda)

library(tidyverse)
library(nimble)
library(coda)
# This index is used to run these multiple hold out datasets across multiple cores. By 
# providing the number of the chain, you will create different output 
# files that you can then combine into a single mcmc.list
curr.set <- as.numeric(commandArgs(trailingOnly = TRUE))
if(length(curr.set) == 0) base::stop('Need to tell NIMBLE the hold out number')

# curr.set <- 1

# Load in all the data. 
load("data/final-bird-data.R")
source("code/cross-val-icom/nimble-code/icom-BBS.R")

# Data prep for analysis --------------------------------------------------
# Number of species
I <- dim(hb.dat)[4]
# Number of BBS sites
J.bbs <- dim(v.2)[2]
# Total number of sites
J <- J.bbs
# Number of years
n.years <- dim(hb.dat)[3]
# Specific years
years <- as.numeric(attr(hb.day, 'dimnames')[[3]])
# Standardize year covariate
years <- c(scale(years))
# Get elevation for all sites
# Sites are organized by DET, NEON, and BBS
elev <- c(bbs.elev)
elev.real <- elev
elev <- c(scale(elev))
percent.for <- c(bbs.for)
percent.for.real <- percent.for
percent.for <- c(scale(percent.for))
# Standardize BBS Data
day.bbs <- array(scale(day.bbs), dim = c(J.bbs, n.years))
n.obsv.bbs <- length(unique(c(obsv.bbs))) 
obsv.bbs <- array(as.numeric(factor(obsv.bbs)), dim = c(J.bbs, n.years))

# Convert all count data to hbection/non-hbection data. 
v.2 <- ifelse(v.2 > 0, 1, v.2)
# Reorder v.2
v.2 <- aperm(v.2, c(2, 1, 3))


# Create Cross Validation Data Sets ---------------------------------------
load("data/cross-val-indices.R")

v.2.fit <- v.2[-bbs.x.indices[[curr.set]], , ]
v.2.pred <- v.2[bbs.x.indices[[curr.set]], , ]

# Get occupancy covariates ------------
# Elevation
elev.fit <- c(bbs.elev[-bbs.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(bbs.elev[bbs.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(bbs.for[-bbs.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(bbs.for[bbs.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)

# Get deection covariates ------------
# BBS
day.bbs.fit <- day.bbs[-bbs.x.indices[[curr.set]], ]
obsv.bbs.fit <- obsv.bbs[-bbs.x.indices[[curr.set]], ]
day.bbs.pred <- day.bbs[bbs.x.indices[[curr.set]], ]
obsv.bbs.pred <- obsv.bbs[bbs.x.indices[[curr.set]], ]

# Get All Data in Long Format ---------------------------------------------
# Modeling the data in this long format drastically speeds up run times
# BBS ---------------------------------
# Fit 
v.2.fit.df <- as.data.frame.table(v.2.fit)
names(v.2.fit.df) <- c('Site', 'Species', 'Year', 'Count')
v.2.fit.df <- v.2.fit.df %>%
  mutate(Site = as.numeric(Site),
	 Species = as.numeric(Species),
	 Year = as.numeric(Year))

# Add DAY to data frame
day.bbs.fit.df <- as.data.frame.table(day.bbs.fit)
names(day.bbs.fit.df) <- c('Site', 'Year', 'Day')
day.bbs.fit.df <- day.bbs.fit.df %>%
  mutate(Site = as.numeric(Site),
	 Year = as.numeric(Year))
# Join to full data
v.2.fit.df <- left_join(v.2.fit.df, day.bbs.fit.df, by = c('Site', 'Year'))
# Add observer identity to data frame
obsv.bbs.fit.df <- as.data.frame.table(obsv.bbs.fit)
names(obsv.bbs.fit.df) <- c('Site', 'Year', 'obsv')
obsv.bbs.fit.df <- obsv.bbs.fit.df %>%
  mutate(Site = as.numeric(Site),
	 Year = as.numeric(Year))
# Join to full data
v.2.fit.df <- left_join(v.2.fit.df, obsv.bbs.fit.df, by = c('Site', 'Year'))
# Remove rows without data
v.2.fit.df <- v.2.fit.df %>%
  filter(!is.na(Count))
n.vals.bbs.fit <- nrow(v.2.fit.df)
v.2.fit.df$pixel <- v.2.fit.df$Site 
v.2.fit.df <- v.2.fit.df %>% arrange(Species, Year)
# Determine observed species 
v.2.fit.sp.obs <- v.2.fit.df %>%
  group_by(Species) %>%
  summarize(counts = sum(Count)) %>%
  filter(counts > 0) %>%
  pull(Species)
v.2.fit.df <- v.2.fit.df %>%
  filter(Species %in% v.2.fit.sp.obs)
n.sp.bbs <- length(v.2.fit.sp.obs)
v.2.fit.df$p.sp <- as.numeric(factor(v.2.fit.df$Species))
n.vals.bbs.fit <- nrow(v.2.fit.df)
# Pred
v.2.pred.df <- as.data.frame.table(v.2.pred)
names(v.2.pred.df) <- c('Site', 'Species', 'Year', 'Count')
v.2.pred.df <- v.2.pred.df %>%
  mutate(Site = as.numeric(Site),
	 Species = as.numeric(Species),
	 Year = as.numeric(Year))

# Add DAY to data frame
day.bbs.pred.df <- as.data.frame.table(day.bbs.pred)
names(day.bbs.pred.df) <- c('Site', 'Year', 'Day')
day.bbs.pred.df <- day.bbs.pred.df %>%
  mutate(Site = as.numeric(Site),
	 Year = as.numeric(Year))
# Join to full data
v.2.pred.df <- left_join(v.2.pred.df, day.bbs.pred.df, by = c('Site', 'Year'))
# Add observer identity to data frame
obsv.bbs.pred.df <- as.data.frame.table(obsv.bbs.pred)
names(obsv.bbs.pred.df) <- c('Site', 'Year', 'obsv')
obsv.bbs.pred.df <- obsv.bbs.pred.df %>%
  mutate(Site = as.numeric(Site),
	 Year = as.numeric(Year))
# Join to full data
v.2.pred.df <- left_join(v.2.pred.df, obsv.bbs.pred.df, by = c('Site', 'Year'))
# Remove rows without data
v.2.pred.df <- v.2.pred.df %>%
  filter(!is.na(Count), !is.na(Day), !is.na(obsv))
n.vals.bbs.pred <- nrow(v.2.pred.df)
v.2.pred.df$pixel <- v.2.pred.df$Site 
v.2.pred.df <- v.2.pred.df %>% arrange(Species, Year)
v.2.pred.df <- v.2.pred.df %>%
  filter(Species %in% v.2.fit.sp.obs)
v.2.pred.df$p.sp <- as.numeric(factor(v.2.pred.df$Species))
n.vals.bbs.pred <- nrow(v.2.pred.df)

# Redfine all constants ---------------------------------------------------
J.fit <- sum(n_distinct(v.2.fit.df$Site))
J.pred <- sum(n_distinct(v.2.pred.df$Site))

# Constants -------------------------------------------------------------
icom.consts <- list(n.years = n.years, I = I, 
		   J = J.fit, n.vals.bbs = n.vals.bbs.fit,
		   years = years, 
		   sp.indx.bbs = v.2.fit.df$Species, 
		   year.indx.bbs = v.2.fit.df$Year, obsv.bbs = v.2.fit.df$obsv, 
		   site.bbs = v.2.fit.df$Site, n.obsv.bbs = n.obsv.bbs, 
                   I.bbs = n.sp.bbs, 
                   sp.indx.bbs.p = v.2.fit.df$p.sp)

# Data ------------------------------------------------------------------
icom.data <- list(ELEV = elev.fit, FOREST = for.fit,
      	         DAY.bbs = v.2.fit.df$Day, v.2 = v.2.fit.df$Count)
# Initial values --------------------------------------------------------
z.init <- array(1, dim = c(I, J.fit, n.years))
icom.inits <- list(z.bbs = z.init, int.beta.mean = runif(n.years, 0.3, 0.9), 
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
						  'phi.mean', 'z.bbs'))
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

file.name <- paste('results/icom-bbs-results-', n.iter, '-iterations-', curr.set, '-holdout-', 
		   date, '.R', sep = '')

save(samples, file = file.name)

