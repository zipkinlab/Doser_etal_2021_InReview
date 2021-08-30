# main-HBEF-NEON-BBS.R: code to run integrated species distribution model with
#                       HBEF, NEON, and BBS data for cross validation. 
#                       See code/wmnf/main-HBEF-NEON-BBS.R for more 
#                       detailed commenting regarding different objects. 
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
inputs <- as.numeric(commandArgs(trailingOnly = TRUE))
if(length(inputs) != 2) {
  base::stop('Need to tell NIMBLE the species number and hold out number')
}
curr.sp <- inputs[1]
curr.set <- inputs[2]

# curr.sp <- 12
# curr.set <- 1

# Load in all the data. 
load("data/final-bird-data.R")
source("code/cross-val-idm/nimble-code/idm-HBEF-NEON-BBS.R")

# Data prep for analysis --------------------------------------------------
# See scripts in wmnf for more detailed commenting. 
# Number of species
I <- dim(hb.dat)[4]
# Number of HBEF sites
J.hbef <- dim(hb.dat)[1]
# Repeat visits at HBEF
K.hbef <- dim(hb.dat)[2]
# Number of BBS sites
J.bbs <- dim(v.2)[2]
# Number of NEON sites
J.neon <- dim(v.1)[1]
# Total number of sites
J <- J.hbef + J.neon + J.bbs
# Number of years
n.years <- dim(hb.dat)[3]
# Specific years
years <- as.numeric(attr(hb.day, 'dimnames')[[3]])
# Standardize year covariate
years <- c(scale(years))
# Get elevation for all sites
# Sites are organized by HBEF, NEON, and BBS
elev <- c(hb.elev, neon.elev, bbs.elev)
elev.real <- elev
elev <- c(scale(elev))
# Percent forest
percent.for <- c(hb.for, neon.for, bbs.for)
percent.for.real <- percent.for
percent.for <- c(scale(percent.for))
# Standardize hubbard brook covariates
hbef.day <- array(scale(hb.day), dim = c(J.hbef, K.hbef, n.years))
hbef.tod <- array(scale(hb.tod), dim = c(J.hbef, K.hbef, n.years))
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
hbef.dat <- aperm(hb.dat, c(1, 2, 4, 3))
v.1 <- aperm(v.1, c(1, 4, 2, 3))
v.2 <- aperm(v.2, c(2, 1, 3))

# Convert all count data to hbefection/non-hbefection data. 
v.1 <- ifelse(v.1 > 0, 1, v.1)
y <- ifelse(hbef.dat > 0, 1, hbef.dat)
v.2 <- ifelse(v.2 > 0, 1, v.2)
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

# Create Cross Validation Data Sets ---------------------------------------
set.seed(2021)
chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
hbef.nums <- sample(1:J.hbef, replace = FALSE)
n.x.val <- 2
hbef.x.indices <- chunk2(hbef.nums, n.x.val)
bbs.nums <- sample(1:J.bbs, replace = FALSE)
bbs.x.indices <- chunk2(bbs.nums, n.x.val)
neon.nums <- sample(1:J.neon, replace = FALSE)
neon.x.indices <- chunk2(neon.nums, n.x.val)
# Save cross validation sets to maintain consistency across models
save(hbef.x.indices, bbs.x.indices, neon.x.indices, file = 'data/cross-val-indices.R')

y.fit <- y[-hbef.x.indices[[curr.set]], , , ]
y.pred <- y[hbef.x.indices[[curr.set]], , , ]
v.2.fit <- v.2[-bbs.x.indices[[curr.set]], , ]
v.2.pred <- v.2[bbs.x.indices[[curr.set]], , ]
v.1.fit <- v.1[-neon.x.indices[[curr.set]], , , ]
v.1.pred <- v.1[neon.x.indices[[curr.set]], , , ]

# Get occupancy covariates ------------
# Elevation
elev.fit <- c(hb.elev[-hbef.x.indices[[curr.set]]], 
	      bbs.elev[-bbs.x.indices[[curr.set]]], 
	      neon.elev[-neon.x.indices[[curr.set]]])
elev.fit.copy <- elev.fit
elev.fit <- (elev.fit - mean(elev.fit)) / sd(elev.fit)
elev.pred <- c(hb.elev[hbef.x.indices[[curr.set]]], 
	       bbs.elev[bbs.x.indices[[curr.set]]], 
	       neon.elev[neon.x.indices[[curr.set]]])
elev.pred <- (elev.pred - mean(elev.fit.copy)) / sd(elev.fit.copy)
# Percent forest
for.fit <- c(hb.for[-hbef.x.indices[[curr.set]]], 
	      bbs.for[-bbs.x.indices[[curr.set]]], 
	      neon.for[-neon.x.indices[[curr.set]]])
for.fit.copy <- for.fit
for.fit <- (for.fit - mean(for.fit)) / sd(for.fit)
for.pred <- c(hb.for[hbef.x.indices[[curr.set]]], 
	       bbs.for[bbs.x.indices[[curr.set]]], 
	       neon.for[neon.x.indices[[curr.set]]])
for.pred <- (for.pred - mean(for.fit.copy)) / sd(for.fit.copy)

# Get hbefection covariates ------------
# HBEF
hbef.day.fit <- hbef.day[-hbef.x.indices[[curr.set]], , ]
hbef.tod.fit <- hbef.tod[-hbef.x.indices[[curr.set]], , ]
hbef.day.pred <- hbef.day[hbef.x.indices[[curr.set]], , ]
hbef.tod.pred <- hbef.tod[hbef.x.indices[[curr.set]], , ]
# BBS
day.bbs.fit <- day.bbs[-bbs.x.indices[[curr.set]], ]
obsv.bbs.fit <- obsv.bbs[-bbs.x.indices[[curr.set]], ]
day.bbs.pred <- day.bbs[bbs.x.indices[[curr.set]], ]
obsv.bbs.pred <- obsv.bbs[bbs.x.indices[[curr.set]], ]
# NEON
neon.day.fit <- neon.day[-neon.x.indices[[curr.set]], ]
neon.hour.fit <- neon.hour[-neon.x.indices[[curr.set]], ]
neon.day.pred <- neon.day[neon.x.indices[[curr.set]], ]
neon.hour.pred <- neon.hour[neon.x.indices[[curr.set]], ]

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
v.1.fit.df$Year.all <- neon.years[v.1.fit.df$Year]
n.vals.neon.fit <- nrow(v.1.fit.df)
v.1.fit.df$pixel <- v.1.fit.df$Site + n_distinct(y.fit.df$Site)
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
v.1.pred.df$Year.all <- neon.years[v.1.pred.df$Year]
v.1.pred.df$pixel <- v.1.pred.df$Site + n_distinct(y.pred.df$Site)
v.1.pred.df <- v.1.pred.df %>% arrange(Species, Year)
v.1.pred.df <- v.1.pred.df %>%
  filter(Species %in% v.1.fit.sp.obs)
v.1.pred.df$p.sp <- as.numeric(factor(v.1.pred.df$Species))
n.vals.neon.pred <- nrow(v.1.pred.df)
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
v.2.fit.df$pixel <- v.2.fit.df$Site + n_distinct(y.fit.df$Site) + 
                      n_distinct(v.1.fit.df$Site)
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
v.2.pred.df$pixel <- v.2.pred.df$Site + n_distinct(y.pred.df$Site) + 
	               n_distinct(v.1.pred.df$Site)
v.2.pred.df <- v.2.pred.df %>% arrange(Species, Year)
v.2.pred.df <- v.2.pred.df %>%
  filter(Species %in% v.2.fit.sp.obs)
v.2.pred.df$p.sp <- as.numeric(factor(v.2.pred.df$Species))
n.vals.bbs.pred <- nrow(v.2.pred.df)


# Work with a single species
sp <- c('AMRE', 'BAWW', 'BHVI', 'BLBW', 'BLPW', 'BTBW',
        'BTNW', 'CAWA', 'MAWA', 'NAWA', 'OVEN', 'REVI')
y.fit.df <- y.fit.df %>%
  filter(Species == which(sp == sp[curr.sp])) %>%
  mutate(Year = as.numeric(factor(Year)))
y.pred.df <- y.pred.df %>%
  filter(Species == which(sp == sp[curr.sp])) %>%
  mutate(Year = as.numeric(factor(Year)))
v.1.fit.df <- v.1.fit.df %>%
  filter(Species == which(sp == sp[curr.sp])) %>%
  mutate(Year = as.numeric(factor(Year)))
v.1.pred.df <- v.1.pred.df %>%
  filter(Species == which(sp == sp[curr.sp])) %>%
  mutate(Year = as.numeric(factor(Year)))
v.2.fit.df <- v.2.fit.df %>%
  filter(Species == which(sp == sp[curr.sp])) %>%
  mutate(Year = as.numeric(factor(Year)))
v.2.pred.df <- v.2.pred.df %>%
  filter(Species == which(sp == sp[curr.sp])) %>%
  mutate(Year = as.numeric(factor(Year)))

# Redfine all constants ---------------------------------------------------
J.hbef.fit <- n_distinct(y.fit.df$Site)
J.neon.fit <- n_distinct(v.1.fit.df$Site)
J.bbs.fit <- n_distinct(v.2.fit.df$Site)
J.fit <- sum(n_distinct(y.fit.df$Site), n_distinct(v.2.fit.df$Site), 
	 n_distinct(v.1.fit.df$Site))
J.pred <- sum(n_distinct(y.pred.df$Site), n_distinct(v.2.pred.df$Site), 
	      n_distinct(v.1.pred.df$Site))
elev.hbef.fit <- elev.fit[1:J.hbef.fit]
elev.neon.fit <- elev.fit[1:J.neon.fit + J.hbef.fit]
elev.bbs.fit <- elev.fit[1:J.bbs.fit + J.hbef.fit + J.neon.fit]
for.hbef.fit <- for.fit[1:J.hbef.fit]
for.neon.fit <- for.fit[1:J.neon.fit + J.hbef.fit]
for.bbs.fit <- for.fit[1:J.bbs.fit + J.hbef.fit + J.neon.fit]
years.neon <- 6:9
n.vals.hbef.fit <- nrow(y.fit.df)
n.vals.hbef.pred <- nrow(y.pred.df)
n.vals.bbs.fit <- nrow(v.2.fit.df)
n.vals.bbs.pred <- nrow(v.2.pred.df)
n.vals.neon.fit <- nrow(v.1.fit.df)
n.vals.neon.pred <- nrow(v.1.pred.df)

# Constants -------------------------------------------------------------
idm.consts <- list(n.years = n.years, n.years.neon = n.years.neon, years.neon = years.neon,
		   J.hbef = J.hbef.fit, n.vals.hbef = n.vals.hbef.fit, 
		   n.vals.neon = n.vals.neon.fit, n.vals.bbs = n.vals.bbs.fit,
		   year.indx.hbef = y.fit.df$Year, 
		   site.hbef = y.fit.df$Site, years = years, 
		   year.indx.bbs = v.2.fit.df$Year, obsv.bbs = v.2.fit.df$obsv, 
		   site.bbs = v.2.fit.df$Site,  
		   year.indx.neon = v.1.fit.df$Year, site.neon = v.1.fit.df$Site, 
		   n.obsv.bbs = n.obsv.bbs, 
		   J.neon = J.neon.fit, J.bbs = J.bbs.fit)
# Data ------------------------------------------------------------------
idm.data <- list(ELEV.hbef = elev.hbef.fit, FOREST.hbef = for.hbef.fit,
		 ELEV.neon = elev.neon.fit, FOREST.neon = for.neon.fit, 
		 ELEV.bbs = elev.bbs.fit, FOREST.bbs = for.bbs.fit, 
		 y = y.fit.df$Count, v.1 = v.1.fit.df$Count,
		 DAY.neon = v.1.fit.df$Day, HOUR.neon = v.1.fit.df$TOD,
		 DAY.hbef = y.fit.df$Day, TOD.hbef = y.fit.df$TOD, 
      	         DAY.bbs = v.2.fit.df$Day, v.2 = v.2.fit.df$Count)
# Initial values --------------------------------------------------------
z.hbef.init <- array(1, dim = c(J.hbef.fit, n.years))
z.neon.init <- array(1, dim = c(J.neon.fit, n.years.neon))
z.bbs.init <- array(1, dim = c(J.bbs.fit, n.years))
idm.inits <- list(z.hbef = z.hbef.init, 
		  z.neon = z.neon.init, 
		  z.bbs = z.bbs.init, 
		  int.beta = runif(n.years, 0.3, 0.9), 
		  beta.1 = rnorm(1), beta.2 = rnorm(1), 
		  beta.3 = rnorm(1), phi = rnorm(1),
		  int.alpha.mean = runif(1, 0.3, 0.9), alpha.1 = rnorm(1), 
		  alpha.2 = rnorm(1), 
		  alpha.3 = rnorm(1), int.gamma.1.mean = runif(1, 0.3, 0.8), 
      	          int.gamma.2.mean = runif(1, 0.2, 0.8), gamma.1.1 = rnorm(1), 
      	          gamma.1.2 = rnorm(1), tau.gamma.2.3 = runif(1, 0.1, 2), 
		  gamma.2.1 = rnorm(1), gamma.2.2 = rnorm(1), 
		  gamma.1.3 = rnorm(1), tau.alpha.0 = runif(1, 0.1, 1), 
		  tau.gamma.1.0 = runif(1, 0.1, 1), tau.gamma.2.0 = runif(1, 0.1, 1))
# Create the model ------------------------------------------------------
idm.model <- nimbleModel(code = idm.code, name = 'idm', constants = idm.consts,
		          data = idm.data, inits = idm.inits)

# Configure MCMC --------------------------------------------------------
idm.conf <- configureMCMC(idm.model, monitors = c('beta.0', 'beta.1', 
      					          'beta.2', 'beta.3', 'phi', 				
						  'int.alpha.mean', 
						  'alpha.1', 'alpha.2', 
						  'alpha.3', 'int.gamma.1.mean', 
						  'gamma.1.1', 'gamma.1.2', 
						  'gamma.1.3', 'int.gamma.2.mean', 
						  'gamma.2.1', 'gamma.2.2', 
						  'gamma.2.3', 'z.hbef', 'z.bbs', 
						  'z.neon'))

# Create an MCMC function -------------------------------------------
idm.mcmc <- buildMCMC(idm.conf)
# Compile model
idm.c.model <- compileNimble(idm.model)
idm.c.mcmc <- compileNimble(idm.mcmc, project = idm.model)
# Number of iterations --------------------------------------------------
n.iter <- 250000
n.burn <- 100000
n.thin <- 20
n.chain <- 1
#set.seed(seed)
samples <- runMCMC(idm.c.mcmc, niter = n.iter, nburnin = n.burn,
	           thin = n.thin, nchains = n.chain, samplesAsCodaMCMC = TRUE)

date <- Sys.Date()
file.name <- paste('results/idm-results-', n.iter, '-iterations-', sp[curr.sp], '-species-', 
		   curr.set, '-holdout-', date, '.R', sep = '')

save(samples, file = file.name)
