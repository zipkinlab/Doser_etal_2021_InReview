# data-prep.R: this code extracts raw data from HBEF, BBS, and NEON data files
#              into a format suitable for analyses. Data files used in this 
#              script are not included on GitHub as a result of file size limits. 
#              Only summarized files are included on GitHub. 
#              Please contact Jeff Doser for a link to a Dropbox repository if 
#              raw data files are desired.
# Author: Jeffrey W. Doser (doserjef@msu.edu)
# Citation: 

rm(list = ls())
library(raster)
library(auk)
library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)
library(FedData)
library(sp)
library(rgdal)
library(abind)
library(sf)

# Get BBS data in proper format -------------------------------------------
# Read in BBS Data and filter for only New Hampshire
bbs.dat <- read.csv("data/BBS/50-StopData/1997ToPresent_SurveyWide/fifty6.csv")
# New Hampshire is 58
bbs.dat <- bbs.dat %>%
  filter(StateNum == 58)
# Route information
route.dat <- read.csv("data/BBS/routes.csv")
# Only grab routes of interest in the WMNF
route.dat <- route.dat %>%
  filter(RouteName %in% c('FRANCONIA', 'JEFFERS NTCH', 'JEFFERSON HL',
			  'MT CHOCORUA'))
# Weather Data
weather.dat <- read.csv("data/BBS/weather.csv")
# Only grab routes of interest during years of interest
weather.dat <- weather.dat %>%
  filter(StateNum == 58,
	 Route %in% route.dat$Route,
	 Year %in% 2010:2018)
# Only get BBS data from these four routes
bbs.dat <- inner_join(bbs.dat, route.dat, by = c('Route', 'CountryNum', 'StateNum'))
# Only grab BBS data for years of interest
bbs.dat <- bbs.dat %>%
  filter(Year %in% 2010:2018) %>%
  select(Route, Year, AOU, starts_with('Stop'), Latitude, Longitude)

# Fill in missing rows for all combos of species, year, and route
bbs.dat <- bbs.dat %>%
  complete(AOU, nesting(Route, Year, Latitude, Longitude))
# Replace NAs with 0s for all columns at once.
bbs.dat <- bbs.dat %>%
  replace(is.na(.), 0)
# Filter out for species of interest
# American Redstart = 6870, Black and White Warbler = 6360, BHVI = 6290,
# Blackburnian warbler = 6620, Blackpoll warbler = 6610,
# BTBW = 6540, BTNW = 6670, Canada Warbler = 6860, MAWA = 6570,
# NAWA = 6450, OVEN = 6740, REVI = 6240
# Create a table for species information
sp.info <- data.frame(AOU = c(6870, 6360, 6290, 6620, 6610, 6540, 6670, 6860, 6570,
		              6450, 6740, 6240),
		      code = c('AMRE', 'BAWW', 'BHVI', 'BLBW', 'BLPW',
			       'BTBW', 'BTNW', 'CAWA', 'MAWA', 'NAWA',
			       'OVEN', 'REVI'),
		      name = c('American Redstart', 'Black-and-white Warbler',
			       'Blue-headed Vireo', 'Blackburnian Warbler',
			       'Blackpoll Warbler', 'Black-throated Blue Warbler',
			       'Black-throated Green Warbler', 'Canada Warbler',
			       'Magnolia Warbler', 'Nashville Warbler', 'Ovenbird',
			       'Red-eyed Vireo'))
# Join species-specific info to full data set
bbs.dat <- inner_join(bbs.dat, sp.info, by = 'AOU')
# Select columns of interest
bbs.dat <- bbs.dat %>%
  select(AOU, Route, Year, starts_with('Stop'), Latitude, Longitude, code, name)
# Complete data set with missing combos of Route, Year, and AOU
bbs.dat <- bbs.dat %>%
  complete(AOU, Year, nesting(Route, Latitude, Longitude))

# Get BBS Stop coordinates from Shapefile ---------------------------------
# BBS stop coordinates are not provided by USGS, and only the beginning of 
# each route and the route shapefile is provided. Each stop is approximately 
# 800m apart, for a total of 50 stops per route. In this analysis, we assume 
# stops are evenly distributed across the route, and assign each stop to this
# specific geographic location. This induces a small degree of spatial error
# in the locations, which we believe is negligible for our current purposes.
# Read in route shapefile. 
bbs.locs <- st_read('data/BBS/coordinates/bbsrtsl020.shp')
bbs.locs <- bbs.locs %>%
  filter(RTENAME %in% c('FRANCONIA', 'JEFFERS NTCH', 'JEFFERSON HL',
			'MT CHOCORUA'))
# Convert to sf
bbs.sf <- bbs.locs %>% st_geometry()

# Reset the projection for use with st_line_sample
bbs.sf <- st_transform(bbs.sf, 26919)
# Getting 50 evenly spaced coordinates along the BBS Route
bbs.stop.locs <- st_line_sample(bbs.sf, n = 50)
# Get coordinates in matrix format.
stop.coords <- as.data.frame(st_coordinates(bbs.stop.locs))
# Assign route numbers
stop.coords$Route <- 12
stop.coords$Route[stop.coords$L1 == 1] <- 17
stop.coords$Route[stop.coords$L1 == 2] <- 116
stop.coords$Route[stop.coords$L1 == 3] <- 18
stop.coords$stop <- rep(paste('Stop', 1:50, sep = ''), times = 4)

# Get data in long format
bbs.dat.long <- pivot_longer(bbs.dat, cols = starts_with('Stop'), names_to = 'stop',
			     values_to = 'count')
# Join actual BBS data with the derived stop coordinates. 
bbs.dat.long <- inner_join(bbs.dat.long, stop.coords, by = c('Route', 'stop'))
bbs.dat.long <- bbs.dat.long %>%
  unite('site', Route, stop, remove = FALSE)
# Join BBS data to weather data. 
weather.dat <- weather.dat %>% unite('date', sep = '-', Year, Month, Day, remove = FALSE)
# Get Julian date
weather.dat$date <- as.Date(weather.dat$date, tz = "America/New_York")
weather.dat$julian <- as.numeric(format(weather.dat$date, '%j'))
bbs.dat.long <- left_join(bbs.dat.long, weather.dat, by = c('Route', 'Year'))
# Convert data frame to multi-dimensional array
J.bbs <- 200 # number of route x stops.
# Number of species
I <- nrow(sp.info)
# Species names
sp <- sp.info$AOU
# Specific years
years <- unique(bbs.dat.long$Year)
# Number of years
n.years <- length(years)
# Site ids
sites <- unique(bbs.dat.long$site)
# Arrays for data
v.2 <- array(NA, dim = c(I, J.bbs, n.years))
day.bbs <- array(NA, dim = c(J.bbs, n.years))
stop.bbs <- as.numeric(word(sites, 2, sep = 'Stop'))
obsv.bbs <- array(NA, dim = c(J.bbs, n.years))
for (j in 1:J.bbs) {
  print(j)
  for (t in 1:n.years) {
      day.bbs[j, t] <- bbs.dat.long %>%
        filter(Year == years[t], site == sites[j]) %>%
	pull(julian) %>%
	unique()
      obsv.bbs[j, t] <- bbs.dat.long %>%
	filter(Year == years[t], site == sites[j]) %>%
	pull(ObsN) %>%
	unique()
    for (i in 1:I) {
      tmp <- bbs.dat.long %>%
        filter(AOU == sp[i],
	       Year == years[t],
	       site == sites[j])
      v.2[i, j, t] <- tmp$count
    } # k
  } # t
} # i

# Read in the rest of the data ----------------------------------------------
# Hubbard Brook data  -----------------
hb.dat <- readRDS("data/hubbard-brook/Abundance_array_Zipkin_species.rds")
# Years of interest
hb.dat <- hb.dat[, , 12:20, ]
# Remove Myrtle Warbler, not part of case study. 
hb.dat <- hb.dat[, , , -which(attr(hb.dat, 'dimnames')[[4]] == 'MYWA')]

# White Mountain National Forest Data -------------------------------------
# Note: this is NAD83
white.mtn <- readOGR("data/wmnf/", "wmnf_boundary")

# Hubbard Brook Data ------------------------------------------------------
# Covariates for HB data
hb.time <- readRDS('data/hubbard-brook/TimeOfCount_Zipkin_species.rds')
# Subset to correct years
hb.time <- hb.time[, , 12:20]
# A copy for later on. 
hb.time.copy <- hb.time
# Time of day for HB data
hb.day <- readRDS('data/hubbard-brook/DateOfCount_Zipkin_species.rds')
hb.day <- hb.day[, , 12:20]
# Additional covariates (contains lat and lon values, and elevation) 
hb.covs <- readRDS('data/hubbard-brook/plot_attributes_Zipkin.rds')

# Get time of day in proper format
hb.time <- ifelse(str_length(hb.time) == 3, paste('0', hb.time, sep = ''), hb.time)
hb.time <- ifelse(str_length(hb.time) == 4, paste(substr(hb.time, 1, 2), ':', 
						  substr(hb.time, 3, 4), 
						  sep = ''), hb.time)

# Convert time to number of minutes since midnight. 
hb.time.sec <- hb.time %>% hm() %>% period_to_seconds()
hb.time.min <- hb.time.sec / 60
hb.tod <- array(hb.time.min, dim = dim(hb.time.copy))

# NEON Data ---------------------------------------------------------------
neon.dat <- read.delim("data/NEON-Birds2018/filesToStack10003/stackedFiles/brd_countdata.csv", 
		       sep = ",")

#Manipulating time and date data 
#Specifying the format of time and date data
neon.dat$startDate <- as.POSIXct(neon.dat$startDate,
                                 format="%Y-%m-%d T %H Z",
                                 tz = 'GMT' 
)
# Convert time to EST
neon.dat$time.EST <- with_tz(neon.dat$startDate)
# Get other covariate information
neon.dat$year <- year(neon.dat$time.EST)
neon.dat$month <- month(neon.dat$time.EST)
neon.dat$day <- day(neon.dat$time.EST)
neon.dat$hour <- hour(neon.dat$time.EST)
neon.dat$julian.day <- as.numeric(strftime(neon.dat$time.EST, format = '%j'))

# Just a copy
neon.full <- neon.dat
# Only select data from Bartlett Forest
neon.dat <- neon.dat %>% filter(siteID == 'BART', taxonID %in% sp.info$code)
# Grab columns of interest
neon.dat <- neon.dat %>%
  select(plotID, pointID, pointCountMinute, taxonID, clusterSize, time.EST, 
	 year, month, day, hour)

# NEON data are collected following a removal sampling approach. This part of the 
# code extracts this information. 
# Number of sites
J.neon <- nrow(unique(neon.dat[, c('plotID', 'pointID')])) 
# Unique NEON sites
sites.neon <- unique(neon.dat[, c('plotID', 'pointID')]) %>% arrange(plotID, pointID)
# Number of time intervals for removal sampling
K.neon <- 3
intervals.neon <- list(one = c(1, 2), two = c(3, 4), three = c(5, 6))
# Years
T.neon <- length(unique(neon.dat$year))
years.neon <- sort(unique(neon.dat$year))
# Species codes
sp <- sort(sp.info$code)
v.1 <- array(NA, c(J.neon, K.neon, T.neon, I))
hour.neon <- array(NA, c(J.neon, T.neon))
day.neon <- array(NA, c(J.neon, T.neon))

# Fill in data set
for (j in 1:J.neon) {
  print(paste('Currently on site ', j, ' out of ', J.neon, sep = ''))
  for (t in 1:T.neon) {
    # To determine if point was actually sampled in given year
    tmp.1 <- neon.full %>%
      filter(plotID == sites.neon$plotID[j], 
             pointID == sites.neon$pointID[j], 
	     year == years.neon[t])
    hour.neon[j, t] <- ifelse(length(unique(tmp.1$hour)) > 0, unique(tmp.1$hour), NA)
    day.neon[j, t] <- ifelse(length(unique(tmp.1$julian.day)) > 0, unique(tmp.1$julian.day), NA)
    for (k in 1:K.neon) {
      for (i in 1:I) {
        tmp <- neon.dat %>%
	  filter(plotID == sites.neon$plotID[j], 
		 pointID == sites.neon$pointID[j], 
		 pointCountMinute %in% intervals.neon[[k]], 
		 year == years.neon[t], 
		 taxonID == sp[i])
        v.1[j, k, t, i] <- ifelse(nrow(tmp.1) == 0, NA, 
			             ifelse(nrow(tmp) == 0, 0, sum(tmp$clusterSize)))
      } # i
    } # j
  } # t
} # r

# Get lat-longitude for neon data
neon.site.dat <- read.csv("data/NEON-Birds2018/filesToStack10003/stackedFiles/brd_perpoint.csv")
# Contains other potentially useful covariates as well. Only grab data from Bartlett Forest
neon.site.dat <- neon.site.dat %>% filter(siteID == 'BART')

# Download NED and NLCD Data -------------------------------------------------------
# These are at a resolution of 30 x 30 m
ned.dat <- get_ned(template = white.mtn, label = 'WM')
nlcd.dat <- get_nlcd(template = white.mtn, label = 'WM')

# Get elevation and Percent Forest for each point count location
# Hubbard Brook -----------------------
hb.covs.sp <- hb.covs
# Put the correct coordinates in for two locations missing in the original data set
hb.covs.sp[4, c(2, 3)] <- c(-71.74122, 43.94012)
hb.covs.sp[320, c(2, 3)] <- c(-71.70435, 43.95123)
# Convert to spatial data frame
coordinates(hb.covs.sp) <- ~Longitude + Latitude
proj4string(hb.covs.sp) <- CRS("+proj=longlat +datum=WGS84")
hb.covs.sp <- spTransform(hb.covs.sp, "+proj=utm +zone=19 +units=m +datum=NAD83")
# Extract elevation at HBEF locations
hb.elev <- raster::extract(ned.dat, hb.covs.sp)
# Extract percent forest within a 250m radius at HBEF locations
hb.buffer <- st_buffer(st_as_sf(hb.covs.sp), 250)
props <- function(a, na.rm = TRUE) {
  my.sum <- sum(!is.na(a))	
  prop.for <- sum(a %in% c(41, 42, 43), na.rm = na.rm) / my.sum
  return(prop.for)
}

hb.for <- raster::extract(nlcd.dat, hb.buffer, props)


# NEON --------------------------------
# NEON only provides coordinates of the center plot, thus need to do some 
# manipulation to get coordinates for each specific point count location
neon.covs <- neon.site.dat %>%
  group_by(plotID, pointID) %>%
  summarize(latitude = unique(decimalLatitude), 
	    longitude = unique(decimalLongitude)) %>%
  ungroup()
coordinates(neon.covs) <- ~longitude + latitude
proj4string(neon.covs) <- CRS("+proj=longlat +datum=WGS84")
neon.covs <- spTransform(neon.covs, "+proj=utm +zone=19 +units=m +datum=NAD83")
neon.covs <- as.data.frame(neon.covs) %>%
  mutate(easting = longitude, 
	 northing = latitude)
# Get location of each point count based on NEON protocol. 
neon.covs$easting <- ifelse(neon.covs$pointID %in% c('A1', 'B1', 'C1'), neon.covs$easting - 250, 
			    ifelse(neon.covs$pointID %in% c('A3', 'B3', 'C3'), neon.covs$easting + 250, 
							    neon.covs$easting))
neon.covs$northing <- ifelse(neon.covs$pointID %in% c('A1', 'A2', 'A3'), neon.covs$northing + 250, 
			     ifelse(neon.covs$pointID %in% c('C1', 'C2', 'C3'), neon.covs$northing - 250, 
				    neon.covs$northing))
# Convert NEON data to spatial data frame
coordinates(neon.covs) <- ~easting + northing
proj4string(neon.covs) <- CRS("+proj=utm +zone=19 +units=m +datum=NAD83")
# Extract elevation at NEON sites
neon.elev <- raster::extract(ned.dat, neon.covs)
# Extract local forest cover at NEON sites
neon.buffer <- st_buffer(st_as_sf(neon.covs), 250)
neon.for <- raster::extract(nlcd.dat, neon.buffer, props)

# BBS ---------------------------------
# Convert to spatial data frame
coordinates(stop.coords) <- ~X + Y
proj4string(stop.coords) <- proj4string(white.mtn)
# Extract elevation at BBS sites
bbs.elev <- raster::extract(ned.dat, stop.coords)
# Extract local forest cover at BBS sites
bbs.buffer <- st_buffer(st_as_sf(stop.coords), 250)
bbs.for <- raster::extract(nlcd.dat, bbs.buffer, props)

# Save results ------------------------------------------------------------
save(hb.covs, hb.dat, hb.tod, hb.day, hb.elev, hb.for,
     day.neon, hour.neon, neon.elev, neon.for, neon.covs, v.1, 
     v.2, day.bbs, stop.bbs, obsv.bbs, bbs.elev, bbs.for, 
     file = 'data/final-bird-data.R') 

