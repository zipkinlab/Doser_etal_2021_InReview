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
# For temperature data
library(prism)
load("bird-500m-3km-model-data.R")

white.mtn <- readOGR("wmnf/", "wmnf_boundary")

# Get climate variables ---------------------------------------------------
# These data are at a resolution of 4km x 4km
# Set directory for prism data to go
prism_set_dl_dir('.')

# Temperature: monthly average temp ---
get_prism_monthlys(type = 'tmean', year = 2010:2018, mon = 5, keepZip = FALSE)
tmeans <- prism_archive_subset(type = 'tmean', temp_period = 'monthly', mon = 5, 
				 years = 2010:2018)
path.tmeans <- sapply(tmeans, pd_to_file)
tmeans.rast <- lapply(path.tmeans, raster)

# Temperature: 30 year normals --------
get_prism_normals(type = 'tmean', resolution = '4km', mon = 5, 
		  keepZip = FALSE)
tnorms <- prism_archive_subset(type = 'tmean', temp_period = 'monthly normals', mon = 5, 
			       resolution = '4km')
path.tnorms <- sapply(tnorms, pd_to_file)
tnorms.rast <- raster(path.tnorms)

# Precipitation. This is total monthly precipitation
get_prism_monthlys(type = 'ppt', year = 2010:2018, mon = 5, keepZip = FALSE)
ppts <- prism_archive_subset(type = 'ppt', temp_period = 'monthly', mon = 5, 
				 years = 2010:2018)
path.ppts <- sapply(ppts, pd_to_file)
ppts.rast <- lapply(path.ppts, raster)

# Precipitation: 30 year normals ------
get_prism_normals(type = 'ppt', resolution = '4km', mon = 5, 
		  keepZip = FALSE)
pptnorms <- prism_archive_subset(type = 'ppt', temp_period = 'monthly normals', mon = 5, 
			       resolution = '4km')
path.pptnorms <- sapply(pptnorms, pd_to_file)
pptnorms.rast <- raster(path.pptnorms)

# Clip climate variables to WMNF
white.mtn.grs80 <- spTransform(white.mtn, proj4string(tmeans.rast[[1]]))
tmean.rast.wmnf <- lapply(tmeans.rast, crop, extent(white.mtn.grs80))
ppt.rast.wmnf <- lapply(ppts.rast, crop, extent(white.mtn.grs80))
tnorm.rast.wmnf <- crop(tnorms.rast, extent(white.mtn.grs80))
pptnorm.rast.wmnf <- crop(pptnorms.rast, extent(white.mtn.grs80))

# Convert climate variables to NAD 83, UTM 19
tmean.wmnf.nad <- lapply(tmean.rast.wmnf, projectRaster, crs = proj4string(white.mtn))
ppt.wmnf.nad <- lapply(ppt.rast.wmnf, projectRaster, crs = proj4string(white.mtn))
tnorm.wmnf.nad <- projectRaster(tnorm.rast.wmnf, crs = proj4string(white.mtn))
pptnorm.wmnf.nad <- projectRaster(pptnorm.rast.wmnf, crs = proj4string(white.mtn))
# Compute centroid of each pixel for association with temp and rain
pixels$centroid.x <- (pixels$xmin + pixels$xmax) / 2
pixels$centroid.y <- (pixels$ymin + pixels$ymax) / 2

coordinates(pixels) <- ~centroid.x + centroid.y
proj4string(pixels) <- proj4string(white.mtn)

ppt.vals <- sapply(ppt.wmnf.nad, raster::extract, pixels)
colnames(ppt.vals) <- NULL
tmean.vals <- sapply(tmean.wmnf.nad, raster::extract, pixels)
colnames(tmean.vals) <- NULL
ppt.norm.vals <- raster::extract(pptnorm.wmnf.nad, pixels)
t.norm.vals <- raster::extract(tnorm.wmnf.nad, pixels)

# Subtract 30 year normals from ppt and tmean values
ppt.vals <- ppt.vals - ppt.norm.vals
tmean.vals <- tmean.vals - t.norm.vals

pixels <- as.data.frame(pixels)

# Save output 
save(hb.covs, hb.dat, hb.tod, hb.day, y.ebird, cells, ppt.vals, tmean.vals,
     pixels, obsv, dist.trav, day.eb, time.eb, length.eb, ppt.norm.vals, t.norm.vals,
     neon.covs, day.neon, hour.neon, y.neon, file = 'bird-500m-3km-model-data.R')
