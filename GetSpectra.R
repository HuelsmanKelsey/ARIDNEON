# Overview: ---------------------------------------------------------------

# Get Spectra lets you

# Choose a site/location and request from API to download
# open and extract necessary info, interpret the data
# Visualize reflectances and save a simple (5-band) preliminary map


# Load packages -----------------------------------------------------------

#basic packages
library('readr')
library('tidyr') 
library('dplyr')
library('ggplot2')
library('lubridate')

#spatial data
library('sf')                 
library('terra')    
library('mapview')

#for NEON API:
library('neonUtilities')    
library('jsonlite')   

# reading HDF5 files (AOP data)
library('rhdf5')              

# to extract metadata from shape files
library('xml2')

#for vegetation data
library('vegan')

#for parallelization
library('doParallel')
library('foreach')

# Choose a site/location and request from API to download
#can use all_NEON_plots w LL / EN centroids or shapefiles
all_NEON_plots <- read.csv(file = 'All_NEON_TOS_Plot_Centroids_V11.csv')
#NEON shape files: choose one to assign to LL; one is plots, one is subplots
plots_shp <- terra::vect('/Users/khuelsma/Desktop/NEON Spectral Variability/Relevant NEON Materials/NEON TOS Plots/All_NEON_TOS_Plot_Polygons_V11.shp')
#CPER plots 
plots_shp_CPER_plants <- subset(plots_shp, stringr::str_detect(plots_shp$plotID, 'CPER') &
                                  'div' %in% plots_shp$appMods)
subplots_shp <- terra::vect('/Users/khuelsma/Desktop/NEON Spectral Variability/Relevant NEON Materials/NEON TOS Plots/All_NEON_TOS_Plot_Subplots_V11.shp')
#CPER subplots
subplot_shp_CPER_plants <- subset(subplots_shp, stringr::str_detect(subplots_shp$plotID, 'CPER') &
                                    'div' %in% subplots_shp$appMods)
#subplot shapes we can use to extract values
subplot_shp_100 <- subset(subplot_shp_CPER_plants, stringr::str_detect(subplot_shp_CPER_plants$subplotID, '_100'))
subplot_shp_10 <- subset(subplot_shp_CPER_plants, stringr::str_detect(subplot_shp_CPER_plants$subplotID, '_10_'))
subplot_shp_1 <- subset(subplot_shp_CPER_plants, stringr::str_detect(subplot_shp_CPER_plants$subplotID, '_1_'))

