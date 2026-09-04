# Overview: ---------------------------------------------------------------

# Get Spectra lets you

# Choose a site/location and open or request from API to download
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

# Choose a site/location and open or request from API to download
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

#if you have the data already: open it
# using files list: 1) make maps, 2) extract data:

#load file names:
CPER_2024_files <- read_csv('CPER_filelist_2024.csv')
#extract rgb vs. hsi
rgb_files_2024 <- CPER_2024_files$rgb_path
hsi_files_2024 <- CPER_2024_files$path

#start w rgb tiles: for each list item (file path), create a list of tilemaps to merge
tilemap_rgb_list <- foreach (img = rgb_files_2024) %do% {
  tilemap <- terra::rast(img)
  #terra::plotRGB(tilemap)
}
full_site_2024 <-  terra::merge(terra::sprc(tilemap_rgb_list))
#terra::plotRGB(full_site_2021)

#terra::plot(LL_poly$geometry, add = TRUE)
#terra::crs(LL_poly)

#pulling out subplots
subplots_rgb_list <- foreach (img = tilemap_rgb_list) %do% { 
  #individual img is a row from the tilemap rgb list (or whatever is entered above)
  tilemap <- img
  #establish the extent of the rgb tile so you can clip the list of all plots later to just those in the tile
  tile_extent <- terra::ext(tilemap)
  
  #project subplot shape file to the tile map so we can extract
  subplots_shp_aligned <- terra::project(subplot_shp_1, tilemap)
  #DA_shp_aligned <- terra::project(poly_vect, tilemap)
  
  #filter all the subplots to just the ones in the tile's extent
  plots_in_tile <- subplots_shp_aligned[tile_extent]
  #plots_in_tile <- DA_shp_aligned[tile_extent]
  
  if (nrow(plots_in_tile) > 0) { #as long as there are plots in the tile,
    plot_lookup <- data.frame(
      ID = 1:nrow(plots_in_tile), #we will make a dataframe with plot ID (just a number)
      subplotID = plots_in_tile$subplotID, #the subplot and plot ID (from the dataframe)
      plotID = plots_in_tile$plotID #so that we can reference which plot is which
    )
    
    foreach(plot = nrow(plots_in_tile),
            .combine = rbind) %do% {
              subplot <- plots_in_tile[plot,] #for each plot
              subplot_df <- as.data.frame(subplot)
              subplotID = subplot_df$subplotID
              buffer_pt <- terra::buffer(subplot, width = 2) #buffer so we account for spatial uncertainty?
              cropped_plot <- terra::crop(tilemap, buffer_pt) #crop the map to just the buffered point
              plotRGB(cropped_plot) #plot the cropped plot
              #terra::plot(buffer_pt, add = TRUE)
            }
  }
}

#same for HSI
path <- img

subplots_hsi_1m <- foreach (
  img = hsi_files_2024,
  .combine = rbind) %do% {
    
    which_site = 'CPER'
    which_year = 2024
    path <- img
    
    # Safe read to avoid crashes from corrupted files
    file_is_readable <- tryCatch({
      md1 <- rhdf5::h5readAttributes(path, paste0("/", which_site, "/Reflectance/Reflectance_Data"))
      TRUE
    }, error = function(e) {
      warning(paste("\nSkipping corrupted or inaccessible file:", path))
      rhdf5::h5closeAll()
      FALSE
    })
    
    if (file_is_readable) {  #if the file isn't readable go to the next one; if it is, it will continue.
      md2 <- rhdf5::h5readAttributes(path, paste0('/', which_site, '/Reflectance'))
      wv <- rhdf5::h5read(path, paste0('/', which_site, '/Reflectance/Metadata/Spectral_Data'))
      omit_windows <- data.frame(
        omit_1_0 = md2$Band_Window_1_Nanometers[1], omit_1_f = md2$Band_Window_1_Nanometers[2],
        omit_2_0 = md2$Band_Window_2_Nanometers[1], omit_2_f = md2$Band_Window_2_Nanometers[2]
      )
      #put metadata into wv df
      file_wv_df <- data.frame(year = as.numeric(which_year), wv = round(wv$Wavelength), band = paste0('B', sprintf("%03d", 1:426))) %>%
        mutate(
          omit_band = (wv >= omit_windows$omit_1_0 & wv <= omit_windows$omit_1_f) | (wv >= omit_windows$omit_2_0 & wv <= omit_windows$omit_2_f),
          keep_band = (wv > 450) & (wv < 2150),
          data_ignore = md1$Data_Ignore_Value, 
          SF = md1$Scale_Factor,
          special_band = case_when(
            wv == wv[which.min(abs(wv - 630))] ~ 'red', 
            wv == wv[which.min(abs(wv - 800))] ~ 'NIR',
            wv == wv[which.min(abs(wv - 570))] ~ 'green', 
            wv == wv[which.min(abs(wv - 480))] ~ 'blue',
            wv == wv[which.min(abs(wv - 531))] ~ 'PRI'
          ))
      
      kept_bands <- file_wv_df %>% filter(keep_band == TRUE & omit_band == FALSE)
      raw_data <- rhdf5::h5read(path, paste0("/", which_site, "/Reflectance/Reflectance_Data"))
      reordered_data <- aperm(raw_data, c(3, 2, 1)) #make sure rows and columns are properly transposed in this step
      epsg_code <- rhdf5::h5read(path, paste0("/", which_site, "/Reflectance/Metadata/Coordinate_System/EPSG Code"))
      map_info <- rhdf5::h5read(path, paste0("/", which_site, "/Reflectance/Metadata/Coordinate_System/Map_Info"))
      map_easting <- as.numeric(strsplit(map_info, ',')[[1]][4])
      map_northing <- as.numeric(strsplit(map_info, ',')[[1]][5]) #which is NOT the same as naming convention
      E_min <- map_easting
      E_max <- map_easting + 1000
      N_min <- map_northing - 1000
      N_max <- map_northing
      
      hsi_rast_raw <- terra::rast(reordered_data, crs = paste0('EPSG:', epsg_code))
      rm(raw_data, reordered_data) #once this is turned into a raster, remove it to save memory
      gc() # Force garbage collection
      
      terra::ext(hsi_rast_raw) <- c(E_min, 
                                    E_max,
                                    N_min, 
                                    N_max)
      names(hsi_rast_raw) <- file_wv_df$band
      
      #clean raster
      hsi_rast_keep <- hsi_rast_raw[[unique(kept_bands$band)]] 
      hsi_rast_ignore <- terra::subst(hsi_rast_keep, unique(kept_bands$data_ignore), NA) 
      hsi_rast <- hsi_rast_ignore / unique(kept_bands$SF)
      tile_extent <- terra::ext(hsi_rast)
      
      subplots_shp_aligned <- terra::project(subplot_shp_1, hsi_rast) #can adjust to any subplot input
      
      plots_in_tile <- subplots_shp_aligned[tile_extent]
      
    } #if readable
    if (nrow(plots_in_tile) > 0) {
      plot_lookup <- data.frame(
        ID = 1:nrow(plots_in_tile),
        plotID = plots_in_tile$plotID,
        subplotID = plots_in_tile$subplotID
      )
      #extract the data for a nice lil dataframe
      extracted_plot_data <- terra::extract(hsi_rast, plots_in_tile, cells = TRUE, xy = TRUE)
      extracted_data_labeled <- merge(extracted_plot_data, plot_lookup, by = "ID") %>%
        mutate(cell_easting = x,
               cell_northing = y)
    }
    #extracted_data_labeled
  }
extracted_1m <- subplots_hsi_1m

#extract wavelengths: only need one image from a year
path <- hsi_files_2024[1]
md1 <- rhdf5::h5readAttributes(path, paste0("/", which_site, "/Reflectance/Reflectance_Data"))
md2 <- rhdf5::h5readAttributes(path, paste0('/', which_site, '/Reflectance'))
wv <- rhdf5::h5read(path, paste0('/', which_site, '/Reflectance/Metadata/Spectral_Data'))
omit_windows <- data.frame(
  omit_1_0 = md2$Band_Window_1_Nanometers[1], omit_1_f = md2$Band_Window_1_Nanometers[2],
  omit_2_0 = md2$Band_Window_2_Nanometers[1], omit_2_f = md2$Band_Window_2_Nanometers[2]
)
#put metadata into wv df
file_wv_df <- data.frame(year = as.numeric(which_year), wv = round(wv$Wavelength), band = paste0('B', sprintf("%03d", 1:426))) %>%
  mutate(
    omit_band = (wv >= omit_windows$omit_1_0 & wv <= omit_windows$omit_1_f) | (wv >= omit_windows$omit_2_0 & wv <= omit_windows$omit_2_f),
    keep_band = (wv > 450) & (wv < 2150),
    data_ignore = md1$Data_Ignore_Value, 
    SF = md1$Scale_Factor,
    special_band = case_when(
      wv == wv[which.min(abs(wv - 630))] ~ 'red', 
      wv == wv[which.min(abs(wv - 800))] ~ 'NIR',
      wv == wv[which.min(abs(wv - 570))] ~ 'green', 
      wv == wv[which.min(abs(wv - 480))] ~ 'blue',
      wv == wv[which.min(abs(wv - 531))] ~ 'PRI'
    ))
kept_bands <- file_wv_df %>% filter(keep_band == TRUE & omit_band == FALSE)

extracted_1m_long <- extracted_1m %>%
  pivot_longer(cols = unique(kept_bands$band), names_to = 'band', values_to = 'refl') %>%
  left_join(kept_bands) %>%
  select(cell, cell_easting, cell_northing,
         plotID, subplotID,
         band, year, wv, refl) %>%
  mutate(eventID = paste0(plotID, '_', subplotID, '_', year))

summary_1m <- extracted_1m_long %>%
  group_by(eventID, wv) %>%
  summarise(mean_refl = mean(refl))
#if used buffer and sample more than 1 cell
#            sd_refl = sd(refl),
#          se_refl = sd_refl/sqrt(n()),
#           CV_refl = sd_refl/mean_refl)

#if you need to download it:

