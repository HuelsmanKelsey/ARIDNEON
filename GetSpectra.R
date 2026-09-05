# Overview: ---------------------------------------------------------------

# Get Spectra lets you

# Choose a site/location and open or request from API to download; 
#     Ends w check files exist and outputs a list of files
# open files from list and extract necessary info, interpret the data
#     Save summaries and map


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


# Choose a site/location and open or request from API to download  -------
#     Ends w check files exist and outputs a list of files

which_site <- 'CPER' #4 letter NEON code
which_year <- 2024 #there will be more... see next section

# First see what's available ----------------------------------------------
get_prod_avail <- function(which_site, which_dp, BRDF = FALSE) {
  
  ID = case_when(
    which_dp == 'veg' ~ 'DP1.10058.001', #one and only veg
    which_dp == 'refl' & BRDF == TRUE ~ 'DP3.30006.002', #BRDF corrected; unusual case?
    which_dp == 'refl' ~ 'DP3.30006.001') #NOT BRDF corrected
  
  Info <- neonUtilities::getProductInfo(dpID = ID) 
  
  site_avail <- Info$siteCodes %>% 
    select(siteCode,availableMonths) %>%
    group_by(siteCode) %>%
    filter(siteCode == which_site) %>%
    reframe(m = unlist(availableMonths)) %>%
    mutate(date_full = ymd(paste0(m, "-01")), #adds day so lubridate associates w date
           year = year(date_full),
           month = month(date_full)
    ) %>%
    dplyr::select(siteCode, year, month) %>%
    distinct()
}

#this will tell you which years have both veg + AOP data and whether the AOP data are brdf corrected
#NOTE: if the NEON API is down it won't work.
both_avail <- function(which_site) { #default is BRDF = FALSE
  veg <- get_prod_avail(which_site, 'veg')
  nonbrdf <- get_prod_avail(which_site, 'refl', BRDF = FALSE) %>%
    mutate(brdf = 'no')
  brdf <- get_prod_avail(which_site, 'refl', BRDF = TRUE) %>%
    mutate(brdf = 'yes')
  
  both <- rbind(nonbrdf, brdf) %>%
    rename(specmonth = month) %>%
    left_join(veg, by = c('siteCode', 'year'),
              relationship = 'many-to-many') %>%
    as.data.frame() %>%
    group_by(siteCode, year, specmonth, brdf) %>%
    summarise(vegmos = list(unique(month)))
  
  return(both)
}

# Downloading data --------------------------------------------------------
# function to get all easting and northing values for all tiles at *a site*
get_site_EN <- function(which_site, 
                        tiles_only = FALSE) {   #this is for file checking
  
  #filter NEON's plot centroids from the repo to site, get unique values
  plot_polygons <- all_NEON_plots %>% #created above
    filter(grepl('div', appMods)) %>% #filter to just diversity plots for now
    filter(siteID == which_site) %>% #filter to the site
    mutate(
      #fixes Blandy's 17/18 N border (I think)
      utmZone_final = ifelse(which_site == 'BLAN', '17N', utmZone),
      epsg_target = 32600 + as.numeric(gsub("[^0-9]", "", utmZone_final)),
      E_tile = 1000*floor(easting/1000),
      N_tile = 1000*floor(northing/1000),
      crs = unique(epsg_target)
    ) %>%
    select(-utmZone)
  
  if (which_site == 'LL') {
    #LL should be loaded above; need to troubleshoot this.
    
    plot_polygons <- plots_shp
  }
  #if tiles only (for file checking and opening)
  if (tiles_only == TRUE) {
    plot_polygons <- plot_polygons %>%
      distinct(E_tile, N_tile, domainID, crs)
    return(plot_polygons)
  }
  
  return(plot_polygons)
}
#just for fun
check_dims <- function(which_site) {
  obs <- get_site_EN(which_site) # will give you ALL THE DETAILS of obs plots but also the tiles
  tiles_fr_obs <- obs %>%
    distinct(E_tile, N_tile)
  n_tiles_fr_obs <- dim(tiles_fr_obs)[1]
  
  tiles <- get_site_EN(which_site, tiles_only = TRUE) #this will not give you all the details; just tiles
  n_tiles <- dim(tiles)[1]
  if (n_tiles_fr_obs != n_tiles) {
    print('check unique tiles from plots; check n tiles available to download')
  }
  if (n_tiles_fr_obs == n_tiles) {
    print(paste('all', n_tiles_fr_obs, 'tiles have TOS plots :)'))
  }
}

check_files_exist <- function(which_site, 
                              which_year = NULL,
                              RGB = FALSE, 
                              download = FALSE) {  
  
  site_tiles <- get_site_EN(which_site, tiles_only = TRUE)
  #site_year_combos <- both_avail(which_site) 
  
  if(nrow(site_tiles) == 0) {
    warning("No plots found for site: ", which_site)
    return(data.frame(site = character(), year = numeric(), path = character()))
  }
  
  # setwd(home) # Make sure 'home' exists in your global environment!
  
  if (is.null(which_year)) {
    years <- unique(site_year_combos$year)
  } else {
    years <- which_year
  }
  
  all_files_list <- list()
  
  for(i in 1:nrow(site_tiles)) {
    for (y in seq_along(years)) {
      
      which_year <- years[y]
      year <- which_year
      cat(sprintf("\nProcessing year: %s\n", which_year))
      
      #brdf_df <- site_year_combos %>% filter(year == which_year) %>% distinct(brdf)
      #yr_brdf <- unique(brdf_df$brdf)
      #set manually because the API is down
      yr_brdf = 'yes'
      
      h5_final_path <- NA
      rgb_final_path <- NA
      file_found <- FALSE
      
      this_tile <- site_tiles[i,]
      domain = this_tile$domainID
      easting = this_tile$E_tile
      northing = this_tile$N_tile
      DP = ifelse(yr_brdf == 'yes', 'DP3.30006.002', 'DP3.30006.001')
      
      
      for(tile in 1:10) { 
        # 1. Check H5
        h5_path <- ifelse(
          yr_brdf == 'yes', 
          paste0(home, DP, '/neon-aop-provisional-products/', year, '/FullSite/', domain, '/', year, '_', which_site, '_', tile, '/L3/Spectrometer/Reflectance/'),
          paste0(home, DP, '/neon-aop-products/', year, '/FullSite/', domain, '/', year, '_', which_site, '_', tile, '/L3/Spectrometer/Reflectance/'))
        
        h5_base <- paste0('NEON_', domain, '_', which_site, '_DP3_', easting, '_', northing)
        BRDF_filename_append <- '_bidirectional_reflectance'
        h5filename <- ifelse(yr_brdf == 'yes', paste0(h5_base, BRDF_filename_append, '.h5'), paste0(h5_base, '_reflectance.h5'))
        h5_filepath <- paste0(h5_path, h5filename)
        
        if(file.exists(h5_filepath)) {
          if (isTRUE(rhdf5::H5Fis_hdf5(h5_filepath))) {
            h5_final_path <- h5_filepath
            file_found <- TRUE
          } else {
            warning(paste("\nCorrupted HDF5 file detected (will attempt to re-download):", h5_filepath))
          }
        } 
        # 2. Check RGB
        if (RGB == TRUE) {
          rgb_DP <- 'DP3.30010.001'
          rgb_path <- paste0(home, rgb_DP, '/neon-aop-products/', year, '/FullSite/', domain, '/', year, '_', which_site, '_', tile, '/L3/Camera/Mosaic/')
          rgbfilename <- paste0(year, '_', which_site, '_', tile, '_', easting, '_', northing, '_image.tif')
          rgb_filepath <- paste0(rgb_path, rgbfilename)
          
          if(file.exists(rgb_filepath)) {
            rgb_final_path <- rgb_filepath
            # Don't overwrite file_found here to ensure HDF5 triggers download if missing
          } else {
            file_found <- FALSE 
          }
        }
        
        if (file_found == TRUE) break
      } # 1-10 folder tile loop
      
      if (file_found == TRUE) {
        tile_df <- data.frame(
          site = which_site, year = which_year, brdf = yr_brdf,
          E_tile = easting, N_tile = northing, path = h5_final_path,
          stringsAsFactors = FALSE
        )
        if (RGB == TRUE) tile_df$rgb_path <- rgb_final_path
        all_files_list[[paste0(which_year, "_", easting, "_", northing)]] <- tile_df
      } 
      
      # 3. Download if missing
      if(!file_found && download) {
        cat(sprintf("  File not found. Downloading tile: E=%d, N=%d, Year=%s\n", easting, northing, which_year))
        
        tryCatch({
          neonUtilities::byTileAOP(dpID = DP, site = which_site, easting = easting, northing = northing, buffer = 0, check.size = FALSE, include.provisional = TRUE, year = which_year, token = neon_token, progress = TRUE)
        }, error = function(e) warning("HDF5 Download failed"))
        
        if (RGB == TRUE) {
          tryCatch({
            # FIX: Explicitly use RGB DP ID here!
            neonUtilities::byTileAOP(dpID = 'DP3.30010.001', site = which_site, easting = easting, northing = northing, buffer = 0, check.size = FALSE, include.provisional = TRUE, year = which_year, token = neon_token, progress = TRUE)
          }, error = function(e) warning("RGB Download failed"))
        }
      } 
    } 
  } 
  return(do.call(rbind, all_files_list))
}

#hyperspectral: this creates a list of hyperspectral raster or special band rasters of each tile
CPER_spbands_tile_list <- c()
CPER_hsi_tile_list <- c()

#I think this code is the one that doesn't work well?
extract_tile_plots <- function(input) {
  if(nrow(input) == 0) return(data.frame())
  all_plot_data <- list()
  
  for (f in 1:nrow(input)) {
    thisfile <- input[f,] #this file
    
    # Safe read to avoid crashes from corrupted files
    file_is_readable <- tryCatch({
      md1 <- rhdf5::h5readAttributes(thisfile$path, paste0("/", which_site, "/Reflectance/Reflectance_Data"))
      TRUE
    }, error = function(e) {
      warning(paste("\nSkipping corrupted or inaccessible file:", thisfile$path))
      rhdf5::h5closeAll()
      FALSE
    })
    
    if (!file_is_readable) next
    
    md2 <- rhdf5::h5readAttributes(thisfile$path, paste0('/', which_site, '/Reflectance'))
    wv <- rhdf5::h5read(thisfile$path, paste0('/', which_site, '/Reflectance/Metadata/Spectral_Data'))
    
    omit_windows <- data.frame(
      omit_1_0 = md2$Band_Window_1_Nanometers[1], omit_1_f = md2$Band_Window_1_Nanometers[2],
      omit_2_0 = md2$Band_Window_2_Nanometers[1], omit_2_f = md2$Band_Window_2_Nanometers[2]
    )
    
    file_wv_df <- data.frame(year = as.numeric(thisfile$year), wv = round(wv$Wavelength), band = paste0('B', sprintf("%03d", 1:426))) %>%
      mutate(
        omit_band = (wv >= omit_windows$omit_1_0 & wv <= omit_windows$omit_1_f) | (wv >= omit_windows$omit_2_0 & wv <= omit_windows$omit_2_f),
        keep_band = (wv > 450) & (wv < 2150),
        data_ignore = md1$Data_Ignore_Value, SF = md1$Scale_Factor
      )
    kept_bands <- file_wv_df %>% filter(keep_band == TRUE & omit_band == FALSE)
    
    raw_data <- rhdf5::h5read(thisfile$path, paste0("/", which_site, "/Reflectance/Reflectance_Data"))
    reordered_data <- aperm(raw_data, c(2, 3, 1))
    
    map_info <- rhdf5::h5read(thisfile$path, paste0("/", which_site, "/Reflectance/Metadata/Coordinate_System/Map_Info"))
    site_crs <- unique(site_tiles$crs)
    hsi_rast_raw <- terra::rast(reordered_data, crs = paste0('EPSG:', site_crs))
    
    map_easting <- as.numeric(strsplit(map_info, ',')[[1]][4])
    map_northing <- as.numeric(strsplit(map_info, ',')[[1]][5]) #which is NOT the same as naming convention
    
    E_min <- map_easting
    E_max <- map_easting + 1000
    N_min <- map_northing - 1000
    N_max <- map_northing
    
    terra::ext(hsi_rast_raw) <- c(E_min, 
                                  E_max,
                                  N_min, 
                                  N_max)
    names(hsi_rast_raw) <- file_wv_df$band
    
    #clean raster
    hsi_rast_keep <- hsi_rast_raw[[unique(kept_bands$band)]] 
    hsi_rast_ignore <- subst(hsi_rast_keep, unique(kept_bands$data_ignore), NA) 
    hsi_rast <- hsi_rast_ignore / unique(kept_bands$SF)     
    
    thisfile_plots <- tile_plots %>% #created above
      filter(E_tile == thisfile$E_tile & N_tile == thisfile$N_tile)
    
    if (nrow(thisfile_plots) > 0) {
      for (t in 1:nrow(thisfile_plots)) {
        thisplot <- thisfile_plots[t,]
        thisplot_vect <- vect(thisplot, geom = c('easting', 'northing'), crs = crs(hsi_rast))
        square_plot <- as.polygons(ext(buffer(thisplot_vect, width = 10)))
        extracted_plot_data <- crop(hsi_rast, terra::ext(square_plot))
        
        extracted_df <- terra::as.data.frame(extracted_plot_data, xy = TRUE, cells = TRUE)
        rc <- terra::rowColFromCell(extracted_plot_data, extracted_df$cell)
        
        extracted_plot_df <- extracted_df %>%
          mutate(sample_row = rc[,1], sample_col = rc[,2], plotID = thisplot$plotID) %>%
          pivot_longer(cols = unique(kept_bands$band), names_to = 'band', values_to = 'refl') %>%
          left_join(kept_bands, by = "band") %>%
          left_join(thisplot %>% select(-any_of("year")), by = "plotID")
        
        all_plot_data[[paste0(thisplot$plotID, "_", thisfile$year, "_", f)]] <- extracted_plot_df
      }
    }
  } 
  rhdf5::h5closeAll()
  return(dplyr::bind_rows(all_plot_data))
}

map_site_tiles <- function(input) {
  if(nrow(input) == 0) return(list())
  
  site_tile_list <- list()
  
  for (f in 1:nrow(input)) {
    thisfile <- input[f,]
    
    file_is_readable <- tryCatch({
      md1 <- rhdf5::h5readAttributes(thisfile$path, paste0("/", which_site, "/Reflectance/Reflectance_Data"))
      TRUE
    }, error = function(e) {
      warning(paste("\nSkipping corrupted or inaccessible file:", thisfile$path))
      rhdf5::h5closeAll()
      FALSE
    })
    
    if (!file_is_readable) next
    
    wv <- rhdf5::h5read(thisfile$path, paste0('/', which_site, '/Reflectance/Metadata/Spectral_Data'))
    
    file_wv_df <- data.frame(wv = round(wv$Wavelength), band = paste0('B', sprintf("%03d", 1:426))) %>%
      mutate(
        data_ignore = md1$Data_Ignore_Value, SF = md1$Scale_Factor,
        special_band = case_when(
          wv == wv[which.min(abs(wv - 630))] ~ 'red', wv == wv[which.min(abs(wv - 800))] ~ 'NIR',
          wv == wv[which.min(abs(wv - 570))] ~ 'green', wv == wv[which.min(abs(wv - 480))] ~ 'blue',
          wv == wv[which.min(abs(wv - 531))] ~ 'PRI'
        )
      )
    
    special_bands <- file_wv_df %>% filter(!is.na(special_band)) %>% distinct()
    
    raw_data <- rhdf5::h5read(thisfile$path, paste0("/", which_site, "/Reflectance/Reflectance_Data"))
    reordered_data <- aperm(raw_data, c(2, 3, 1))
    
    map_info <- rhdf5::h5read(thisfile$path, paste0("/", which_site, "/Reflectance/Metadata/Coordinate_System/Map_Info"))
    site_crs <- unique(site_tiles$crs)
    hsi_rast_raw <- terra::rast(reordered_data, crs = paste0('EPSG:', site_crs))
    
    map_easting <- as.numeric(strsplit(map_info, ',')[[1]][4])
    map_northing <- as.numeric(strsplit(map_info, ',')[[1]][5])
    terra::ext(hsi_rast_raw) <- c(map_easting, map_easting + 1000, map_northing - 1000, map_northing)
    names(hsi_rast_raw) <- file_wv_df$band
    
    hsi_rast_keep <- hsi_rast_raw[[unique(special_bands$band)]] 
    hsi_rast_ignore <- subst(hsi_rast_keep, unique(special_bands$data_ignore), NA) 
    indices_subset <- hsi_rast_ignore / unique(special_bands$SF)     
    
    names(indices_subset) <- special_bands$special_band
    site_tile_list[[paste0(thisfile$year, "_", thisfile$E_tile, "_", thisfile$N_tile)]] <- indices_subset
  } 
  rhdf5::h5closeAll()
  return(site_tile_list)
}


input <- hsi_files_2024
foreach( 
  img = 1:length(input)) %do% {
    
    path <- input[[img]]
    
    # Safe read to avoid crashes from corrupted files
    file_is_readable <- tryCatch({
      md1 <- rhdf5::h5readAttributes(path, paste0("/", which_site, "/Reflectance/Reflectance_Data"))
      TRUE
    }, error = function(e) {
      warning(paste("\nSkipping corrupted or inaccessible file:", path))
      rhdf5::h5closeAll()
      FALSE
    })
    #if the file isn't readable go to the next one; if it is, it will continue.
    if (!file_is_readable) next
    
    if (file_is_readable) {
      
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
      #once a raster is made from raw data, clean up big df's to save memory:
      rm(raw_data, reordered_data)
      gc() # Force garbage collection
      
      #set extent and names of raster:
      terra::ext(hsi_rast_raw) <- c(E_min, 
                                    E_max,
                                    N_min, 
                                    N_max)
      names(hsi_rast_raw) <- file_wv_df$band
      #clean raster: remove atmospheric absorption windows; ignore values; divide by scale factor
      hsi_rast_keep <- hsi_rast_raw[[unique(kept_bands$band)]] 
      hsi_rast_ignore <- terra::subst(hsi_rast_keep, unique(kept_bands$data_ignore), NA) 
      hsi_rast <- hsi_rast_ignore / unique(kept_bands$SF)
      #export the list of regular rasters here:
      CPER_hsi_tile_list[[paste0(img)]] <- hsi_rast
      
      # ^this is the raster with ALL the reflectances in it, which should be used for extracting.
      
      # special bands raster:
      sp_bands_rast <- hsi_rast[[!is.na(kept_bands$special_band)]]
      names(sp_bands_rast) <- unique(kept_bands$special_band[!is.na(kept_bands$special_band)])
      #export the list of special band rasters here:
      CPER_spbands_tile_list[[paste0(img)]] <- sp_bands_rast
    }
  }

#RGB: 


#pulling out subplots from RGB imagery
subplots_rgb_list <- foreach (img = tilemap_rgb_list) %do% { 
  #individual img is a row from the tilemap rgb list (or whatever is entered above)
  tilemap <- img
  #establish the extent of the rgb tile so you can clip the list of all plots later to just those in the tile
  tile_extent <- terra::ext(tilemap)
  #project subplot shape file to the tile map so we can extract
  subplots_shp_aligned <- terra::project(which_subplots, tilemap)
  #filter all the subplots to just the ones in the tile's extent
  plots_in_tile <- subplots_shp_aligned[tile_extent]
  
  if (nrow(plots_in_tile) > 0) { #as long as there are plots in the tile,
    plot_lookup <- data.frame(
      ID = 1:nrow(plots_in_tile), #we will make a dataframe with plot ID (just a number)
      subplotID = plots_in_tile$subplotID, #the subplot and plot ID (from the dataframe)
      plotID = plots_in_tile$plotID #so that we can reference which plot is which
    )
    
    foreach(plot = nrow(plots_in_tile),
            .combine = rbind) %do% {
              subplot <- plots_in_tile[plot,] #for each plot, numbered 1 to nrow
              subplot_df <- as.data.frame(subplot)
              subplotID = subplot_df$subplotID
              buffer_pt <- terra::buffer(subplot, width = 2) #buffer so we account for spatial uncertainty?
              cropped_plot <- terra::crop(tilemap, buffer_pt) #crop the map to just the buffered point
              plotRGB(cropped_plot) #plot the cropped plot
            }
  } #plots in tile
} #each tile from the list

# Open files from list and extract necessary info, interpret the data
#input is a tile list: spbands or hsi (CPER_spbands_tile_list or CPER_hsi_tile_list)

#     Save summaries and HSI map
#output is summary and map list
# location info for plots:
#can use all_NEON_plots w LL / EN centroids or shapefiles, both are from NEON's site
all_NEON_plots <- read.csv(file = 'All_NEON_TOS_Plot_Centroids_V11.csv')
plots_shp <- terra::vect('/Users/khuelsma/Desktop/NEON Spectral Variability/Relevant NEON Materials/NEON TOS Plots/All_NEON_TOS_Plot_Polygons_V11.shp')
#CPER plots 
plots_shp_CPER_plants <- subset(plots_shp, stringr::str_detect(plots_shp$plotID, 'CPER') &
                                  'div' %in% plots_shp$appMods)
subplots_shp <- terra::vect('/Users/khuelsma/Desktop/NEON Spectral Variability/Relevant NEON Materials/NEON TOS Plots/All_NEON_TOS_Plot_Subplots_V11.shp')
#CPER subplots
subplot_shp_CPER_plants <- subset(subplots_shp, stringr::str_detect(subplots_shp$plotID, 'CPER') &
                                   'div' %in% subplots_shp$appMods)

  

#subplot shapes we can use to extract values
subplot_shp_400 <- subset(subplot_shp_CPER_plants, stringr::str_detect(subplot_shp_CPER_plants$subplotID, '_400'))
subplot_shp_100 <- subset(subplot_shp_CPER_plants, stringr::str_detect(subplot_shp_CPER_plants$subplotID, '_100'))
subplot_shp_10 <- subset(subplot_shp_CPER_plants, stringr::str_detect(subplot_shp_CPER_plants$subplotID, '_10_'))
subplot_shp_1 <- subset(subplot_shp_CPER_plants, stringr::str_detect(subplot_shp_CPER_plants$subplotID, '_1_'))



input <- CPER_hsi_tile_list
which_subplots <- subplot_shp_1
#which_subplots <- subplot_shp_10
#which_subplots <- subplot_shp_100
#which_subplots <- subplot_shp_400 #(or 100, 10 or 1)

subplots_hsi_1m <- foreach( 
  img = 1:length(input),
  .combine = rbind) %do% {
    
    hsi_rast <- input[[img]]
    #if you need to extract HSI: align hsi_rast to whatever you're using to extract (shape file or proj)
    tile_extent <- terra::ext(hsi_rast)
    subplots_shp_aligned <- terra::project(which_subplots, hsi_rast) #can adjust to any subplot input
    #restrict to just the subplots in that tile
    plots_in_tile <- subplots_shp_aligned[tile_extent]
    if (nrow(plots_in_tile) > 0) { #if there are plots in the tile:
      #create a lookup table with plotID and subplotID
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
      
      extracted_data_labeled
    }
  } #end each tile

extracted_1m <- subplots_hsi_1m
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

summary_1m %>%
  ggplot(aes(x = wv, y = mean_refl)) +
  geom_point()





#code that might be useful: 
#loading a single file to get relevant metadata; just in case?
which_site <- 'CPER'
which_year = 2024
path = hsi_files_2024[1]
raw_data <- rhdf5::h5read(hsi_files_2024[1], paste0("/", which_site, "/Reflectance/Reflectance_Data"))
reordered_data <- aperm(raw_data, c(3, 2, 1)) #make sure rows and columns are properly transposed in this step
epsg_code <- rhdf5::h5read(path, paste0("/", which_site, "/Reflectance/Metadata/Coordinate_System/EPSG Code"))
#needs wv df
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
hsi_rast_raw <- terra::rast(reordered_data, crs = paste0('EPSG:', epsg_code))
poly_proj <- terra::project(poly_vect, hsi_rast)
hsi_rast_keep <- hsi_rast_raw[[unique(kept_bands$band)]] 
hsi_rast_ignore <- terra::subst(hsi_rast_keep, unique(kept_bands$data_ignore), NA) 
hsi_rast <- hsi_rast_ignore / unique(kept_bands$SF)
tile_extent <- terra::ext(hsi_rast)
subplots_shp_aligned <- terra::project(subplot_shp_1, hsi_rast) #can adjust to any subplot input










