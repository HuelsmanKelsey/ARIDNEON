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
library('stringr')

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



# First see what's available ----------------------------------------------
get_prod_avail <- function(which_site, which_dp, BRDF = FALSE) {
  
  ID = case_when(
    which_dp == 'veg' ~ 'DP1.10058.001', #one and only veg
    which_dp == 'refl' & BRDF == TRUE ~ 'DP3.30006.002', #BRDF corrected; unusual case?
    which_dp == 'refl' ~ 'DP3.30006.001') #NOT BRDF corrected
  
  Info <- neonUtilities::getProductInfo(dpID = ID) 
  site_avail_data <- Info$siteCodes %>%
    filter(siteCode == which_site) %>%
    select(siteCode, availableMonths) %>%
    group_by(siteCode) %>%
    reframe(m = unlist(availableMonths)) %>%
    mutate(date_full = ymd(paste0(m, "-01")), #adds day so lubridate associates w date
           year = year(date_full),
           month = month(date_full)
    ) %>%
    distinct(siteCode, year, month) 
}

#just spec, just to see if different from both_avail (i.e., no concurrent veg)
spec_avail <- function(which_site) {
  
  nonbrdf <- get_prod_avail(which_site, 'refl', BRDF = FALSE) %>%
    mutate(brdf = 'no')
  
  brdf <- get_prod_avail(which_site, 'refl', BRDF = TRUE) %>%
    mutate(brdf = 'yes')
  
  both <- rbind(nonbrdf, brdf) %>%
    rename(specmonth = month)
}

#this will tell you which years have both veg + AOP data and whether the AOP data are brdf corrected
#NOTE: if the NEON API is down it won't work.
both_avail <- function(which_site) { #default is BRDF = FALSE
  veg <- get_prod_avail(which_site, 'veg')
  
  both <- spec_avail(which_site) %>%
    left_join(veg, by = c('siteCode', 'year'),
              relationship = 'many-to-many') %>%
    as.data.frame() %>%
    group_by(siteCode, year, specmonth, brdf) %>%
    summarise(vegmos = list(unique(month)))
  return(both)
}

# Downloading data --------------------------------------------------------
# function to get all easting and northing values for all tiles at *a site*
#req's locations of NEON plots from their website:
#can use all_NEON_plots w LL / EN centroids or shapefiles, both are from NEON's site
all_NEON_plots <- read.csv(file = 'All_NEON_TOS_Plot_Centroids_V11.csv')
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
home <- '/Users/khuelsma/'
avail_file_list <- function(which_site, 
                           which_year = NULL,
                           download = FALSE) {
  #tiles at the site
  site_tiles <- get_site_EN(which_site, tiles_only = TRUE)
  
  if(nrow(site_tiles) == 0) {
    warning("No plots found for site: ", which_site)
    #return(data.frame(site = character(), year = numeric(), path = character()))
  }

  #years at the site
  site_year_combos <- both_avail(which_site) 
  
  if (nrow(site_year_combos) == 0) {
    warning('No data found for site: ', which_site)
    # return(data.frame(siteCode = character(), 
    #                   year = numeric(), 
    #                   brdf = character(),
    #                   folder = cnumeric()))
  }
  
  #this will give us the number of bouts and therefore number of folders we need to find for each year
  minbouts <- site_year_combos %>%
    group_by(year) %>%
    summarise(minbouts = n_distinct(specmonth))
  
  site_year_tile_combos <- site_year_combos %>%
    left_join(minbouts) %>%
    cross_join(site_tiles) %>%
    rename(which_site = siteCode, 
           easting = E_tile, 
           northing = N_tile, 
           domain = domainID) %>%
    mutate(DP = ifelse(brdf == 'yes', 'DP3.30006.002', 'DP3.30006.001'),
           rgb_DP = 'DP3.30010.001') %>%
    select(which_site, year, minbouts, DP, brdf, easting, northing, domain, crs, rgb_DP)
  
  if (!is.null(which_year)) {
    site_year_tile_combos <- site_year_tile_combos %>%
      filter(year == which_year)
  } else { #all years
    site_year_tile_combos
  }
  #go thru each avail site file and see if we have it and can open it
  foreach(
    #for each row in site_year_tile_combos
    TYrow = 1:nrow(site_year_tile_combos),
    .combine = rbind
  ) %do% {
    tile_year_base <- site_year_tile_combos[TYrow,] %>%
      mutate(h5_filepath = NA,
             rgb_filepath = NA)
    rgb_DP = tile_year_base$rgb_DP
    DP = tile_year_base$DP
    year = tile_year_base$year 
    domain = tile_year_base$domain
    brdf = tile_year_base$brdf
    easting = tile_year_base$easting
    northing = tile_year_base$northing
    
    h5_base = paste0('NEON_', domain, '_', which_site, '_DP3_', easting, '_', northing)
    BRDF_filename_append = '_bidirectional_reflectance'
    h5filename = ifelse(brdf == 'yes', paste0(h5_base, BRDF_filename_append, '.h5'), paste0(h5_base, '_reflectance.h5'))
    
    # INNER LOOP
    found_tiles <- foreach(
      folder = 1:10, 
      .combine = rbind) %do% {
        
        # Make a fresh copy for this specific folder iteration
        current_tile <- tile_year_base 
        h5_path <- ifelse(
          brdf == 'yes', 
          paste0(home, DP, '/neon-aop-provisional-products/', year, '/FullSite/', domain, '/', year, '_', which_site, '_', folder, '/L3/Spectrometer/Reflectance/'),
          paste0(home, DP, '/neon-aop-products/', year, '/FullSite/', domain, '/', year, '_', which_site, '_', folder, '/L3/Spectrometer/Reflectance/')
        )
        h5_filepath = paste0(h5_path, h5filename)
        
        # 1. If file exists
        if (file.exists(h5_filepath)) {
          
          # safe read using tryCatch
          file_is_readable <- tryCatch({
            md1 <- rhdf5::h5readAttributes(h5_filepath, paste0("/", which_site, "/Reflectance/Reflectance_Data"))
            TRUE
          }, error = function(e) {
            warning(paste("\nSkipping corrupted or inaccessible file:", h5_filepath))
            rhdf5::h5closeAll()
            FALSE
          })
          
          # 2. If file is readable
          if (file_is_readable) { 
            current_tile$h5_filepath <- h5_filepath
            
            rgb_path = paste0(home, rgb_DP, '/neon-aop-products/', year, '/FullSite/', domain, '/', year, '_', which_site, '_', folder, '/L3/Camera/Mosaic/')
            rgbfilename = paste0(year, '_', which_site, '_', folder, '_', easting, '_', northing, '_image.tif')
            rgb_filepath_full = paste0(rgb_path, rgbfilename)
            
            # 3. Check RGB file
            if (file.exists(rgb_filepath_full)) {
              current_tile$rgb_filepath <- rgb_filepath_full
            } else { 
              current_tile$rgb_filepath <- NA
              print('download needed')
            }
            
            return(current_tile) # Return this valid row to rbind
            
          } else {
            return(NULL) # Unreadable, skip appending
          }
        } else {
          return(NULL) # Doesn't exist, skip appending
        }
      } # end inner foreach from folders
    
    # Return whatever valid rows were found for this tile/year
    found_tiles 
  } #each tile year
}


#make a wv dataframe if needed:
get_wvs <- function(which_site, which_year) {
  path <- avail_file_list(which_site, which_year)[1,]$h5_filepath
  md1 <- rhdf5::h5readAttributes(path, paste0("/", which_site, "/Reflectance/Reflectance_Data"))
  md2 <- rhdf5::h5readAttributes(path, paste0('/', which_site, '/Reflectance'))
  wv <- rhdf5::h5read(path, paste0('/', which_site, '/Reflectance/Metadata/Spectral_Data'))
  omit_windows <- data.frame(
    omit_1_0 = md2$Band_Window_1_Nanometers[1], omit_1_f = md2$Band_Window_1_Nanometers[2],
    omit_2_0 = md2$Band_Window_2_Nanometers[1], omit_2_f = md2$Band_Window_2_Nanometers[2]
  )
  # put metadata into wv df
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
  kept_bands <- kept_bands %>%
    select(year, wv, band)
  return(kept_bands)
}

home <- '/Users/khuelsma/'

#map_site_tiles will return a dataframe that gives you easy to open .tifs of full tiles: 
#using input of a list of downloaded tiles
# hsi, special bands, and rgb
map_site_tiles <- function(which_site, 
                           which_year = NULL, 
                           download = FALSE, 
                           output_dir = home) {
  
  filepath <- paste0(
    '/Users/khuelsma/Desktop/ARIDNEON/',
    which_site, '_', which_year, '_', 'tiles.csv')
  
  if (file.exists(filepath)) {
    
    input_with_paths <- read.csv(filepath)
    return(input_with_paths)
  }
  
  if (!file.exists(filepath)) {
    print('file does not exist')
  
  # 0. Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 1. Get the dataframe of available files
  input <- avail_file_list(which_site, which_year, download)
  
  if(nrow(input) == 0) {
    warning("No valid files found in avail_file_list.")
    return(input)
  }
  
  # 2. Iterate through rows, process rasters, SAVE to disk, and return paths
  input_with_paths <- foreach(
    img = 1:nrow(input), 
    .combine = rbind) %do% {
    
    this_row <- input[img,]
    path <- this_row$h5_filepath
    easting <- this_row$easting
    northing <- this_row$northing
    
    # Initialize output paths as NA
    hsi_out_path <- NA
    sp_out_path <- NA
    
    # --- LOAD H5 (HSI & SP BANDS) ---
    if (!is.na(path) && file.exists(path)) {
      
      # Define what the saved filenames should be
      hsi_target_file <- paste0(output_dir, which_site, "_", which_year, "_", easting, "_", northing, "_HSI.tif")
      sp_target_file <- paste0(output_dir, which_site, "_", which_year, "_", easting, "_", northing, "_SpBands.tif")
      
      # Check if we already processed this tile previously
      if (file.exists(hsi_target_file) && file.exists(sp_target_file)) {
        print(paste("Files exist, saving paths for tile:", easting, northing))
        hsi_out_path <- hsi_target_file
        sp_out_path <- sp_target_file
      } else { #if the tifs don't exist yet, make them:
        
        # Safe read to avoid crashes from corrupted files
        file_is_readable <- tryCatch({
          md1 <- rhdf5::h5readAttributes(path, paste0("/", which_site, "/Reflectance/Reflectance_Data"))
          TRUE
        }, error = function(e) {
          warning(paste("\nSkipping corrupted or inaccessible file:", path))
          rhdf5::h5closeAll()
          FALSE
        })
        
        if (file_is_readable) {
          
          md2 <- rhdf5::h5readAttributes(path, paste0('/', which_site, '/Reflectance'))
          wv <- rhdf5::h5read(path, paste0('/', which_site, '/Reflectance/Metadata/Spectral_Data'))
          omit_windows <- data.frame(
            omit_1_0 = md2$Band_Window_1_Nanometers[1], omit_1_f = md2$Band_Window_1_Nanometers[2],
            omit_2_0 = md2$Band_Window_2_Nanometers[1], omit_2_f = md2$Band_Window_2_Nanometers[2]
          )
          
          # put metadata into wv df
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
          reordered_data <- aperm(raw_data, c(3, 2, 1))
          epsg_code <- rhdf5::h5read(path, paste0("/", which_site, "/Reflectance/Metadata/Coordinate_System/EPSG Code"))
          map_info <- rhdf5::h5read(path, paste0("/", which_site, "/Reflectance/Metadata/Coordinate_System/Map_Info"))
          
          map_easting <- as.numeric(strsplit(map_info, ',')[[1]][4])
          map_northing <- as.numeric(strsplit(map_info, ',')[[1]][5])
          E_min <- map_easting
          E_max <- map_easting + 1000
          N_min <- map_northing - 1000
          N_max <- map_northing
          
          hsi_rast_raw <- terra::rast(reordered_data, crs = paste0('EPSG:', epsg_code))
          
          # Clean up big raw array immediately
          rm(raw_data, reordered_data)
          gc() 
          
          # set extent and names of raster
          terra::ext(hsi_rast_raw) <- c(E_min, E_max, N_min, N_max)
          names(hsi_rast_raw) <- file_wv_df$band
          
          # clean raster: remove atmospheric absorption windows; ignore values; divide by SF
          hsi_rast_keep <- hsi_rast_raw[[unique(kept_bands$band)]] 
          hsi_rast_ignore <- terra::subst(hsi_rast_keep, unique(kept_bands$data_ignore), NA) 
          hsi_rast <- hsi_rast_ignore / unique(kept_bands$SF)
          
          # special bands raster
          sp_bands_rast <- hsi_rast[[!is.na(kept_bands$special_band)]]
          names(sp_bands_rast) <- unique(kept_bands$special_band[!is.na(kept_bands$special_band)])
          
          # --- SAVE RASTERS TO DISK ---
          print(paste("Writing TIFs for tile:", easting, northing))
          terra::writeRaster(hsi_rast, hsi_target_file, overwrite = TRUE)
          terra::writeRaster(sp_bands_rast, sp_target_file, overwrite = TRUE)
          
          # Record the successful paths
          hsi_out_path <- hsi_target_file
          sp_out_path <- sp_target_file
          
          
          # Clean up raster memory from this iteration
          rm(hsi_rast_raw, hsi_rast_keep, hsi_rast_ignore, hsi_rast, sp_bands_rast)
          gc()
          rhdf5::h5closeAll()
          
        } # end if file is readable
      } # end if file exists check
    } # end if H5 path is valid
    
    # Return a 1-row data frame with the paths for this specific iteration
    saved_paths <-   data.frame(
      hsi_saved_path = hsi_out_path, 
      spbands_saved_path = sp_out_path, 
      stringsAsFactors = FALSE
    )
    
    cbind(this_row, saved_paths)
    
    } #each row 
  
  #export this dataframe:
  setwd('/Users/khuelsma/Desktop/ARIDNEON/')
  write.csv(input_with_paths, file = paste0(which_site, '_', which_year, '_', 'tiles.csv'))
  
  return(input_with_paths)
  } #this creates the file
}
input_with_paths <- map_site_tiles('CPER', 2024)

#IF YOU MADE IT THIS FAR, THE TILES ARE DOWNLOADED AND YOU CAN OPEN THEM

#the function version:
get_site_subplots <- function(which_site, which_size) {
  
  #all subplots (points)
  subplots_shp <- terra::vect(
    '/Users/khuelsma/Desktop/NEON Spectral Variability/Relevant NEON Materials/NEON TOS Plots/All_NEON_TOS_Plot_Subplots_V11.shp')

  #subset the shape using only plots at the site and "div" (diversity) applications
  site_subplots_shp <- subset(subplots_shp, 
                        stringr::str_detect(subplots_shp$plotID, which_site ) & 
                          stringr::str_detect(subplots_shp$appMods, 'div'))
  
  if (which_size == 100) {
    subs <- subset(site_subplots_shp, stringr::str_detect(site_subplots_shp$subplotID, paste0('_', which_size)))
    return(subs)
  }
  
  if (which_size == 1 | which_size == 10) {
    subs <- subset(site_subplots_shp, stringr::str_detect(site_subplots_shp$subplotID, paste0('_', which_size, '_')))
    return(subs)
  }
}

#function to open or to create (and save) a full site raster,
# which will be used in both extract_subplot_spectra and mapping!
create_site_raster <- function(which_rast, which_site, which_year) {
  input_dir = '/Users/khuelsma/Desktop/ARIDNEON/'
  csv_path <- paste0(input_dir, which_site, '_', which_year, '_tiles.csv')
  if(!file.exists(csv_path)) stop("Tiles CSV not found. Run map_site_tiles first.")
  input_with_paths <- read.csv(csv_path)
  
  merge_tiles <- function(filepaths) {
    valid_paths <- na.omit(unique(filepaths))
    rast_list <- lapply(valid_paths, terra::rast)
    return(terra::merge(terra::sprc(rast_list)))
  }

  if (which_rast == 'hsi') {
    # look for already created files or create them.
    # Import plots shape files to crop and save
    plots_shp <- terra::vect('/Users/khuelsma/Desktop/NEON Spectral Variability/Relevant NEON Materials/NEON TOS Plots/All_NEON_TOS_Plot_Polygons_V11.shp')
    site_plot_polygons <- subset(plots_shp, 
                                 stringr::str_detect(plots_shp$plotID, which_site ) & 
                                   stringr::str_detect(plots_shp$appMods, 'div'))
    # Project plots to match the raster CRS
    desired_projection <- terra::crs(terra::rast(na.omit(input_with_paths$hsi_saved_path)[1]))
    site_plots_projected <- terra::project(site_plot_polygons, desired_projection) #can also use David's data here
    
    #look for exported tile directory:
    hsi_dir <- '/Users/khuelsma/Desktop/ARIDNEON/hsi_plots/'
    
    if(!dir.exists(hsi_dir)) dir.create(hsi_dir, recursive = TRUE)

    expected_files <- paste0(hsi_dir, which_site, "_", which_year, "_", site_plots_projected$plotID, "_hsi.tif")
    # If ANY of the expected plots are missing, we run the extraction
    if (!all(file.exists(expected_files))) {
      print('Missing cached HSI plots. Generating them now...')
      hsi_merged <- merge_tiles(input_with_paths$hsi_saved_path)
      
      #open the files:
      for (ID in 1:nrow(site_plots_projected)) {
        this_plot <- site_plots_projected[ID, ]
        plotID <- as.character(this_plot$plotID)
        plot_ext <- terra::ext(this_plot)
        
        #doing tryCatch so if there's an error it's chill; if not null, continues
        cropped_plot <- tryCatch({
          terra::crop(hsi_merged, plot_ext)
        }, error = function(e) NULL)
        
        if (!is.null(cropped_plot)) {
          #so when you need to open: 
          file_name <- paste0(hsi_dir, which_site, "_", which_year, "_", plotID, "_hsi.tif")
          terra::writeRaster(cropped_plot, file_name, overwrite = TRUE)
          }
      }
    }
    valid_files <- expected_files[file.exists(expected_files)]
    #return a virtual map of the plots
    return(terra::vrt(valid_files))
  } #if which_rast == hsi
  
  if (which_rast == 'rgb') {
    save_dir <- '/Users/khuelsma/Desktop/ARIDNEON/'
    file_name <- paste0(save_dir, which_site, "_", which_year, "_rgb.tif")
    if (file.exists(file_name)) {
      return(terra::rast(file_name))
    } else {
      print('Missing cached RGB')
      rgb_merged <- merge_tiles(input_with_paths$rgb_filepath)
      terra::writeRaster(rgb_merged, file_name, overwrite = TRUE)
      return(rgb_merged)
    }
  }
  if (which_rast == 'vi') {
    save_dir <- '/Users/khuelsma/Desktop/ARIDNEON/'
    file_name <- paste0(save_dir, which_site, "_", which_year, "_vi.tif")
    if (file.exists(file_name)) {
      return(terra::rast(file_name))
    } else {
      print('Missing cached vi map')
      # Merge SpBands first, then calculate VIs globally across the whole site!
      sp_merged <- merge_tiles(input_with_paths$spbands_saved_path)
      
      # Calculate VIs on the merged raster (much faster/cleaner than tile-by-tile)
      NDVI <- (sp_merged[['NIR']] - sp_merged[['red']]) / (sp_merged[['NIR']] + sp_merged[['red']] + 0.0001)
      PRI  <- (sp_merged[['NIR']] - sp_merged[['green']]) / (sp_merged[['NIR']] + sp_merged[['green']])
      NIRv <- sp_merged[['NIR']] * NDVI
      CCI  <- (sp_merged[['PRI']] - sp_merged[['red']]) / (sp_merged[['PRI']] + sp_merged[['red']])
      
      vi_merged <- c(CCI, NIRv, PRI)
      names(vi_merged) <- c("CCI", "NIRv", "PRI")
      terra::writeRaster(vi_merged, file_name, overwrite = TRUE)
      return(vi_merged)
    }
  }
}

CPER_hsi_2021 <- create_site_raster('hsi', 'CPER', 2021)
CPER_hsi_2024 <- create_site_raster('hsi', 'CPER', 2024)

create_square_polygons <- function(input) { #where input is the projected subplots
  coords <- terra::crds(input)
  
  poly_list <- lapply(1:nrow(coords), function(i) {
    x <- coords[i, "x"]
    y <- coords[i, "y"]
    side <- sqrt(input$subpltSize[i])
    if (is.na(side) || is.null(side)) side <- 20 # Fallback for 400m2
  # Create extent (xmin, xmax, ymin, ymax) and cast to polygon
  sq_ext <- terra::ext(x, x + side, y, y + side)
  sq_poly <- terra::vect(sq_ext)
  terra::crs(sq_poly) <- terra::crs(input)
  return(sq_poly)
  }) 
  return(do.call(rbind, poly_list))
}

#using the tile list, we can make a full site map or extract plots and subplots.
#so if I say make_maps, I want to load every tile, but extract the plot [[plot_ID]] and then subplots [[subplot_ID]]
#for make_maps, save a list of them. For extract_subplot_spectra, 
make_maps <- function(which_rast, 
                      which_site, 
                      which_year, 
                      full_site = FALSE,
                      which_buff) {
  
  site_raster <- create_site_raster(which_rast, 
                                    which_site, 
                                    which_year)
  if(is.null(site_raster)) stop("site_raster not found. Run create_site_raster(which_rast, which_site, which_year) first.")
  
  # Import plots shape files
  plots_shp <- terra::vect('/Users/khuelsma/Desktop/NEON Spectral Variability/Relevant NEON Materials/NEON TOS Plots/All_NEON_TOS_Plot_Polygons_V11.shp')
  site_plot_polygons <- subset(plots_shp, 
                               stringr::str_detect(plots_shp$plotID, which_site ) & 
                                 stringr::str_detect(plots_shp$appMods, 'div'))
  # Project plots to match the raster CRS
  site_plots_projected <- terra::project(site_plot_polygons, terra::crs(site_raster)) #can also use David's data here
  
  if (full_site == TRUE) {
    terra::plotRGB(site_raster, r = 1, g = 2, b = 3, stretch = 'lin')
    terra::plot(site_plots_projected, add = TRUE, border = "yellow", lwd = 4)
    terra::text(site_plots_projected, site_plots_projected$plotID, col = "black", halo = TRUE, hc = "white", cex = 1.2, font = 2, pos = 3)
    
    return(site_raster)
  }
  
  #moved get subplots
  site_subplots <- get_site_subplots(which_site, 1)
  subplots_projected <- terra::project(site_subplots, site_raster)
  
  site_output_list <- list() # Initialize the primary named list
  
  # Loop over every PLOT in the site
  for (ID in 1:nrow(site_plots_projected)) {
    this_plot <- site_plots_projected[ID, ]
    plotID <- as.character(this_plot$plotID)
    plot_ext <- terra::ext(this_plot)
    
    #doing tryCatch so if there's an error it's chill; if not null, continues
    cropped_plot <- tryCatch({
      terra::crop(site_raster, plot_ext)
    }, error = function(e) NULL)
    
    if (!is.null(cropped_plot)) {
      save_dir <- '/Users/khuelsma/Desktop/ARIDNEON/hsi_plots/'
      if(!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
      
      file_name <- paste0(save_dir, which_site, "_", which_year, "_", plotID, "_", which_rast, ".tif")
      
      if(!file.exists(file_name)) { #only if the file doesn't already exist, make a raster of it
        terra::writeRaster(cropped_plot, file_name, overwrite = TRUE)
      }
      
      # Initialize a nested list for THIS specific plot
      plot_list <- list()
      plot_list[['plot_full']] <- cropped_plot
      
      # Find all subplots that belong to this plot
      subplots_in_plot <- subset(subplots_projected, subplots_projected$plotID == plotID)
      
      # If subplots exist, extract them
      if (nrow(subplots_in_plot) > 0) {
        
        #create square geometries for subplots:
        sub_polys <- create_square_polygons(subplots_in_plot)
        #for each one,
        for (sub in 1:nrow(sub_polys)) {
          this_sub <- subplots_in_plot[sub, ]
          subplotID <- as.character(this_sub$subplotID)
          
          # Buffer and crop from the already-cropped plot raster for speed
          buff_ext <- terra::ext(terra::buffer(sub_polys[sub, ], which_buff))
          cropped_subplot <- tryCatch({ terra::crop(cropped_plot, buff_ext) }, error = function(e) NULL)
          
          if (!is.null(cropped_subplot)) plot_list[[subplotID]] <- cropped_subplot
        }
      }
      site_output_list[[plotID]] <- plot_list
    }
  }
  return(site_output_list)
}


extract_subplot_spectra <- function(which_rast, which_site, which_year, which_size, which_buff) {
 # replacing input_with_paths <- map_site_tiles(which_site, which_year) with create_site_raster outputs
  site_raster <- create_site_raster(which_rast, which_site, which_year)
  #input <- raster_list
  which_subplots <- get_site_subplots(which_site, which_size)
  #project to site_raster projection
  subplots_proj <- terra::project(which_subplots, terra::crs(site_raster))
  
  #create polygons using create_square_polygons helper and buffer.
  subplot_polys <- create_square_polygons(subplots_proj)
  buffed_subplots <- terra::buffer(subplot_polys, which_buff)
  
  extracted_spectra <- data.frame() #replacing c() with data.frame()

  for (subplot in 1:nrow(buffed_subplots)) {
    this_sub <- buffed_subplots[subplot, ]
    this_subplot_df <- as.data.frame(this_sub)
  
    subplot_metrics <- terra::extract(site_raster, 
                                        this_sub,
                                        xy = TRUE,
                                        cells = TRUE) %>%
        cbind(this_subplot_df) %>% #changed cross_join to cbind
      
        #making up bout number for this one... will need to do something diff for 2020 data though
      
        mutate(eventID = paste0(plotID, '_', subplotID, '_', which_year, '_', 1))
      
    extracted_spectra <- rbind(extracted_spectra, subplot_metrics)
    }
       
  return(extracted_spectra)
}


CPER_2021_hsi <- create_site_raster('hsi', 'CPER', 2021)
CPER_2024_hsi <- create_site_raster('hsi', 'CPER', 2024)
CPER_2021_rgb <- create_site_raster('rgb', 'CPER', 2021)
CPER_2024_rgb <- create_site_raster('rgb', 'CPER', 2024)
CPER_2021_vi <- create_site_raster('vi', 'CPER', 2021)
CPER_2024_vi <- create_site_raster('vi', 'CPER', 2024)

extract_subplot_spectra('hsi', 'CPER', 2024, 1, 0)
