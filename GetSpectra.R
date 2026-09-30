# Overview: ---------------------------------------------------------------

# Get Spectra lets you

# Choose a site/location and open or request from API to download; 
#     Ends w check files exist and outputs a list of files
# open files from list and extract necessary info,
#     Save summaries and map: full site vegetation indices and 10 cm rgb; plot level (big plots) snapshots of rgb; and plot level (big plots) snapshots of 1m hsi

neon_token <- if (file.exists('/Users/khuelsma/NEON_token.txt')) {
  neon_token <- readLines('/Users/khuelsma/NEON_token.txt')
} else {
  neon_token <- rstudioapi::askForPassword(prompt = 'enter NEON token')
}
options(timeout = 3600) #increase timeout to 1 hour
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
library('purrr')


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
repo_dir <- '/Users/khuelsma/Desktop/ARIDNEON/'
#need to add these shapefiles to repo
DA_vegpt <-  "/Users/khuelsma/ARIDNEON_Working/Augustine_VegPoint_2026.shp/VegPOINT_2026.shp"
DA_vegpoly <- "/Users/khuelsma/ARIDNEON_Working/Augustine_VegPolygon_2026/VegPolygon2026.shp"
other_locations <- list(DA_vegpt, DA_vegpoly)


#get site EN gets NEON's veg plot locations for a given site;
# default for other locations = NULL, but you can add them by setting above.
get_site_EN <- function(which_site, 
                        other_locations = NULL, #default
                        tiles_only = FALSE) {   #this is for file checking
  
  #filter NEON's plot centroids from the repo to site, get crs for the cases when we have LL inputs
  neon_plot_polygons <- read.csv(file = paste0(repo_dir, 'All_NEON_TOS_Plot_Centroids_V11.csv')) %>%
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
    select(plotID, easting, northing, E_tile, N_tile, domainID, crs)
  
  target_crs <- unique(neon_plot_polygons$crs)
  domainID <- unique(neon_plot_polygons$domainID)
  
  plot_polygons <- neon_plot_polygons %>% 
    ungroup()

  #other_locations is set to NULL by default, but you can add other_locations = whatever list
  
  if (!is.null(other_locations)) { #not default
    for (this_geom in other_locations) {  #go through each set of location files
      input_geom <- sf::st_zm(sf::st_read(this_geom, quiet = TRUE))
      geom_type <- as.character(sf::st_geometry_type(input_geom, by_geometry = FALSE))
      geom_proj <- terra::project(terra::vect(input_geom), paste0("EPSG:", target_crs))

      other_tiles <- terra::geom(geom_proj) %>%
        as.data.frame() %>%
        mutate(
          this_row = row_number(),
          plotID = paste0('other_DA_', 
                          geom_type, '_',
                          this_row),
          easting = x, northing = y, 
          E_tile = 1000*floor(x/1000),
          N_tile = 1000*floor(y/1000),
          domainID = domainID,
          crs = target_crs) %>%
        select(plotID, easting, northing, E_tile, N_tile, domainID, crs)
      plot_polygons <- rbind(plot_polygons, other_tiles) %>% distinct()
    }
  } #other locations == TRUE
  
  #if tiles only (for file checking and opening)
  if (tiles_only == TRUE) {
    return(plot_polygons %>% distinct(E_tile, N_tile, domainID, crs))
  }
  # if want details about each plot, tiles_only == FALSE, which is default
  return(plot_polygons)
}
all_plots <- get_site_EN('CPER', other_locations = other_locations)


avail_file_list <- function(which_site, 
                           which_year = NULL,
                           other_locations = NULL, # default; inherited from get_site_EN
                           download = FALSE,
                           home_dir = '/Users/khuelsma/',
                           repo_dir = "/Users/khuelsma/Desktop/ARIDNEON/"
                           ) {
  #tiles at the site
  site_tiles <- get_site_EN(which_site, other_locations, tiles_only = TRUE)
  if (nrow(site_tiles) ==0) return(data.frame())
  #years at the site
  site_year_combos <- both_avail(which_site) 
  if (nrow(site_year_combos) == 0) return(data.frame())

  #this will give us the number of bouts and therefore number of folders we need to find for each year
  minbouts <- site_year_combos %>%
    group_by(year) %>%
    summarise(minbouts = n_distinct(specmonth))
  
  #get any years and tiles with observations in the same year:
  site_year_tile_combos <- site_year_combos %>%
    left_join(minbouts, by = "year") %>%
    cross_join(site_tiles) %>%
    rename(siteID = siteCode, 
           easting = E_tile, 
           northing = N_tile, 
           domain = domainID) %>%
    mutate(DP = ifelse(brdf == 'yes', 'DP3.30006.002', 'DP3.30006.001'),
           rgb_DP = 'DP3.30010.001') %>%
    distinct(siteID, year, minbouts, DP, brdf, easting, northing, domain, crs, rgb_DP)
  
  #if you select a year, it will filter, otherwise just defaults to all (which_year == NULL)
  if (!is.null(which_year)) {
    site_year_tile_combos <- site_year_tile_combos %>% filter(year == which_year)
  }
  if (nrow(site_year_tile_combos) == 0) {
    warning(paste("No matching files for", which_site, which_year))
    return(data.frame())
  }
  
  #go thru each avail site file and see if we have it and can open it
  final_list <- foreach(
    TYrow = 1:nrow(site_year_tile_combos), 
    .combine = rbind) %do% {
      
      current_tile <- site_year_tile_combos[TYrow,]

      h5_base <- paste0('NEON_', current_tile$domain, '_', which_site, '_DP3_', current_tile$easting, '_', current_tile$northing)
      BRDF_append <- '_bidirectional_reflectance'
      h5filename <- ifelse(current_tile$brdf == 'yes', 
                           paste0(h5_base, BRDF_append, '.h5'), 
                           paste0(h5_base, '_reflectance.h5'))
             
      #helper function, find_local_files, which establishes an empty list for h5, rgb, and folder number(s)             
      find_local_files <- function() {
        file_list <- list(h5 = NA, 
                          rgb = NA, 
                          found_folder = NA)
        
        for (folder in 1:10) { #for 1:10, look for the file. If it exists, that's the folder number
          
          # Check H5
          h5_path <- ifelse(current_tile$brdf == 'yes', 
                            paste0(home_dir, current_tile$DP, '/neon-aop-provisional-products/', current_tile$year, '/FullSite/', current_tile$domain, '/', current_tile$year, '_', which_site, '_', folder, '/L3/Spectrometer/Reflectance/'),
                            paste0(home_dir, current_tile$DP, '/neon-aop-products/', current_tile$year, '/FullSite/', current_tile$domain, '/', current_tile$year, '_', which_site, '_', folder, '/L3/Spectrometer/Reflectance/'))
          
          test_h5 <- paste0(h5_path, h5filename)
          
          if (is.na(file_list$h5) && file.exists(test_h5)) {
            file_list$found_folder <- folder #save the folder number regardless of if it's readable.
            
            #see if the file is readable, and if it is, add it to the file list.
            is_readable <- tryCatch({
              rhdf5::h5readAttributes(test_h5, paste0("/", which_site, "/Reflectance/Reflectance_Data"))
              rhdf5::h5closeAll() # close after testing
              TRUE
            }, error = function(e) {
              rhdf5::h5closeAll()
              FALSE
            })
            
            if (is_readable) file_list$h5 <- test_h5
          } # check for h5 filename 
          
          } #all folders
        
        # Check RGB independently, but use the folder optimization structure created above.
                            # If we found the H5 folder, only check that one. Otherwise, check all 10.
                            
        folders_to_check <- if (!is.na(file_list$found_folder)) unique(file_list$found_folder) else 1:10

        for (folder in folders_to_check) {
          rgb_path <- paste0(home_dir, current_tile$rgb_DP, '/neon-aop-products/', current_tile$year, '/FullSite/', current_tile$domain, '/', current_tile$year, '_', which_site, '_', folder, '/L3/Camera/Mosaic/')
          rgbfilename <- paste0(current_tile$year, '_', which_site, '_', folder, '_', current_tile$easting, '_', current_tile$northing, '_image.tif')
          test_rgb <- paste0(rgb_path, rgbfilename)
          #if test_rgb exists, put it in the list
          if (is.na(file_list$rgb) && file.exists(test_rgb)) {
            file_list$rgb <- test_rgb
            }
          } #for all folders
        return(file_list)
      } #found local files 
      
      #run the function created above to find all local files; 
      found <- find_local_files()
      download_success <- FALSE #set download success as false ahead of determining needs to download
      #because we will only check for files again if download_success switches to TRUE
      
      if (is.na(found$h5) && download == TRUE) {
        tryCatch({
          neonUtilities::byTileAOP(dpID = current_tile$DP, 
                                   site = which_site, 
                                   easting = current_tile$easting, 
                                   northing = current_tile$northing, 
                                   buffer = 0, check.size = FALSE, 
                                   include.provisional = TRUE, 
                                   year = current_tile$year, 
                                   token = neon_token, 
                                   savepath = home_dir, 
                                   progress = TRUE)
          download_success <- TRUE 
        }, error = function(e) warning("HDF5 Download failed for tile: ", current_tile$easting, "_", current_tile$northing))
      }
      
      if (is.na(found$rgb) && download == TRUE) {
        tryCatch({
          neonUtilities::byTileAOP(dpID = current_tile$rgb_DP, 
                                   site = which_site, 
                                   easting = current_tile$easting, 
                                   northing = current_tile$northing, 
                                   buffer = 0, check.size = FALSE, 
                                   include.provisional = TRUE, 
                                   year = current_tile$year, 
                                   token = neon_token, 
                                   savepath = home_dir, 
                                   progress = TRUE)
          download_success <- TRUE 
        }, error = function(e) warning("RGB Download failed for tile: ", current_tile$easting, "_", current_tile$northing))
      }
      
      if (download_success) { #should only have download_success if needed to download
        found <- find_local_files() 
      }
      current_tile$h5_filepath <- found$h5
      current_tile$rgb_filepath <- found$rgb
      current_tile$folderNumber <- found$found_folder 
      current_tile
    } #for each tile / year combo row
  
  final_list <- final_list %>% 
    mutate(boutNumber = match(folderNumber, unique(final_list$folderNumber)))
  return(final_list)
}

file_list <- avail_file_list('CPER', 2020, other_locations = other_locations, download = TRUE)

#make a wv dataframe if needed:#make a wv dataframe if neededother_locations = :
get_wvs <- function(which_site, which_year) { #which_site is included, because it will look for a filepath for reference
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
    select(year, wv, band, special_band)
  return(kept_bands)
}

# input = rgb
# which_year requires an input
# other locations defaults to NULL, meaning only NEON plots are used
# download = FALSE is default
get_site_raster <- function(which_site, 
                           which_year, #no more NULL default 
                           which_bout = 1, #default
                           other_locations = NULL,
                           download = FALSE,  #default
                           home_dir = '/Users/khuelsma/',
                           repo_dir = "/Users/khuelsma/Desktop/ARIDNEON/") { #default
  
  if (is.null(which_year)) stop('specify a year')
  print(paste("Building Virtual RGB Site Map for", which_year, "..."))
  input <- avail_file_list(which_site,
                           which_year, 
                           other_locations,
                           download,
                           home_dir, 
                           repo_dir)
  
  if(nrow(input) == 0) {
    warning("No valid files found")
    return(NULL)
  }
  
  file_input <- input %>%
    filter(boutNumber == which_bout)
    
  paths <- unique(file_input$rgb_filepath) #may need to make sure the files exist
    if (length(paths) == 0) return(NULL)
    site_raster <- terra::vrt(paths) #rgb is already tiffs
    names(site_raster) <- c('red', 'green', 'blue')
    return(site_raster)
}

rgb_rast <- get_site_raster('CPER', 2021, other_locations = other_locations)

cache_hsi_plots <- function(which_site, 
                            which_year, 
                            which_bout = 1, #default
                            other_locations = NULL,
                            download = FALSE,
                            home_dir = '/Users/khuelsma/',
                            repo_dir = "/Users/khuelsma/Desktop/ARIDNEON/") {
  
  save_dir <- paste0(repo_dir, 'plots/')
  if(!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  
  # 1. Prepare Plots
  NEON_plots_shp <- terra::vect('/Users/khuelsma/Desktop/NEON Spectral Variability/Relevant NEON Materials/NEON TOS Plots/All_NEON_TOS_Plot_Polygons_V11.shp')
  NEON_plots <- subset(NEON_plots_shp, stringr::str_detect(NEON_plots_shp$plotID, which_site) & stringr::str_detect(NEON_plots_shp$appMods, 'div'))
  plot_inputs <- list(NEON_plots)
  
  if (!is.null(other_locations)) {
    
    for (l in 1:length(other_locations)) {
      sf_obj <- sf::st_zm(sf::st_read(other_locations[[l]], quiet = TRUE))
      geom_type <- as.character(sf::st_geometry_type(sf_obj, by_geometry = FALSE))
      
      other_plots <- terra::vect(sf_obj)
      other_plots$row <- 1:nrow(other_plots)
      other_plots$plotID <- paste0('other_DA_', 
                                   geom_type, '_',
                                   1:nrow(other_plots))
      #can't rbind these, leave them in a list: 
      plot_inputs[[length(plot_inputs) + 1]] <- other_plots[, "plotID"]
      
    }
  } #so plots_inputs is a list of multiple spat vectors
  
  # 2. Check which plots are missing from saved files
  all_plotIDs <- unlist(lapply(plot_inputs, function(x) x$plotID))
  
  file_list <- avail_file_list(which_site, which_year, other_locations, download, home_dir, repo_dir)
  
  if(nrow(file_list) == 0) {
    warning("No valid files found")
    return(NULL)
  }
  
  file_input <- file_list %>%
    filter(boutNumber == which_bout) %>%
    print()
  
  #by default, it won't say a bout number; if there are multiple boutnumbers, it will include before _hsi.tif
  expected_files <- ifelse(which_bout == 1, 
                           paste0(save_dir, which_site, "_", which_year, "_", all_plotIDs, "_hsi.tif"),
                           paste0(save_dir, which_site, "_", which_year, "_", all_plotIDs, "_", which_bout, "_hsi.tif"))
  #was expected_files <- paste0(save_dir, which_site, "_", which_year, "_", all_plotIDs, "_hsi.tif")
  if (all(file.exists(expected_files))) {
    print("All HSI plots already cached! Skipping H5 processing.")
    return(TRUE)
  }
  # 3. If missing plots, use the available files list
  if (nrow(file_input) == 0) stop("No files available to process.")
  
  # 4. Process ONLY the tiles we need
  for (img in 1:nrow(file_input)) {
    h5_path <- file_input$h5_filepath[img]
    if (is.na(h5_path) || !file.exists(h5_path)) next
    
    easting <- file_input$easting[img]
    northing <- file_input$northing[img]
    
    bout <- file_input$boutNumber[img]
    
    # Peek at the EPSG code first so we can project our plots correctly!
    epsg_code <- tryCatch({ 
      rhdf5::h5read(h5_path, paste0("/", which_site, "/Reflectance/Metadata/Coordinate_System/EPSG Code")) }, 
      error = function(e) NULL)
    if (is.null(epsg_code)) next
    
    tile_crs <- paste0('EPSG:', epsg_code)
    tile_poly <- terra::vect(terra::ext(easting, easting + 1000, northing - 1000, northing))
    
    # Loop through our inputs list to find ANY intersecting plots
    plots_to_crop <- list()
    
    for (i in 1:length(plot_inputs)) {
      proj_input <- terra::project(plot_inputs[[i]], tile_crs)
      intersects_logical <- as.vector(terra::relate(proj_input, tile_poly, "intersects"))
      intersecting <- proj_input[intersects_logical, ]
      
      if (nrow(intersecting) > 0) {
        plots_to_crop[[length(plots_to_crop) + 1]] <- intersecting
      }
    }
    
    if (length(plots_to_crop) == 0) next # Skip this massive H5 file entirely if no plots are in it!
    
    # If there are plots, open the H5
    file_is_readable <- tryCatch({ 
      rhdf5::h5readAttributes(h5_path, paste0("/", which_site, "/Reflectance/Reflectance_Data")); TRUE }, 
      error = function(e) { rhdf5::h5closeAll(); FALSE })
    
    if (file_is_readable) {
      print(paste("Extracting plots from tile:", easting, northing))
      
      md1 <- rhdf5::h5readAttributes(h5_path, paste0("/", which_site, "/Reflectance/Reflectance_Data"))
      md2 <- rhdf5::h5readAttributes(h5_path, paste0('/', which_site, '/Reflectance'))
      wv <- rhdf5::h5read(h5_path, paste0('/', which_site, '/Reflectance/Metadata/Spectral_Data'))
      
      omit_windows <- data.frame(
        omit_1_0 = md2$Band_Window_1_Nanometers[1], omit_1_f = md2$Band_Window_1_Nanometers[2],
        omit_2_0 = md2$Band_Window_2_Nanometers[1], omit_2_f = md2$Band_Window_2_Nanometers[2]
      ) 
      
      file_wv_df <- data.frame(wv = round(wv$Wavelength), band = paste0('B', sprintf("%03d", 1:426))) %>%
        mutate(
          omit_band = (wv >= omit_windows$omit_1_0 & wv <= omit_windows$omit_1_f) | (wv >= omit_windows$omit_2_0 & wv <= omit_windows$omit_2_f),
          keep_band = (wv > 450) & (wv < 2150)
        )
      kept_bands <- file_wv_df %>% filter(keep_band == TRUE & omit_band == FALSE)
      
      raw_data <- rhdf5::h5read(h5_path, paste0("/", which_site, "/Reflectance/Reflectance_Data"))
      reordered_data <- aperm(raw_data, c(3, 2, 1))
      epsg_code <- rhdf5::h5read(h5_path, paste0("/", which_site, "/Reflectance/Metadata/Coordinate_System/EPSG Code"))
      
      # Keep the raster RAW (No math yet to save memory and temp space!)
      hsi_rast_raw <- terra::rast(reordered_data, crs = paste0('EPSG:', epsg_code))
      terra::ext(hsi_rast_raw) <- terra::ext(tile_poly)
      
      rm(raw_data, reordered_data)
      gc()
      # Loop through the successful intersecting plots and crop them!
      for (i in 1:length(plots_to_crop)) {
        this_group <- plots_to_crop[[i]]
        
        for (p in 1:nrow(this_group)) {
          this_plot <- this_group[p, ]
          plotID <- as.character(this_plot$plotID)
          
          #by default, it won't say a bout number; if there are multiple boutnumbers, it will include before _hsi.tif
          file_name <- ifelse(bout == 1, 
                              paste0(save_dir, which_site, "_", which_year, "_", plotID, "_hsi.tif"),
                              paste0(save_dir, which_site, "_", which_year, "_", plotID, "_", boutNumber, "_hsi.tif"))

          if (!file.exists(file_name)) {
            cropped_raw <- tryCatch({ 
              terra::crop(hsi_rast_raw, terra::ext(this_plot)) }, error = function(e) NULL)
            
            if (!is.null(cropped_raw)) {
              names(cropped_raw) <- file_wv_df$band
              cropped_keep <- cropped_raw[[unique(kept_bands$band)]]
              cropped_clean <- terra::subst(cropped_keep, md1$Data_Ignore_Value, NA) / md1$Scale_Factor
              terra::writeRaster(cropped_clean, file_name, overwrite = TRUE)
            }
          }
        }
      } #plots to crop
      
      rm(hsi_rast_raw)
      gc()
      rhdf5::h5closeAll()
    } #file readable
  } #each tile
  return(TRUE) # This just tells extract_from_maplist() that the caching succeeded!
}

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

#make maps will create maps for hsi or rgb (with different settings)
# choose a site, a year, 
# whether you want to include other_locations
# which size and which_buff
# download is set to FALSE by default, 
# full_site is set to FALSE by default
make_maps <- function(which_site, 
                      which_year, #no default since maps need to coincide temporally
                      
                      which_bout = 1, #default, but will need which_bout for 2020
                      other_locations = NULL, 
                      download = FALSE, #default, but can set to TRUE which will trigger site vrt
                      
                      full_site = FALSE,
                      
                      which_size,
                      which_buff,
                      
                      home_dir = "/Users/khuelsma/",
                      repo_dir = "/Users/khuelsma/Desktop/ARIDNEON/") {
  
  if (is.null(which_year)) stop('specify a year to make maps')
  
  save_dir <- paste0(repo_dir,'plots/')
  if(!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  
  # create names for these rasters
  rgb_names <- c('red', 'green', 'blue')
  
  # Import sampling plots shape files:
  NEON_plots_shp <- terra::vect('/Users/khuelsma/Desktop/NEON Spectral Variability/Relevant NEON Materials/NEON TOS Plots/All_NEON_TOS_Plot_Polygons_V11.shp')
  NEON_site_plot_polygons <- subset(NEON_plots_shp, 
                                    stringr::str_detect(NEON_plots_shp$plotID, which_site) & 
                                      stringr::str_detect(NEON_plots_shp$appMods, 'div'))
  
  shapefile_plots <- list()  #create an empty list for other locations
  if (!is.null(other_locations)) { #if other_locations isn't NULL (NULL is default, so you have to add them)
    for (l in 1:length(other_locations)) {
      #read the shape file
      sf_obj <- sf::st_zm(sf::st_read(other_locations[[l]], quiet = TRUE))
      geom_type <- as.character(sf::st_geometry_type(sf_obj, by_geometry = FALSE))
      terra_obj <- terra::vect(sf_obj) 
      terra_obj$row <- 1:nrow(terra_obj)
      terra_obj$plotID <- paste0('other_DA_', 
                                 geom_type, '_',
                                 1:nrow(terra_obj))
      
      # Check geometry type on each; if it's a POINT we need to turn it into a polygon for cropping
      if (str_detect(tolower(geom_type), 'point')) {
        for (i in 1:nrow(terra_obj)) {
          this_poly <- terra_obj[i,]
          this_polyid <- as.character(this_poly$plotID) 
          # establish extent, then turn it into a vector / polygon
          poly_from_pt_ext <- terra::ext(terra::buffer(this_poly, which_buff))
          poly_from_pt <- terra::vect(poly_from_pt_ext)
        } #for each row in the shapefile, turn it into a polygon
        } #if it's a point
      #can't rbind these, leave them in a list: 
      shapefile_plots[[l]] <- terra_obj
    }
  }
  
  site_raster <- get_site_raster(which_site, which_year, which_bout, other_locations, download, home_dir, repo_dir)
  if (is.null(site_raster)) stop("Failed to build site raster. Check if files exist or set download = TRUE.")
  
  NEON_site_plots_projected <- terra::project(NEON_site_plot_polygons, terra::crs(site_raster))
  #for each shapefile (other_locations)    
  other_plots_projected <- list() #make an empty list with default of NULL to avoid errors later
  if (length(shapefile_plots) > 0) {
    other_plots_projected <- lapply(shapefile_plots, function(x) terra::project(x, terra::crs(site_raster)))
  }
  
  #combine them all
  plot_inputs <- list(NEON_site_plots_projected)
  if (length(other_plots_projected) > 0) {
    plot_inputs <- append(plot_inputs, other_plots_projected)
  }
    
  site_subplots <- get_site_subplots(which_site, which_size)
  #need site_subplots on same projection as cropped rgb plot
  site_subplots_projected <- terra::project(site_subplots, terra::crs(site_raster))
  
  if (full_site == TRUE) {
    #first, look for the pdf file name: FullSite_Map_ ... _rgb.pdf
    pdf_filename <- ifelse(which_bout == 1,
                           paste0(repo_dir, "FullSite_Map_", which_site, "_", which_year, "_rgb.pdf"),
                           paste0(repo_dir, "FullSite_Map_", which_site, "_", which_year, "_rgb_", which_bout, ".pdf"))
    if (file.exists(pdf_filename)) {
      print(paste("High-Res PDF already exists:", pdf_filename))
      return(NULL)
    } # have the pdf
    
    # If it doesn't exist, make the raster, draw the shapefiles
    print("Drawing High-Res PDF...")
    pdf(pdf_filename, width = 15, height = 15)
    terra::plotRGB(site_raster, r = 1, g = 2, b = 3, stretch = 'lin',
                   maxcell = 1e7,  #saves space
                   main = paste(which_site, which_year, "Full Site RGB"))
    terra::plot(NEON_site_plots_projected, add = TRUE, border = "black", lwd = 6)
    terra::plot(NEON_site_plots_projected, add = TRUE, border = "yellow", lwd = 3)

    if (length(other_plots_projected) > 0) { #if we have other plots, add them
      
      for (op in other_plots_projected) {
        gtype <- as.character(terra::geomtype(op))
        if (str_detect(tolower(gtype), "point")) {
          terra::plot(poly_from_pt, add = TRUE, border = "black", lwd = 6)
          terra::plot(poly_from_pt, add = TRUE, border = "turquoise", lwd = 3)
          } else {
            #if it's a polygon just map it
          terra::plot(op, add = TRUE, border = "black", lwd = 6)
          terra::plot(op, add = TRUE, border = "pink", lwd = 3)
        }
      }
    }
    dev.off()
    print(paste("Successfully saved High-Res PDF to:", pdf_filename))
    return(NULL)
  } #if full site = true
  
  #if not full site,
  site_output_list <- list() #establish site_output_list for individual plots:
  for (i in 1:length(plot_inputs)) {
    this_input <- plot_inputs[[i]] 
    if (nrow(this_input) == 0) next 
    
    #this input is a particular plot. We loop over every PLOT, NEON's and other plots if added.
    for (ID in 1:nrow(this_input)) {
      this_plot <- this_input[ID,]
      plotID <- as.character(this_plot$plotID)
      
      plot_output_list <- list()
      cropped_rgb_plot <- tryCatch({ 
        terra::crop(site_raster, terra::ext(this_plot)) }, error = function(e) NULL)
    #if we successfully cropped the plot, prep its names, save crop, and put it in the list for the plot and the site.
    if(!is.null(cropped_rgb_plot)) {
      #give cropped plot appropriate names
      names(cropped_rgb_plot) <- rgb_names
      plot_output_list[['plot_full_rgb']] <- cropped_rgb_plot
      
      #subplots in plot: 
        # Check if this plot has 1m2 subplots
        subplots_in_plot <- subset(site_subplots_projected, site_subplots_projected$plotID == plotID)
        
        if (nrow(subplots_in_plot) > 0) {
          # Swapped purrr::map_dfr for a standard loop, which is much better for lists
          for (sub in 1:nrow(subplots_in_plot)) {
            this_sub <- subplots_in_plot[sub, ]
            subplotID <- as.character(this_sub$subplotID)         
          # establish extent, then turn it into a vector / polygon
          poly_from_pt_ext <- terra::ext(terra::buffer(this_sub, which_buff))
          # Crop the 1m2 subplot directly from the 20x20m raster
          sub_crop <- tryCatch({ terra::crop(cropped_rgb_plot, poly_from_pt_ext) }, error = function(e) NULL)
          
          if(!is.null(sub_crop)) {
            #give cropped plot appropriate names
            names(sub_crop) <- rgb_names
            plot_output_list[[subplotID]] <- sub_crop
          
            } #if we have a cropped plot, add it to the plot list
          }
        }
        site_output_list[[plotID]] <- plot_output_list
        
    }
      }
    }
  return(site_output_list)
}


CPER_2020 <- make_maps('CPER', 2020, other_locations = other_locations, which_size = 1, which_buff = 2)
CPER_2021 <- make_maps('CPER', 2021, other_locations = other_locations, which_size = 1, which_buff = 2)
CPER_2024 <- make_maps('CPER', 2024, other_locations = other_locations, which_size = 1, which_buff = 2)

plotRGB(CPER_2024$CPER_001$plot_full_rgb, stretch = 'lin')
plotRGB(CPER_2024$CPER_002$plot_full_rgb, stretch = 'lin')
plotRGB(CPER_2024$CPER_003$plot_full_rgb, stretch = 'lin')

CPER2020extracted <- extract_from_maplist('CPER', 2020, other_locations = other_locations, which_size = 1, which_buff = 2)
CPER2021extracted <- extract_from_maplist('CPER', 2021, other_locations = other_locations, which_size = 1, which_buff = 2)
CPER2024extracted <- extract_from_maplist('CPER', 2024, other_locations = other_locations, which_size = 1, which_buff = 2)


#extract from maplist will get us a dataframe of the hsi reflectance from maplist 
extract_from_maplist <- function(which_site, 
                                 which_year, #no default
                                 which_bout = 1, #default
                                 other_locations = NULL, #default
                                 download = FALSE, #default
                                 which_size,
                                 which_buff,
                                 repo_dir = "/Users/khuelsma/Desktop/ARIDNEON/",
                                 output_dir = "/Users/khuelsma/") { #default

  if (is.null(which_year)) stop('specify a year to extract spectra')
  
  csv_file <- paste0(output_dir, "extracted_spectra_", which_site, "_", which_year, ".csv")
  
  if (file.exists(csv_file)) {
    print(paste("Spectra already extracted! Loading:", csv_file))
    return(read.csv(csv_file))
  }
  
  print("Extracted spectra not found. Generating maps and extracting now...")
  
  # 1. Make sure all tiny HSI plots are saved to the hard drive (skips instantly if they are)
  cache_hsi_plots(which_site, which_year, which_bout, other_locations, download, 
                  home_dir = output_dir, repo_dir = repo_dir)
  #need to save wavelength names
  wv_df <- get_wvs(which_site, which_year) #site doesn't really matter, but year does
  
  save_dir <- paste0(repo_dir, 'plots/')
  plot_files <- list.files(save_dir, pattern = paste0(which_site, "_", which_year, ".*_hsi\\.tif$"), full.names = TRUE)
  
  if (length(plot_files) == 0) {
    warning("No HSI plots were found in the cache folder!")
    return(NULL)
  }
  site_subplots <- get_site_subplots(which_site, which_size)
  
  # map_dfr loops through the list and binds the resulting dataframes together instantly
  extracted_spectra <- purrr::map_dfr(plot_files, function(file_path) { #like for each
    
    plot_rast <- terra::rast(file_path)
    names(plot_rast) <- unique(wv_df$band)
    #Extract plotID from the filename
    eventID <- str_extract(file_path, "(?<=plots//).*(?=\\_hsi.tif)")
    plot_year_string <- paste0(which_site, "_", which_year)
    plotID <- sub(paste0(plot_year_string, "_"), "", eventID)
    # Check if this plot has 1m2 subplots
    subplots_in_plot <- subset(site_subplots, site_subplots$plotID == plotID)
    
    if (nrow(subplots_in_plot) == 0) {
      # No subplots (e.g., other_DA_ plots) -> Extract full plot
      df <- terra::as.data.frame(plot_rast, xy = TRUE, cells = TRUE, na.rm = TRUE)
      if (nrow(df) == 0) return(NULL)
      
      df$year <- which_year
      df$plotID <- plotID
      df$subplotID <- "plot_full"
      df$eventID <- paste0(plotID, '_plot_full_', which_year)
      df$buff <- which_buff
      return(df)
      
      } else {
      sub_polys <- create_square_polygons(subplots_in_plot)
      buffed_subplots <- terra::buffer(sub_polys, which_buff) # which_buff = 1
      
      purrr::map_dfr(1:nrow(buffed_subplots), function(sub) {
        this_sub <- buffed_subplots[sub, ]
        subplotID <- as.character(subplots_in_plot$subplotID[sub])
        
        # Crop the 1m2 subplot directly from the 20x20m raster
        sub_crop <- tryCatch({ terra::crop(plot_rast, terra::ext(this_sub)) }, error = function(e) NULL)
        if (is.null(sub_crop)) return(NULL)
        
        df <- terra::as.data.frame(sub_crop, xy = TRUE, cells = TRUE, na.rm = TRUE)
        if (nrow(df) == 0) return(NULL)
          df$year <- which_year
          df$plotID <- plotID
          df$subplotID <- subplotID
          df$eventID <- paste0(plotID, '_', subplotID, '_', which_year, '_', bout_num)
          df$buff <- which_buff
          
          return(df)
        })
    }
  }) 
  write.csv(extracted_spectra, csv_file, row.names = FALSE)
  print(paste("Successfully saved extracted spectra to:", csv_file))
  return(extracted_spectra)
}


# End of functions --------------------------------------------------------


# Beginning of analysis ---------------------------------------------------

both_avail('CPER')

#if you need to download: use avail_file_list OR get_site_raster, which uses avail_file_list as input
# OR make_maps, which uses get_site_raster


# want rgb and vi maps of full site and of plots ('rgb' and 'vi' for which_rast; full_site = TRUE and FALSE)
# Define the years you want
my_years <- c(2021, 2024)

# map_dfr will run the function for each year and stitch the dataframes together
site_plots_hsi_df <- map_dfr(my_years, function(yr) {
  extract_from_maplist(which_site = 'CPER', 
                       which_year = yr, 
                       which_bout,
                       other_locations = other_locations,
                       which_size = 1,
                       which_buff = which_buff)})
                        
#rgb maps of full site:
make_maps(which_rast = 'rgb', 
          which_site = 'CPER', 
          which_year = 2024,
          other_locations = other_locations, 
          full_site = TRUE, 
          download = TRUE, which_size = 1, which_buff = 2)
#vi maps of full site (just ran for 2021)

make_maps(which_rast = 'vi', 
          which_site = 'CPER', 
          which_year,
          other_locations = other_locations, 
          full_site = TRUE, 
          download = TRUE)

# rgb maps of plots:
rgb_maps_2021 <- make_maps(which_rast = 'rgb', 
                           which_site = 'CPER', 
                           which_year = 2021, 
                           other_locations = other_locations, 
                           download = TRUE, 
                           full_site = TRUE, 
                           which_size = 1, 
                           which_buff = 2)

plotRGB(rgb_maps_2021$other$plot_full, stretch = 'lin')
plotRGB(rgb_maps_2021$other_DA_MULTIPOLYGON_42$plot_full, stretch = 'lin')

rgb_maps_2024 <- make_maps(which_rast = 'rgb', 
                           which_site = 'CPER', 
                           which_year = 2024, 
                           other_locations = NULL, 
                           download = FALSE, 
                           full_site = FALSE, 
                           which_size = 1, 
                           which_buff = 1)

# want hsi maps created for the sake of extraction but NOT full site. 
site_plots_hsi_df <- extract_from_maplist(which_site = 'CPER', 
                                          which_year = c(2021, 2024), 
                                          which_bout = 1,
                                          other_locations = FALSE,
                                          download = TRUE)

spec_sampled <- site_plots_hsi_df %>%
  dplyr::select(plotID, subplotID, eventID,
                cell, x, y,
                starts_with("B")) %>%
  tidyr::pivot_longer(cols = starts_with("B"), 
                      names_to = "band", 
                      values_to = "refl") %>%
  dplyr::left_join(wvs, by = c("band")) %>%
  dplyr::filter(!is.na(wv)) %>% 
  print()

#how many observations came from each
n_obs <- spec_sampled %>%
  dplyr::group_by(plotID, subplotID, eventID) %>%
  summarize(cells = length(unique(cell))) %>%
  print(n = 600)

#a summary of each -- reflectance, variability, with how many observations
spec_summary <- spec_sampled %>%
  dplyr::group_by(plotID, subplotID, eventID, wv) %>%
  summarize(m_refl = mean(refl),
            sd_refl = sd(refl),
            se_refl = sd_refl/sqrt(n()),
            CV_refl = sd_refl/m_refl) %>%
  left_join(n_obs)

fig_all_refls <- spec_sampled %>%
  ggplot(aes(x = wv, y = refl)) +
  geom_point(aes(colour = plotID), size= 0.5) +
  xlab('wavelength(nm)') + ylab('all reflectances from surveyed space') +
  facet_wrap(~plotID) +
  theme_classic()

fig_mean_sub_refl <- spec_summary %>%
  ggplot(aes(x = wv, y = m_refl)) +
  geom_point(aes(colour = subplotID), size= 0.5) +
  geom_ribbon(aes(x = wv, 
                    ymax = m_refl + se_refl, 
                    ymin = m_refl - se_refl,
                    fill = subplotID),
              alpha = 0.5) +
  xlab('wavelength(nm)') + ylab('mean subplot reflectance') +
  facet_wrap(~plotID) +
  theme_classic()

fig_CV_sub_refl <- spec_summary %>%
  ggplot(aes(x = wv, y = CV_refl)) +
  geom_point(aes(colour = subplotID), size= 0.5) +
  xlab('wavelength(nm)') + ylab('subplot variability in reflectance') +
  facet_wrap(~plotID) +
  theme_classic()

fig_all_refls
fig_mean_sub_refl
fig_CV_sub_refl
#now we can take the sampled spectra and pair them w veg data :) 