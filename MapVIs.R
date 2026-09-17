# Overview: ---------------------------------------------------------------

# MapVIs lets you map rgb or vegetation indices for all tiles at a site for a given year

which_site <- 'CPER'
which_year <- 2024
which_rast <- 'vi'
map_tile_spectra(which_rast, which_site, which_year)

map_tile_spectra <- function(which_rast, #vi or rgb?
                                which_site, 
                                which_year) {
  #bring in the _tiles.csv file paths for downloaded tiles
  input_with_paths <- read.csv(paste0(which_site, '_', which_year, '_tiles.csv'))
  
  plots_shp <- terra::vect('/Users/khuelsma/Desktop/NEON Spectral Variability/Relevant NEON Materials/NEON TOS Plots/All_NEON_TOS_Plot_Polygons_V11.shp')
  site_plot_polygons <- subset(plots_shp, 
                               stringr::str_detect(plots_shp$plotID, which_site ) & 
                                 stringr::str_detect(plots_shp$appMods, 'div'))
  subplots_shp <- terra::vect(
    '/Users/khuelsma/Desktop/NEON Spectral Variability/Relevant NEON Materials/NEON TOS Plots/All_NEON_TOS_Plot_Subplots_V11.shp')
  site_subplot_pts <- subset(subplots_shp, 
                             stringr::str_detect(subplots_shp$plotID, which_site) & 
                               stringr::str_detect(subplots_shp$appMods, 'div'))
  if (which_rast == 'vi') {
    
    site_tilelist <- c() #initialize empty list
    
    spbands_list <- input_with_paths %>%
      distinct(spbands_saved_path)
    spbands_list <- unique(spbands_list$spbands_saved_path) %>% na.omit()

    foreach(t = 1:length(spbands_list)) %do% {
      this_tile <- terra::rast(spbands_list[[t]])
      
      #make rasters
      PRI_rast <- (this_tile[['NIR']] - this_tile[['green']]) / (this_tile[['NIR']] + this_tile[['green']])
      NDVI_rast <- (this_tile[['NIR']] - this_tile[['red']]) / (this_tile[['NIR']] + this_tile[['red']] + 0.0001)
      
      NIRv_rast <- this_tile[['NIR']] * NDVI_rast
      CCI_rast <- (this_tile[['PRI']] - this_tile[['red']]) / (this_tile[['PRI']] + this_tile[['red']])
      
      indices_raster <- c(CCI_rast, NIRv_rast, PRI_rast)
      names(indices_raster) <- c("CCI", "NIRv", "PRI")

      #saves a list of them for full site mapping
      site_tilelist[[t]] <- indices_raster
      
      NIRv_col = colorRampPalette(c("white", "darkgreen"))(20)
      PRI_col = colorRampPalette(c("white", "blue"))(20)
      CCI_col = colorRampPalette(c("white", "red"))(20)
      terra::plot(CCI_rast, col = CCI_col, main = "CCI")
      terra::plot(NIRv_rast, col = NIRv_col, main = "NIRv")
      terra::plot(PRI_rast, col = PRI_col, main = 'PRI')
      terra::plotRGB(indices_raster, r=1, g=2, b=3, stretch="lin")
      
      } #for each file 
    } #end vi
  
  #mapping rgb
  if (which_rast == 'rgb') {
    site_tilelist <- c() #initialize empty list
    
    rgb_list <- input_with_paths %>%
      distinct(rgb_filepath)
    rgb_tilelist <- unique(rgb_list$rgb_filepath) %>% na.omit()
    foreach(t = 1:length(rgb_tilelist)) %do% {
      this_tile <- terra::rast(rgb_tilelist[[t]])
      site_tilelist[[t]] <- this_tile
      
      terra::plotRGB(this_tile, r = 1, g = 2, b = 3, stretch = 'lin')
    }
  }
  return(site_tilelist)
}
#go through all tiles to plot PLOTS and SUBPLOTS

which_site <- 'CPER'
which_year <- 2021
CPER_2021_vi <- map_tile_spectra('vi', which_site, which_year)
CPER_2021_rgb <- map_tile_spectra('rgb', which_site, which_year)

which_year <- 2024
CPER_2024_vi <- map_tile_spectra('vi', which_site, which_year)
CPER_2024_rgb <- map_tile_spectra('rgb', which_site, which_year)


#using the tile list, we can make a full site map or extract plots:

#I want to edit this to create maps of every plot and subplots within
#so if I say make_maps, I want to load every tile, but extract the plot [[plot_ID]] and then subplots [[subplot_ID]]
make_maps <- function(which_rast, 
                      which_site, 
                      which_year, 
                      full_site = FALSE) {
  
  # Get the list of raster tiles
  tile_list <- map_tile_spectra(which_rast, which_site, which_year)
  merged_tiles <- terra::merge(terra::sprc(tile_list))
  
  # Import plots shape files
  plots_shp <- terra::vect('/Users/khuelsma/Desktop/NEON Spectral Variability/Relevant NEON Materials/NEON TOS Plots/All_NEON_TOS_Plot_Polygons_V11.shp')
  site_plot_polygons <- subset(plots_shp, 
                               stringr::str_detect(plots_shp$plotID, which_site ) & 
                                 stringr::str_detect(plots_shp$appMods, 'div'))
  # Project plots to match the raster CRS
  plots_projected <- terra::project(site_plot_polygons, merged_tiles) #can also use David's data here
  
  site_subplots <- get_site_subplots(which_site, which_size)
  subplots_projected <- terra::project(site_subplots, merged_tiles)
  
  if (full_site == TRUE) {
    terra::plotRGB(merged_tiles, r = 1, g = 2, b = 3, stretch = 'lin')
    terra::plot(plots_projected, add = TRUE, lwd = 3, border = 'white')
    return(merged_tiles)
  }
  
  # 5. INDIVIDUAL PLOT & SUBPLOT EXTRACTION
  site_output_list <- list() # Initialize the primary named list
  
  # Loop over every PLOT in the site
  for (ID in 1:nrow(plots_projected)) {
    this_plot <- plots_projected[ID, ]
    plotID <- as.character(this_plot$plotID)
    plot_ext <- terra::ext(this_plot)
    
    # Try cropping (handles cases where a plot might be outside the downloaded tiles)
    cropped_plot <- tryCatch({
      terra::crop(merged_tiles, plot_ext)
    }, error = function(e) NULL)
    
    # If the plot exists inside our raster boundary, proceed:
    if (!is.null(cropped_plot)) {
      
      # Initialize a nested list for THIS specific plot
      plot_list <- list()
      
      # Save the full 20x20m plot raster
      plot_list[['plot_full']] <- cropped_plot
      
      # Find all subplots that belong to this plot
      subplots_in_plot <- subset(subplots_projected, subplots_projected$plotID == plotID)
      
      # If subplots exist, extract them
      if (nrow(subplots_in_plot) > 0) {
        for (sub in 1:nrow(subplots_in_plot)) {
          this_sub <- subplots_in_plot[sub, ]
          subplotID <- as.character(this_sub$subplotID)
          
          # Compute buffer for subplot mapping
          size <- this_sub$subpltSize
          if (is.na(size) || is.null(size)) size <- 400 # Default to 400m2 if missing
          
          subplot_size <- 0.5 * sqrt(as.numeric(size))
          buff_size <- 2 # Expand bounding box by 2m
          
          # Buffer the subplot point and grab its extent
          buff_extent <- terra::ext(terra::buffer(this_sub, width = subplot_size + buff_size))
          
          # Crop from the ALREADY cropped plot for speed
          cropped_subplot <- tryCatch({
            terra::crop(cropped_plot, buff_extent)
          }, error = function(e) NULL)
          
          # Save to the nested list using the subplotID
          if (!is.null(cropped_subplot)) {
            plot_list[[subplotID]] <- cropped_subplot
          }
        }
      }
      
      # Add this plot's completed nested list to the main site list
      site_output_list[[plotID]] <- plot_list
    }
  }
  
  return(site_output_list)
}


plotmaps_CPER_2021_vi <- make_maps('vi', 'CPER', 2021, full_site = FALSE)
plotmaps_CPER_2024_vi <- make_maps('vi', 'CPER', 2024, full_site = FALSE)

plotmaps_CPER_2021_rgb <- make_maps('rgb', 'CPER', 2021, full_site = FALSE)
plotmaps_CPER_2024_rgb <- make_maps('rgb', 'CPER', 2024, full_site = FALSE)


map_plot_subplot <- function(which_site, which_year, which_plot, which_subplot) {
  
  rgb_site_year_list <- get(paste0('plotmaps_', which_site, '_', which_year, '_rgb'))
  vi_site_year_list <- get(paste0('plotmaps_', which_site, '_', which_year, '_vi'))

  subplotrgb <- terra::plotRGB(rgb_site_year_list[[which_plot]][[which_subplot]],
                               r = 1, g = 2, b = 3, stretch = 'lin')
  subplotvi <- terra::plotRGB(vi_site_year_list[[which_plot]][[which_subplot]],
                              r = 1, g = 2, b = 3, stretch = 'lin')
}
  
map_plot_subplot('CPER', 2021, 'CPER_001', 'plot_full')

