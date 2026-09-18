# Overview: ---------------------------------------------------------------

# MapVIs lets you map rgb or vegetation indices for all tiles at a site for a given year

#using the tile list, we can make a full site map or extract plots and subplots.
#so if I say make_maps, I want to load every tile, but extract the plot [[plot_ID]] and then subplots [[subplot_ID]]
#for make_maps, save a list of them. For extract_subplot_spectra, 
make_maps <- function(which_rast, 
                      which_site, 
                      which_year, 
                      full_site = FALSE,
                      which_buff) {
  
  site_raster <- get(paste0(which_site, '_', which_year, '_', which_rast))
  # only if file doesn't exist: site_raster <- create_site_raster(which_rast, which_site, which_year)
  
  # Import plots shape files
  plots_shp <- terra::vect('/Users/khuelsma/Desktop/NEON Spectral Variability/Relevant NEON Materials/NEON TOS Plots/All_NEON_TOS_Plot_Polygons_V11.shp')
  site_plot_polygons <- subset(plots_shp, 
                               stringr::str_detect(plots_shp$plotID, which_site ) & 
                                 stringr::str_detect(plots_shp$appMods, 'div'))
  # Project plots to match the raster CRS
  site_plots_projected <- terra::project(site_plot_polygons, terra::crs(site_raster)) #can also use David's data here
  
  if (full_site == TRUE) {
    terra::plotRGB(site_raster, r = 1, g = 2, b = 3, stretch = 'lin')
    #terra::plot(site_plots_projected, add = TRUE, lwd = 3, border = 'white')
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
      
      # Initialize a nested list for THIS specific plot
      plot_list <- list()
      
      # Save the full 20x20m plot raster
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

plotmaps_CPER_2021_rgb <- make_maps('rgb', 'CPER', 2021, full_site = TRUE, which_buff = 0)

plotmaps_CPER_2021_vi <- make_maps('vi', 'CPER', 2021, full_site = FALSE, which_buff)
plotmaps_CPER_2024_vi <- make_maps('vi', 'CPER', 2024, full_site = FALSE, which_buff)
plotmaps_CPER_2021_rgb <- make_maps('rgb', 'CPER', 2021, full_site = FALSE, which_buff)
plotmaps_CPER_2024_rgb <- make_maps('rgb', 'CPER', 2024, full_site = FALSE, which_buff)