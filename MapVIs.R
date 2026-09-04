# Overview: ---------------------------------------------------------------

# MapVIs lets you

# Load saved tiffs of RGB VIs, OR make the map from scratch
# open and extract necessary info to combine w vegetation


#we can load the VI RGB maps and clip them to subplot / plot sizes
foreach(t = 1:23) %do% { #length(CPER_hsi_tile_list)) %do% {
  #load the saved raster and extract polygons:
  filename <- paste0('/Users/khuelsma/indices_', which_site, '_', which_year, '_', t, '.tif')
  VI_rast <- terra::rast(filename)
  tile_extent <- terra::ext(VI_rast)
  plots_projected <- terra::project(subplot_shp_1, VI_rast) #can also use David's data here
  plots_in_tile <- plots_projected[tile_extent]
  
  #full tile maps of VIs
  terra::plotRGB(VI_rast, r = 1, g = 2, b = 3, stretch = 'lin')
  terra::plot(plots_in_tile, add = TRUE, lwd = 3, col = 'white')
  
  plots_df <- as.data.frame(plots_in_tile)
  
  # (Added a check in case no plots fall in this tile)
  if (nrow(plots_in_tile) > 0) {
    #create a translation between extraction ID and plotID
    plot_lookup <- data.frame(
      ID = 1:nrow(plots_in_tile),
      plotID = plots_in_tile$plotID 
    )
    
    foreach(p = 1:nrow(plots_in_tile)) %do% {
      this_plot <- plots_in_tile[p,]
      buffer_pt <- terra::buffer(this_plot, width = 2)
      
      cropped_plot <- terra::crop(VI_rast, buffer_pt)
      #mapping
      terra::plotRGB(cropped_plot, r = 1, g = 2, b = 3, stretch = 'lin')
      
      #extracting:
      subplot_metrics <- terra::extract(VI_rast, buffer_pt)
      subplot_metrics <- subplot_metrics %>%
        rename(plotnum = ID) %>%
        mutate(tile = filename) %>%
        left_join(plots_df)
      
    }
  }
}

# open and extract necessary info, interpret the data

