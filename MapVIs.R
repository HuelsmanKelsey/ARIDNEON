# Overview: ---------------------------------------------------------------

# MapVIs lets you

# Load saved tiffs of RGB VIs, OR make the map from scratch
# extract necessary info to combine w vegetation


#from the tile list... grab relevant locations:
foreach(t = 1:23) %do% {
  #take the already created raster file: this_tile <- paste0('Users/khuelsma/indices_CPER_2024_', t)
  # or make it now:
  this_tile <- CPER_hsi_tile_list[[img]]
  
  tile_extent <- ext(this_tile)
  plots_projected <- terra::project(plots_shp, this_tile) #can also use David's data here
  plots_in_tile <- plots_projected[tile_extent]
  
  #make rasters
  PRI_rast <- (this_tile[['NIR']] - this_tile[['green']]) / (this_tile[['NIR']] + this_tile[['green']])
  
  NDVI_rast <- (this_tile[['NIR']] - this_tile[['red']]) / (this_tile[['NIR']] + this_tile[['red']] + 0.0001)
  
  NIRv_rast <- this_tile[['NIR']] * NDVI_rast
  
  CCI_rast <- (this_tile[['PRI']] - this_tile[['red']]) / (this_tile[['PRI']] + this_tile[['red']])
  
  indices <- c(CCI_rast, NIRv_rast, PRI_rast)
  
  names(indices) <- c("CCI", "NIRv", "PRI")
  
  indices_rast <- terra::plotRGB(indices, r=1, g=2, b=3, stretch="lin")
  
  #export whole tile: terra::writeRaster(indices, filename = paste0('indices_', which_site, '_', which_year, '_', t, '.tif'))
}



# Load saved tiffs of RGB VIs -------------------------------------------
#we can load the VI RGB maps and clip them to subplot / plot sizes
VIs <- foreach(t = 1:23,
        .combine = rbind) %do% { # or length(CPER_hsi_tile_list)) %do% {
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
      return(subplot_metrics)
    }
  }
}

VIs
# open and extract necessary info, interpret the data

