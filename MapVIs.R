# Overview: ---------------------------------------------------------------

# MapVIs lets you
# Load saved tiffs of RGB VIs
# extract necessary info to combine w vegetation

input_with_paths <- map_site_tiles('CPER', 2024)

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


# Visualize reflectances and save a simple (5-band) preliminary map, which can then be used to calculate VIs

input <- CPER_2024$spbands_saved_path

#from the tile list... grab relevant locations:

foreach(t = 1:length(input)) %do% {
  this_tile <- terra::rast(input[t])
  tile_extent <- ext(this_tile)
  
  plots_projected <- terra::project(plots_shp, this_tile) #can also use David's data here
  plots_in_tile <- plots_projected[tile_extent]
  
  print(plots_in_tile)
}
  #make rasters
  PRI_rast <- (this_tile[['NIR']] - this_tile[['green']]) / (this_tile[['NIR']] + this_tile[['green']])
  
  NDVI_rast <- (this_tile[['NIR']] - this_tile[['red']]) / (this_tile[['NIR']] + this_tile[['red']] + 0.0001)
  
  NIRv_rast <- this_tile[['NIR']] * NDVI_rast
  
  CCI_rast <- (this_tile[['PRI']] - this_tile[['red']]) / (this_tile[['PRI']] + this_tile[['red']])
  
  indices <- c(CCI_rast, NIRv_rast, PRI_rast)
  names(indices) <- c("CCI", "NIRv", "PRI")
  
  #VI_rast <- terra::plotRGB(indices, r=1, g=2, b=3, stretch="lin")
  
  if (nrow(plots_in_tile) > 0) { #if there are plots in a tile, create lookup table
    #create a translation between extraction ID and plotID
    plot_lookup <- data.frame(
      ID = 1:nrow(plots_in_tile),
      plotID = plots_in_tile$plotID 
    ) %>%
      distinct()
  }
}
    
    
    foreach(p = 1:nrow(plot_lookup)) %do% {
      this_plot <- plots_in_tile[p,] #each row
      buffer_pt <- terra::buffer(this_plot, width = 2)
      cropped_plot <- terra::crop(VI_rast, buffer_pt)
      
      terra::plotRGB(cropped_plot, r = 1, g = 2, b = 3, stretch = 'lin')
      
      
    }
  }
}
 
      #extracting:
      subplot_metrics <- terra::extract(VI_rast, buffer_pt)
      subplot_metrics <- subplot_metrics %>%
        rename(plotnum = ID) %>%
        mutate(tile = filename) %>%
        left_join(plots_df) %>%
        print()
      
      subplot_metrics


#full tile maps of VIs
terra::plotRGB(VI_rast, r = 1, g = 2, b = 3, stretch = 'lin')
terra::plot(plots_in_tile, add = TRUE, lwd = 3, col = 'white')

plots_df <- as.data.frame(plots_in_tile)




  plots_projected <- terra::project(subplot_shp_1, VI_rast) #can also use David's data here

    plots_df <- as.data.frame(plots_in_tile)
  

    