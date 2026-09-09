# Overview: ---------------------------------------------------------------

# MapVIs lets you
# Load saved tiffs of RGB VIs
# extract necessary info to combine w vegetation

#input_with_paths <- map_site_tiles('CPER', 2024)
input_with_paths <- read.csv('CPER_2024_tiles.csv')

hsi_list <- unique(input_with_paths$hsi_saved_path)
rgb_list <- unique(input_with_paths$rgb_filepath)
spbands_list <- unique(input_with_paths$spbands_saved_path)

#what we will extract:
subplots_shp <- terra::vect(
  '/Users/khuelsma/Desktop/NEON Spectral Variability/Relevant NEON Materials/NEON TOS Plots/All_NEON_TOS_Plot_Subplots_V11.shp')
#subset the shape using only plots at the site and "div" (diversity) applications
site_subplots_shp <- subset(subplots_shp, 
                            stringr::str_detect(subplots_shp$plotID, which_site ) & 
                              stringr::str_detect(subplots_shp$appMods, 'div'))

input <- spbands_list
VI_list <- c() #initialize empty list
#using spbands_list, make vegetation index rasters and an RGB composite: 
foreach(t = 1:length(input)) %do% {
  this_tile <- terra::rast(input[t])
  tile_extent <- ext(this_tile)
  
    #make rasters
    PRI_rast <- (this_tile[['NIR']] - this_tile[['green']]) / (this_tile[['NIR']] + this_tile[['green']])
    PRI_col = colorRampPalette(c("white", "blue"))(20)
    NDVI_rast <- (this_tile[['NIR']] - this_tile[['red']]) / (this_tile[['NIR']] + this_tile[['red']] + 0.0001)
    
    NIRv_rast <- this_tile[['NIR']] * NDVI_rast
    NIRv_col = colorRampPalette(c("white", "darkgreen"))(20)
    CCI_rast <- (this_tile[['PRI']] - this_tile[['red']]) / (this_tile[['PRI']] + this_tile[['red']])
    CCI_col = colorRampPalette(c("white", "red"))(20)
    
    indices <- c(CCI_rast, NIRv_rast, PRI_rast)
    names(indices) <- c("CCI", "NIRv", "PRI")
    
    terra::plot(CCI_rast, col = CCI_col, main = "CCI")
    terra::plot(NIRv_rast, col = NIRv_col, main = "NIRv")
    terra::plot(PRI_rast, col = PRI_col, main = 'PRI')
    VI_rast <- terra::plotRGB(indices, r=1, g=2, b=3, stretch="lin")
    
    #arrange these:
    
    #save the rasters as tiffs or a list
    #terra::writeRaster(indices, filename = paste0('indices_', which_site, '_', which_year, '_', t, '.tif'))
    VI_list[[paste0(t)]] <- indices  
    return(VI_list)
}

#for every tile file input:
input <- rgb_list
tiles_only <- TRUE
site_tilelist <- c()
#mapping rgb
foreach(t = 1:length(input),
                 .combine = rbind) %do% {
                   
  filename <- input[t]
  this_tile <- terra::rast(input[[t]])
  tile_extent <- ext(this_tile)
  
  #default map:
  tilemap <- terra::plotRGB(this_tile) #r = 1, g = 2, b = 3 is default
  #dramatic map:
  tilemap <- terra::plotRGB(this_tile, r=1, g=2, b=3, stretch="lin")

  if (tiles_only == TRUE) {
    #saves a list of them for full site mapping
    site_tilelist[[t]] <- this_tile
  }
  if (tiles_only == FALSE) {
    
    plots_projected <- terra::project(site_subplots_shp, this_tile) #can also use David's data here
    plots_in_tile <- plots_projected[tile_extent]
    
    if (nrow(plots_in_tile) > 0) { #if there are plots in a tile, create lookup table
      #create a translation between extraction ID and plotID
      
      hello <- foreach(plot = nrow(plots_in_tile)) %do% {
                
                subplot <- plots_in_tile[plot,] #for each plot, numbered 1 to nrow
                
                subplot_df <- as.data.frame(subplot)
                subplotID = subplot_df$subplotID
                
                buffer_pt <- terra::buffer(subplot, width = 2) #buffer so we account for spatial uncertainty?
                cropped_plot <- terra::crop(this_tile, buffer_pt) #crop the map to just the buffered point

                terra::plotRGB(cropped_plot)
                terra::plotRGB(cropped_plot, r = 1, g = 2, b = 3, stretch = 'lin') #plot the cropped plot
                
                subplot_df  #for rbinding
              }
    } #plots in tile
    
  } #each tile from the list
}

#create full site rgb image:
full_site_rgb10 <- terra::merge(terra::sprc(site_tilelist))

full_site_dullrgb <- terra::plotRGB(full_site_rgb10)
full_site_stretchrgb <- terra::plotRGB(full_site_rgb10, stretch = 'lin')


# LOCATION INFO -----------------------------------------------------------
#import shapefiles or lat/lon csv
all_NEON_plots <- read.csv(file = '/Users/khuelsma/ARIDNEON/All_NEON_TOS_Plot_Centroids_V11.csv')
#NEON shape files: choose one to assign to LL
plots_shp <- terra::vect('/Users/khuelsma/Desktop/NEON Spectral Variability/Relevant NEON Materials/NEON TOS Plots/All_NEON_TOS_Plot_Polygons_V11.shp')
#project when needed later
subplots_shp <- terra::vect('/Users/khuelsma/Desktop/NEON Spectral Variability/Relevant NEON Materials/NEON TOS Plots/All_NEON_TOS_Plot_Subplots_V11.shp')

#buff_extent <- terra::ext(terra::buffer(transects, width = 5))

if (tiles_only == FALSE) {
  
  plots_projected <- terra::project(site_subplots_shp, this_tile) #can also use David's data here
  plots_in_tile <- plots_projected[tile_extent]
  
  if (nrow(plots_in_tile) > 0) { #if there are plots in a tile, create lookup table
    #create a translation between extraction ID and plotID
    
    plot_lookup <- data.frame(
      ID = 1:nrow(plots_in_tile), #we will make a dataframe with plot ID (just a number)
      subplotID = plots_in_tile$subplotID, #the subplot and plot ID (from the dataframe)
      plotID = plots_in_tile$plotID #so that we can reference which plot is which
    )
    
    hello <- foreach(plot = nrow(plots_in_tile)) %do% {
      
      subplot <- plots_in_tile[plot,] #for each plot, numbered 1 to nrow
      
      subplot_df <- as.data.frame(subplot)
      subplotID = subplot_df$subplotID
      
      buffer_pt <- terra::buffer(subplot, width = 2) #buffer so we account for spatial uncertainty?
      cropped_plot <- terra::crop(this_tile, buffer_pt) #crop the map to just the buffered point
      
      terra::plotRGB(cropped_plot)
      terra::plotRGB(cropped_plot, r = 1, g = 2, b = 3, stretch = 'lin') #plot the cropped plot
      
      subplot_metrics  #for rbinding
    }
  } #plots in tile
  
  
  
input <- hsi_list
#extracting:
hello <- foreach(t = 1:3,
                 .combine = rbind) %do% {
                   
                   this_tile <- terra::rast(input[[t]])
                   tile_extent <- ext(this_tile)
                   
                   plots_projected <- terra::project(site_subplots_shp, this_tile) #can also use David's data here
                   plots_in_tile <- plots_projected[tile_extent]
                     
                     
                   if (nrow(plots_in_tile) > 0) { #if there are plots in a tile, create lookup table
                       #create a translation between extraction ID and plotID
                       
                     plot_lookup <- data.frame(
                       ID = 1:nrow(plots_in_tile), #we will make a dataframe with plot ID (just a number)

                       subplotID = plots_in_tile$subplotID, #the subplot and plot ID (from the dataframe)
                       plotID = plots_in_tile$plotID #so that we can reference which plot is which
                       )
                       
                       hello <- foreach(
                         plot = nrow(plots_in_tile),
                                        
                         .combine = rbind) %do% {
                           subplot <- plots_in_tile[plot,] #for each plot, numbered 1 to nrow
                           subplot_df <- as.data.frame(subplot)
                           subplotID = subplot_df$subplotID
                           buff_size = 0.5*sqrt(subplot_df$subpltSize) +2
                           subplot_metrics <- terra::extract(this_tile, 
                                                             subplot, 
                                                             buffer = buff_size,
                                                             fun = mean,
                                                             xy = TRUE,
                                                             cells = TRUE)
                         
                           print(subplot_metrics)
                           }
                   }
                 }
}
                           
#each tile from the list

hello %>%
  group_by(plotID, subplotID) %>%
  summarise(hi = n_distinct(cell))

#full tile maps of VIs
terra::plotRGB(VI_rast, r = 1, g = 2, b = 3, stretch = 'lin')
terra::plot(plots_in_tile, add = TRUE, lwd = 3, col = 'white')

    