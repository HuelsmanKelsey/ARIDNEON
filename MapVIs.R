# Overview: ---------------------------------------------------------------

# MapVIs lets you
# Load saved tiffs of VIs

#maps are already made:
input_with_paths <- read.csv('CPER_2024_tiles.csv')

hsi_list <- unique(input_with_paths$hsi_saved_path)
rgb_list <- unique(input_with_paths$rgb_filepath)
spbands_list <- unique(input_with_paths$spbands_saved_path)

#mapping rgb
input <- rgb_list %>% na.omit() #for every rgb tile file input
#tiles_only = just map and save the tiles in site_tilelist; 
# if tiles_only == FALSE, we also clip out PLOTS
tiles_only <- TRUE
site_tilelist <- c()

foreach(t = 1:length(input),
        .combine = rbind) %do% {
          
          this_tile <- terra::rast(input[[t]])
          tile_extent <- ext(this_tile)
            
          if (tiles_only == TRUE) {
            #saves a list of them for full site mapping
            site_tilelist[[t]] <- this_tile
          }
          
          #default map:
          tilemap <- terra::plotRGB(this_tile) #r = 1, g = 2, b = 3 is default
          
          #dramatic map:
          tilemap <- terra::plotRGB(this_tile, r=1, g=2, b=3, stretch="lin")
          
          
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
              
              hello <- foreach(ID = nrow(plots_in_tile),
                               .combine = rbind) %do% {
                                 
                                 this_plot <- plots_in_tile[ID,] #for each plot, numbered 1 to nrow
                                 
                                 this_plot_df <- as.data.frame(this_plot)
                                 subplotID = this_plot_df$subplotID
                                 
                                 buff_size = 0.5*sqrt(this_plot_df$subpltSize) +2
                                 buff_extent <- terra::ext(terra::buffer(this_plot, 
                                                                         width = buff_size))
                                 cropped_plot <- terra::crop(this_tile, buff_extent) #crop the map to just the buffered point
                                 
                                 terra::plotRGB(cropped_plot)
                                 terra::plotRGB(cropped_plot, r = 1, g = 2, b = 3, stretch = 'lin') #plot the cropped plot
                                 
                                 this_plot_df
                               }
            } #plots in tile
            hello
          } #tiles only = FALSE
          
        } #each file

#create full site rgb image:
full_site_rgb10 <- terra::merge(terra::sprc(site_tilelist))

full_site_dullrgb <- terra::plotRGB(full_site_rgb10)
full_site_stretchrgb <- terra::plotRGB(full_site_rgb10, stretch = 'lin')


input <- spbands_list
VI_list <- c() #initialize empty list
#using spbands_list, make vegetation index rasters and an RGB composite: 
foreach(t = 1:length(input)) %do% {
  this_tile <- terra::rast(input[t])
  tile_extent <- ext(this_tile)
  
  #make rasters
  PRI_rast <- (this_tile[['NIR']] - this_tile[['green']]) / (this_tile[['NIR']] + this_tile[['green']])
  NDVI_rast <- (this_tile[['NIR']] - this_tile[['red']]) / (this_tile[['NIR']] + this_tile[['red']] + 0.0001)
  
  NIRv_rast <- this_tile[['NIR']] * NDVI_rast
  CCI_rast <- (this_tile[['PRI']] - this_tile[['red']]) / (this_tile[['PRI']] + this_tile[['red']])
  
  indices <- c(CCI_rast, NIRv_rast, PRI_rast)
  names(indices) <- c("CCI", "NIRv", "PRI")
  
  VI_list[[t]] <- indices  
}

NIRv_col = colorRampPalette(c("white", "darkgreen"))(20)
PRI_col = colorRampPalette(c("white", "blue"))(20)
CCI_col = colorRampPalette(c("white", "red"))(20)

input <- VI_list
foreach(t = 1:length(input)) %do% {
  this_tile <- VI_list[[t]]
  terra::plot(this_tile$CCI, col = CCI_col, main = "CCI")
  terra::plot(this_tile$NIRv, col = NIRv_col, main = "NIRv")
  terra::plot(this_tile$PRI, col = PRI_col, main = 'PRI')
  
  VI_rast <- terra::plotRGB(this_tile, r=1, g=2, b=3, stretch="lin")
}

#create fullsite vi map:
full_site_vis <- terra::merge(terra::sprc(VI_list))
terra::plotRGB(full_site_vis, r = 1, g = 2, b = 3, stretch = 'lin')

#using VI_list and rgb list, crop the POLYGON plots

#1) a file list of tiffs: hsi_list or rgb_list or spbands_list
#which_rast = hsi_list, rgb_list, spbands_list

rgb_map_sampling_plot <- function(which_site) {
  input <- rgb_list %>% na.omit()
  
  plots_shp <- terra::vect('/Users/khuelsma/Desktop/NEON Spectral Variability/Relevant NEON Materials/NEON TOS Plots/All_NEON_TOS_Plot_Polygons_V11.shp')
  site_plot_polygons <- subset(plots_shp, 
                               stringr::str_detect(plots_shp$plotID, which_site ) & 
                                 stringr::str_detect(plots_shp$appMods, 'div'))
  subplots_shp <- terra::vect(
    '/Users/khuelsma/Desktop/NEON Spectral Variability/Relevant NEON Materials/NEON TOS Plots/All_NEON_TOS_Plot_Subplots_V11.shp')
  site_subplot_pts <- subset(subplots_shp, 
                             stringr::str_detect(subplots_shp$plotID, which_site) & 
                               stringr::str_detect(subplots_shp$appMods, 'div'))
  #for each tile from input
  foreach(
    img = 1:length(input)) %do% {
      this_tile <- terra::rast(input[[img]])
      tile_extent <- ext(this_tile)
      
      plots_projected <- terra::project(site_plot_polygons, this_tile) #can also use David's data here
      plots_in_tile <- plots_projected[tile_extent]
      subplots_projected <- terra::project(site_subplot_pts, this_tile)
      subplots_in_tile <- subplots_projected[tile_extent]
      
      if (nrow(plots_in_tile) > 0) {
        foreach( 
          plot = 1:nrow(plots_in_tile)) %do% {
            
            this_plot <- plots_in_tile[plot,] #for each plot, numbered 1 to nrow
            plot_ext <- ext(this_plot)
            subplots_in_plot <- subplots_in_tile[plot_ext]

            this_plot_df <- as.data.frame(this_plot)
            plotID = this_plot_df$plotID
            
            buff_extent <- plot_ext
            cropped_plot <- terra::crop(this_tile, buff_extent) #crop the map to just the buffered point
            terra::plotRGB(cropped_plot, r = 1, g = 2, b = 3, stretch = 'lin') #plot the cropped plot
            terra::plot(subplots_in_plot, add = TRUE, lwd = 3, col = 'white')
            
            }
        
      }
    }
}

            

VI_sampling_plot <- function(which_site) {
  input <- VI_list %>% na.omit()
  
  plots_shp <- terra::vect('/Users/khuelsma/Desktop/NEON Spectral Variability/Relevant NEON Materials/NEON TOS Plots/All_NEON_TOS_Plot_Polygons_V11.shp')
  site_plot_polygons <- subset(plots_shp, 
                               stringr::str_detect(plots_shp$plotID, which_site ) & 
                                 stringr::str_detect(plots_shp$appMods, 'div'))
  subplots_shp <- terra::vect(
    '/Users/khuelsma/Desktop/NEON Spectral Variability/Relevant NEON Materials/NEON TOS Plots/All_NEON_TOS_Plot_Subplots_V11.shp')
  site_subplot_pts <- subset(subplots_shp, 
                             stringr::str_detect(subplots_shp$plotID, which_site) & 
                               stringr::str_detect(subplots_shp$appMods, 'div'))
  #for each tile from input
  foreach(
    img = 1:length(input)) %do% {
      this_tile <- input[[img]]
      tile_extent <- ext(this_tile)
      
      plots_projected <- terra::project(site_plot_polygons, this_tile) #can also use David's data here
      plots_in_tile <- plots_projected[tile_extent]
      subplots_projected <- terra::project(site_subplot_pts, this_tile)
      subplots_in_tile <- subplots_projected[tile_extent]
      
      if (nrow(plots_in_tile) > 0) {
        foreach( 
          plot = 1:nrow(plots_in_tile)) %do% {
            
            this_plot <- plots_in_tile[plot,] #for each plot, numbered 1 to nrow
            plot_ext <- ext(this_plot)
            subplots_in_plot <- subplots_in_tile[plot_ext]
            
            this_plot_df <- as.data.frame(this_plot)
            plotID = this_plot_df$plotID
            
            buff_extent <- plot_ext
            cropped_plot <- terra::crop(this_tile, buff_extent) #crop the map to just the buffered point
            terra::plotRGB(cropped_plot, r = 1, g = 2, b = 3, stretch = 'lin') #plot the cropped plot
            terra::plot(subplots_in_plot, add = TRUE, lwd = 3, col = 'white')
          }
        
      }
    }
}
rgb_map_sampling_plot('CPER')
VI_sampling_plot( 'CPER')
