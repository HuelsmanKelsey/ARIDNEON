
# OTHER DETAILS: SPECTRA --------------------------------------------------
#Reminder to get USGS Endmembers
# OTHER DETAILS: SPECTRA --------------------------------------------------
#Reminder to get USGS Endmembers
EMs <- read.csv('/Users/khuelsma/Downloads/USGS_grass_endmems.csv')
#bad bands will have negative reflectance
EMs %>%
  filter(refl > 0) %>%
  ggplot(aes(x = wv_nm, y = refl)) +
  geom_point(aes(colour = endmem)) +
  theme_classic()

# David Augustine's data:
#cover for 48 polygons from a few days in June 2026:
#plot_surveys <- read.csv('Augustine_Plot_Surveys.csv')
#pure shrub and bare soil points and polygons
#pure_ems <- read.csv('Augustine_Pure_EMs.csv')


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



# Begin Workflow: ---------------------------------------------------------

# load what I already extracted OR 
CPER_1_10 <- readr::read_csv('/Users/khuelsma/CPER_2024_extractedplots1to10.csv')
CPER_11_15 <- readr::read_csv('/Users/khuelsma/CPER_2024_extractedplots11to15.csv')
CPER_16_23 <- readr::read_csv('/Users/khuelsma/CPER_2024_extractedplots16to23.csv')


# Choose a site/location and open or request from API to download
which_site <- 'CPER'

# you will need a NEON token to use the API for this command:
neon_token <- if (file.exists('NEON_token.txt')) {
  neon_token <- readLines('neon_token.txt')
} else {
  neon_token <- rstudioapi::askForPassword(prompt = 'enter NEON token')
}
options(timeout = 3600) #increase timeout to 1 hour

#this gives us a timeline of available data and folder numbers for opening the HSI with a little more control
both <- both_avail(which_site)
both$folder = 1:nrow(both)

#set the site, get all plots and associated tiles at the site
tile_plots <- get_site_EN(which_site, tiles_only = FALSE)

filename <- paste0(which_site, '_filelist_', which_year, '.csv')
#get a list of the files to open (download them if needed)

#if we have a list of files, open them
if (file.exists(filename)) {
  filelist <- read_csv(filename)
  #create separate lists of rgb (rgb_path) vs. hsi (path) 
  rgb_files <- filelist$rgb_path
  hsi_files <- filelist$path
  
}

#if we don't, see if we have the files but not the list.
if (!file.exists(filename)) {
  #check files exist and make it
  site_files <- check_files_exist(which_site)
  filelist <- site_files %>% filter(!is.na(path))

  if (nrow(files) == 0) {
    print('download data')
  }
  
  if (nrow(files) > 0) {
  return(site_files)
  }
}

#using files list: 1) make maps, 2) extract data:
maplist <- map_site_tiles(files)
#save extracted data
extracted_data <- extract_tile_plots(files)
write.csv(extracted_data, file = paste0(which_site, '_refl.csv'))


#start w rgb tiles: for each list item (file path), create a list of tilemaps to merge
tilemap_rgb_list <- foreach (img = rgb_files_2024) %do% {
  tilemap <- terra::rast(img)
  #terra::plotRGB(tilemap)
}
#full site map (where sampling occurred)
full_site_2024 <-  terra::merge(terra::sprc(tilemap_rgb_list))
#terra::plotRGB(full_site_2021)

#terra::plot(LL_poly$geometry, add = TRUE)
#terra::crs(LL_poly)

#make a full site mosaic:
site_hsi_mosaic <- terra::merge(terra::sprc(maplist))
plot(site_hsi_mosaic)

# each map needs: a) bandmath and b) plot extraction. and saving as a png at each step.
maplist



# 2. Do the band math on the FULL SITE mosaic
site_ndvi <- (site_hsi_mosaic[["NIR"]] - site_hsi_mosaic[["red"]]) / (site_hsi_mosaic[["NIR"]] + site_hsi_mosaic[["red"]] + 0.0001)
site_nirv <- site_hsi_mosaic[["NIR"]] * site_ndvi
site_pri <- (site_hsi_mosaic[["PRI"]] - site_hsi_mosaic[["green"]]) / (site_hsi_mosaic[["PRI"]] + site_hsi_mosaic[["green"]])
site_cci <- (site_hsi_mosaic[["PRI"]] - site_hsi_mosaic[["red"]]) / (site_hsi_mosaic[["PRI"]] + site_hsi_mosaic[["red"]])

site_indices <- c(site_cci, site_nirv, site_pri)
names(site_indices) <- c("CCI", "NIRv", "PRI")

# 3. Export the Full Site Images
# Note: Making the width/height larger here since it covers the whole site
# Full Site Pseudo-RGB
png(file.path(out_dir, paste0(which_site, "_", which_year, "_FullSite_pseudoRGB.png")), width = 3000, height = 3000)
terra::plotRGB(site_hsi_mosaic[[c("red", "green", "blue")]], r=1, g=2, b=3, stretch="lin", main=paste(which_site, which_year, "Full Site Pseudo-RGB"))
dev.off()

# Full Site Individual Indices
for (idx in names(site_indices)) {
  png(file.path(out_dir, paste0(which_site, "_", which_year, "_FullSite_", idx, ".png")), width = 3000, height = 3000)
  plot(site_indices[[idx]], main = paste(which_site, which_year, "Full Site -", idx))
  dev.off()
}

# Full Site Composite RGB Indices (R = CCI, G = NIRv, B = PRI)
png(file.path(out_dir, paste0(which_site, "_", which_year, "_FullSite_Indices_RGB.png")), width = 3000, height = 3000)
terra::plotRGB(site_indices, r=1, g=2, b=3, stretch="lin", main=paste(which_site, which_year, "Full Site\n(R=CCI, G=NIRv, B=PRI)"))
dev.off()

#save all tiles and plots:
out_dir <- paste0(which_site, "_", which_year, "_Maps")
if(!dir.exists(out_dir)) dir.create(out_dir)
for (tile_name in names(maplist)) {
  parts <- strsplit(tile_name, "_")[[1]]
  t_year <- as.numeric(parts[1])
  t_E <- as.numeric(parts[2])
  t_N <- as.numeric(parts[3])
  hsi_rast <- maplist[[tile_name]]
  
  ndvi <- (hsi_rast[["NIR"]] - hsi_rast[["red"]]) / (hsi_rast[["NIR"]] + hsi_rast[["red"]] + 0.0001)
  nirv <- hsi_rast[["NIR"]] * ndvi
  pri <- (hsi_rast[["PRI"]] - hsi_rast[["green"]]) / (hsi_rast[["PRI"]] + hsi_rast[["green"]])
  cci <- (hsi_rast[["PRI"]] - hsi_rast[["red"]]) / (hsi_rast[["PRI"]] + hsi_rast[["red"]])
  
  indices <- c(cci, nirv, pri)
  names(indices) <- c("CCI", "NIRv", "PRI")
  
  tile_file_info <- site_files %>% filter(year == t_year & E_tile == t_E & N_tile == t_N)
  
  if (nrow(tile_file_info) > 0 && !is.na(tile_file_info$rgb_path[1])) {
    fine_rgb <- terra::rast(tile_file_info$rgb_path[1])
  } else {
    fine_rgb <- NULL
  }
  cat("Saving Tile Maps:", tile_name, "\n")
  
  png(file.path(out_dir, paste0(tile_name, "_pseudoRGB.png")), width = 1000, height = 1000)
  terra::plotRGB(hsi_rast[[c("red", "green", "blue")]], r=1, g=2, b=3, stretch="lin", main=paste("Tile", tile_name, "Pseudo-RGB"))
  dev.off()
  
  for (idx in names(indices)) {
    png(file.path(out_dir, paste0(tile_name, "_", idx, ".png")), width = 1000, height = 1000)
    plot(indices[[idx]], main = paste("Tile", tile_name, "-", idx))
    dev.off()
  }
  
  png(file.path(out_dir, paste0(tile_name, "_Indices_RGB.png")), width = 1000, height = 1000)
  terra::plotRGB(indices, r=1, g=2, b=3, stretch="lin", main=paste("Tile", tile_name, "\n(R=CCI, G=NIRv, B=PRI)"))
  dev.off()
  
  if (!is.null(fine_rgb)) {
    png(file.path(out_dir, paste0(tile_name, "_10cm_RGB.png")), width = 1000, height = 1000)
    terra::plotRGB(fine_rgb, r=1, g=2, b=3, stretch="lin", main=paste("Tile", tile_name, "10cm RGB"))
    dev.off()
  }
  
  thisfile_plots <- tile_plots %>% filter(E_tile == t_E & N_tile == t_N)
  
  if(nrow(thisfile_plots) > 0) {
    for(t in 1:nrow(thisfile_plots)) {
      thisplot <- thisfile_plots[t,]
      plotID <- thisplot$plotID
      
      cat("  - Extracting Plot:", plotID, "\n")
      thisplot_vect <- vect(thisplot, geom = c('easting', 'northing'), crs = crs(hsi_rast))
      square_plot <- as.polygons(ext(buffer(thisplot_vect, width = 10)))
      
      plot_hsi <- crop(hsi_rast, square_plot)
      plot_indices <- crop(indices, square_plot)
      
      png(file.path(out_dir, paste0(plotID, "_", t_year, "_pseudoRGB.png")), width = 400, height = 400)
      terra::plotRGB(plot_hsi[[c("red", "green", "blue")]], r=1, g=2, b=3, stretch="lin", main=paste(plotID, "Pseudo-RGB"))
      dev.off()
      
      for (idx in names(plot_indices)) {
        png(file.path(out_dir, paste0(plotID, "_", t_year, "_", idx, ".png")), width = 400, height = 400)
        plot(plot_indices[[idx]], main = paste(plotID, "-", idx))
        dev.off()
      }
      
      if(!is.null(fine_rgb)) {
        plot_fine_rgb <- crop(fine_rgb, square_plot)
        png(file.path(out_dir, paste0(plotID, "_", t_year, "_10cm_RGB.png")), width = 400, height = 400)
        terra::plotRGB(plot_fine_rgb, r=1, g=2, b=3, stretch="lin", main=paste(plotID, "10cm RGB"))
        dev.off()
      }
    }
  }
}
cat("All maps successfully exported to", out_dir, "\n")



#can we locate the 1 m2 subplots?


#train on the 8 x 1m2 subplots: 

#can add in details from RGB imagery: 

#generate a field-based map of something:
#generate a map from RGB imagery of site or subplot:

#train HS on just 1 m2 data

#try to geolocate subplots and look at their reflectances:

subplots_shp <- sf::st_read('/Users/khuelsma/Desktop/NEON Spectral Variability/Relevant NEON Materials/NEON TOS Plots/All_NEON_TOS_Plot_Subplots_V11.shp')

plot(subplots_shp$geometry)

subplots_1m <- subplots_shp %>%
  filter(subpltDim == "1m x 1m")
subplots_10m <- subplots_shp %>%
  filter(subpltDim == "10m x 10m")

subplots_shp %>%
  group_by(plotID, subplotID, easting, northing, horzUncert, nlcdClass, subpltSize, geometry) %>%
  filter(grep('CPER', plotID) == TRUE)

CPER_subs %>%
  group_by(plotID, boutNumber, divDataType, otherVariables) %>%
  summarise(n_subplots = n_distinct(subplotID),
            perc_cover = mean(percentCover),
            sd_perc_cover = sd(percentCover)) %>%
  print(n = 100)

ARIDveg <- neon_veg %>%
  mutate(year = lubridate::year(lubridate::ymd(endDate))) %>%
  filter(siteID == 'SRER') #JORN and SRER later

# deal w 1 m subplots later <- read.csv('NEON_veg_subplots_2026.csv')
common_sp <- ARIDveg %>%
  group_by(scientificName, year) %>%
  summarise(nplots = n_distinct(plotID)) %>%
  group_by(scientificName) %>%
  summarise(years = n_distinct(year),
            mean_nplots = mean(nplots)) %>%
  filter(mean_nplots > 10)
common_sp %>% print(n = 70)
```

We can use the common species list (common_sp) to get a feel for how to characterize veg communities. For each plot, see what from the list is present.

Need to filter out Unknown and do other clean up things (sp.)

```{r}
common <- ARIDveg %>%
  mutate(keep = scientificName %in% common_sp$scientificName) %>%
  filter(keep == TRUE) %>%
  group_by(plotID, year, nlcdClass, family, nativeStatusCode) %>%
  #we can do 'count' to get a rough estimate of its common-ness in plot.
  count(scientificName) %>%
  print()
```

What species are in what plots?
  ```{r}

plotsp <- common %>% 
  #filter to a single year
  filter(year == 2016) %>%
  group_by(scientificName, nativeStatusCode, family) %>%
  summarise(nplots = n_distinct(plotID),
            which_plots = list(unique(plotID)))


#wow ok all the common species in are native?? 
#let's take a matrix from a single year with a row for every plot and see how different they are:
#before I figure out scientific name simplification I'm just going to do family
common_mat <- common %>%
  filter(year == 2016) %>%
  ungroup() %>%
  select(plotID, family, n) %>%
  group_by(plotID, family) %>%
  summarise(N = sum(n)) %>%
  tidyr::pivot_wider(names_from = family, values_from = N, values_fill = 0)

dist_matrix <- common_mat %>%
  dist() %>%
  as.matrix()

#this is a visual of the plot similarities / differences
dist_matrix %>%
  as.data.frame() %>%
  mutate(R = rownames(dist_matrix)) %>%
  tidyr::pivot_longer(1:33) %>%
  filter(value > 20) %>%
  distinct(R, name, value) %>%
  ggplot() +
  geom_point(aes(x = R, y = name, colour = log(value)), size = 3) +
  theme_classic()

#20 and 23 are rather different at SRER; 
#early 20s and early 30s are all kind of diff from same ones.

common %>%
  filter(plotID == 'SRER_002') %>%
  filter(year == 2016) %>%
  filter(n > 5)

#SRER 1: Zinnia acerosa, Muhlenbergia porteri, Setaria leucopila
#SRER 2: Chamaesyce arizonica, Chamaesyce florida
#SRER 20: Evolvulus arizonicus, Prosopis velutina from the pea fam
#SRER 23: Cylindropuntia acanthocarpa, Opuntia engelmannii, Evolvulus arizonicus, Prosopis velutina,
# Eragrostis lehmanniana, Aristida ternipes, Digitaria californica, Heteropogon contortus
```