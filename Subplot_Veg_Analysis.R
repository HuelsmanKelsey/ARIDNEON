
#how to combine vegetation details and spectra:

# 1) get perccov_spectra: 
#req raster inputs; choose a cover type, choose a size to focus on subplots or whole plot if a size isn't given.

get_perccov_spectra <- function(which_site, which_year, which_rast, which_size, which_cover) {
  extracted <- extract_subplot_spectra(which_rast, which_site, which_size)
  perccov <- make_perccov_df(which_site, which_year, which_cover, which_size)
  extracted_perccov <- extracted %>%
    left_join(perccov)%>%
    na.omit() %>%
    distinct() %>%
    print()
}
# 2) get_PA_spectra:
# req raster inputs; choose a cover type, choose a size to focus on subplots or whole plot if a size isn't given.

#presence / absence in each subplot or whole plot
get_PA_spectra <- function(which_site, which_year, which_rast, which_size, which_species) {
  #use extract_subplot_spectra with hsilist + which_size;
  extracted <- extract_subplot_spectra(which_rast, which_site, which_size)
  #make a P/A dataframe for the species of interest + which_size for every subplot
  PA <- make_PA_df(which_site, which_year, which_species, which_size, matrix = FALSE)
  
  extracted_PA <- extracted %>%
    left_join(PA) %>%
    na.omit() %>%
    distinct() %>%
    print()
}

PA_full <- make_PA_df('CPER', 2024, 'Opuntia polyacantha')
PA_1m <- make_PA_df('CPER', 2024, 'Opuntia polyacantha', 1)
PA_10m <- make_PA_df('CPER', 2024, 'Opuntia polyacantha', 10)
PA_100m <- make_PA_df('CPER', 2024, 'Opuntia polyacantha', 100)

#presence / absence
OPPO_PA_1 <- get_PA_spectra('CPER', 2024, hsi_list,  1, 'Opuntia polyacantha')
BOGR_PA_1 <- get_PA_spectra('CPER', 2024, hsi_list,  1, 'Bouteloua gracilis')
BODA_PA_1 <- get_PA_spectra('CPER', 2024, hsi_list,  1, 'Bouteloua dactyloides')

#whole plot:
OPPO_PA_full <- get_PA_spectra('CPER', 2024, hsi_list, 'Opuntia polyacantha')

#% cover
OPPO_cov_1 <- get_perccov_spectra('CPER', 2024, hsi_list,  1, 'Opuntia polyacantha')
BOGR_cov_1 <- get_perccov_spectra('CPER', 2024, hsi_list,  1, 'Bouteloua gracilis')
BODA_cov_1 <- get_perccov_spectra('CPER', 2024, hsi_list,  1, 'Bouteloua dactyloides')
soil_cov_1 <- get_perccov_spectra('CPER', 2024, hsi_list,  1, 'soil')

#whole plot 
OPPO_cov_full <- make_perccov_df('CPER', 2024, 'Opuntia polyacantha')
BOGR_cov_full <- make_perccov_df('CPER', 2024, 'Bouteloua gracilis')
soil_cov_full <- make_perccov_df('CPER', 2024, 'soil')

ggplot(soil_cov_1) +
  #geom_density(aes(mean_cov), lwd = 2) +
  geom_histogram(aes(), bins = 20, alpha = 0.5) +
  theme_classic()
#shows us that most 1m2 plots have 10-15% soil cover (nearly 1000?)
#second most: 5-10, then 20-25


soil_cov_1 %>%
  filter(!is.na(mean_cov)) %>%
  pivot_longer(cols = unique(kept_bands$band), names_to = 'band', values_to = 'refl') %>%
  left_join(kept_bands) %>%
  group_by(eventID, plotID, subplotID, nlcdClass, easting, northing, latitude, longitude, horzUncert,
           subpltSize, refPoint, band, year, wv, mean_cov) %>%
  summarise(mean_refl = mean(refl),
            sd_refl = sd(refl),
            se_refl = sd(refl)/sqrt(n()),
            CV_refl = sd_refl/mean_refl) %>%
  mutate(cov_cat = case_when(mean_cov > 70 ~ 'high',
                             mean_cov <= 70 & mean_cov >= 10 ~ 'mod',
                             mean_cov < 10 ~ 'low')) %>%
  ggplot(aes(x = wv, y = mean_refl)) +
  geom_point(aes(colour = cov_cat)) +
  geom_errorbar(aes(x = wv, 
                    ymin = mean_refl + 2*se_refl,
                    ymax = mean_refl + 2*se_refl,
                    colour = cov_cat))+
  theme_classic()


#make long
  extracted_PAlong <- extracted_PA %>%
    pivot_longer(cols = unique(kept_bands$band), names_to = 'band', values_to = 'refl') %>%
    left_join(kept_bands) %>%
    select(eventID, plotID, subplotID, SOI, nlcdClass,
           latitude, longitude, easting, northing, horzUncert, 
           subpltSize, refPoint, band, refl, year, wv) %>%
    mutate(eventID = paste0(plotID, '_', subplotID, '_', year),
           parentEventID = paste0(plotID, '_', year)) %>%
    #group by eventID first, 
    group_by(SOI, eventID, wv) %>%
    summarise(mean_subplot_refl = mean(refl),
              sd_subplot_refl = sd(refl),
              se_subplot_refl = sd(refl)/sqrt(n()),
              CV_subplot_refl = sd_subplot_refl/mean_subplot_refl)
  
  return(extracted_PAlong)
}

extracted_PAlong <- get_PA_spectra(which_site, which_year, which_size, which_species)

#every subplot's reflectance
extracted_PAlong %>%
  mutate(SOI = as.factor(SOI)) %>%
  ggplot(aes(x = wv, y = mean_subplot_refl)) +
  geom_point(aes(colour = SOI)) +
  geom_errorbar(aes(x = wv, 
                    ymin = mean_subplot_refl + 2*se_subplot_refl,
                    ymax = mean_subplot_refl + 2*se_subplot_refl,
                    colour = SOI))+
  theme_classic()



#summarize for each plot
summary_soil_spec_1m <- soil_spec_1m %>%
  dplyr::mutate(soil_cov_cat = case_when(percentCover > 75 ~ 'high',
                                         percentCover <= 75 & percentCover >= 30 ~ 'mod',
                                         percentCover < 30 ~ 'low')
                ) %>%
  group_by(soil_cov_cat, wv) %>%
  summarise(mean_refl = mean(refl),   
            sd_refl = sd(refl),
            se_refl = sd_refl/sqrt(n()),
            CV_refl = sd_refl/mean_refl) %>%
  ggplot(aes(x = wv, y = mean_refl)) +
  geom_point(aes(colour = soil_cov_cat)) +
  geom_errorbar(aes(x = wv, 
                    ymin = mean_refl + 2*se_refl,
                    ymax = mean_refl + 2*se_refl,
                    colour = soil_cov_cat))+
  theme_classic()

summary_soil_spec_1m


#PLSR:
soil_spec_1m





#train on the 8 x 1m2 subplots in each plot: 
plot_summary_1m <- soil_spec_1m %>%
  group_by(plotID, wv) %>%
  summarise(n = n_distinct(subplotID))

  
#can add in details from RGB imagery: 


#generate a field-based map of something:

#generate a map from RGB imagery of site or subplot:

