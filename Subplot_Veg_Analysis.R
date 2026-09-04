
#can we locate the 1 m2 subplots?

#train on the 8 x 1m2 subplots: 

#can add in details from RGB imagery: 

#generate a field-based map of something:
#generate a map from RGB imagery of site or subplot:

#train HS on just 1 m2 data

#try to geolocate subplots and look at their reflectances:

subplots_1m <- sf::st_read('/Users/khuelsma/Desktop/NEON Spectral Variability/Relevant NEON Materials/NEON TOS Plots/All_NEON_TOS_Plot_Subplots_V11.shp')
subplots_1m <- subplots_1m %>%
  filter(subpltDim == "1m x 1m") %>%
  filter(siteID == 'CPER')

CPER_subs %>%
  group_by(plotID, boutNumber, divDataType, otherVariables) %>%
  summarise(n_subplots = n_distinct(subplotID),
            perc_cover = mean(percentCover),
            sd_perc_cover = sd(percentCover)) %>%
  print(n = 100)
  
