
# Overview: ---------------------------------------------------------------

#Get Veg lets you...
# Choose a site/location and filter the vegetation dataset
# Crosswalk / interpret the data
# Visualize and understand the distributions 

# Summarize: 
# Whole plot, subplot size-specific analyses of % cover
# relative frequency for plot and larger subplots

# capturing vegetation “vibe” using PA or % cover
 

# Load packages -----------------------------------------------------------

#basic packages
library('readr')
library('tidyr') 
library('dplyr')
library('ggplot2')
library('lubridate')
library('stringr')

#spatial data
library('sf')                 
library('terra')    
library('mapview')

#for NEON API:
library('neonUtilities')    
library('jsonlite')   

# reading HDF5 files (AOP data)
library('rhdf5')              

# to extract metadata from shape files
library('xml2')

#for vegetation data
library('vegan')

#for parallelization
library('doParallel')
library('foreach')

# Workflow:  --------------------------------------------------------------
# Choose a site/location and filter the vegetation dataset
repo_dir <- '/Users/khuelsma/Desktop/ARIDNEON/'
setwd(repo_dir)

#plots and subplots from ARID sites: CPER, RMNP, STER, JORN, SRER
ARID_veg <- read.csv('NEON_ARID_veg.csv')
subplots <- read.csv('NEON_ARID_subs.csv')
#subplots from CPER only
CPER_subs <- subplots %>%
  filter(siteID == 'CPER')
CPER_veg <- ARID_veg %>%
  filter(siteID == 'CPER')

# Crosswalk / interpret the data
# Creating a single Event dataframe ---------------------------------------
Event_subplots <- subplots %>%
  select(domainID, 
         siteID, 
         endDate, 
         decimalLatitude, 
         decimalLongitude,
         geodeticDatum, 
         coordinateUncertainty, 
         elevation, 
         elevationUncertainty,
         nlcdClass, 
         plotID, 
         subplotID, 
         boutNumber, 
         release, 
         publicationDate, 
         samplingImpractical,
         samplingImpracticalRemarks, 
         samplingProtocolVersion)
# ^ this gets rbinded to the next section
Event <- ARID_veg %>%
  #create variables that are in the subplot dataset so we can combine the two
  mutate(divDataType = NA,
         otherVariablesPresent = NA,
         otherVariables = NA,
         percentCover = NA) %>%
  select(names(Event_subplots)) %>%
  rbind(Event_subplots) %>% 
  
  #now Event_subplots is in here, so we start renaming, etc.
  dplyr::rename(
    coordinateUncertaintyInMeters = coordinateUncertainty,
    minimumElevationInMeters = elevation,
    samplingProtocol = samplingProtocolVersion
  ) %>%
  dplyr::mutate(
    maximumElevationInMeters = minimumElevationInMeters,
    eventDate = lubridate::date(endDate),
    year = format(as.Date(endDate, format="%Y-%m-%d"),"%Y"),
    #location = siteName, #domainID,
    locationID = siteID,
    #event ID will indicate both subplot and bout; if there are multiple bouts, can take annual max.
    eventID = paste0(plotID, '_', subplotID, '_', year, '_', boutNumber),
    eventRemarks = ifelse(is.na(samplingImpractical), NA, paste('sampling issues:', samplingImpractical, ':', samplingImpracticalRemarks)),
    siteNumber = paste0(domainID, '_', siteID),
    habitat = paste("NLCD:", nlcdClass, sep = " "),
    sampleSizeValue = case_when(
      stringr::str_detect(subplotID, '_100') ~ 100,
      stringr::str_detect(subplotID, '_10_') ~ 10,
      stringr::str_detect(subplotID, '_1_') ~ 1,
      stringr::str_detect(subplotID, '_400') ~ 400),
    sampleSizeUnit = 'm',
    #samplingEffort = ' ' #figure out later
    
    #parent_event for annual survey (for combining bouts)
    parentEventID = paste0(plotID, '_', subplotID, '_', year),
    
    #adding all smaller scale occurrences to broader ones
    contained_within_100 = case_when(
      sampleSizeValue == 100 ~ paste0(plotID, '_', subplotID, '_', year),
      sampleSizeValue == 10 ~ paste0(plotID, '_', paste0(strtrim(subplotID, 7), 0), '_', year),
      sampleSizeValue == 1 ~ paste0(plotID, '_', paste0(strtrim(subplotID, 5), 00), '_', year)),
    contained_within_10 = case_when(
      sampleSizeValue == 10 ~ paste0(plotID, '_', subplotID, '_', year),
      sampleSizeValue == 1 ~ paste0(plotID, '_', paste0(strtrim(subplotID, 5), 0), '_', year))
  ) %>% 
  select(
    siteNumber,
    locationID, 
    habitat,
    parentEventID, #temporarily set to single year
    year,
    eventDate, 
    contained_within_10, 
    contained_within_100, 
    eventID, 
    eventRemarks,
    sampleSizeValue, 
    sampleSizeUnit,
    decimalLatitude, 
    decimalLongitude,
    geodeticDatum, 
    maximumElevationInMeters,
    coordinateUncertaintyInMeters,
    minimumElevationInMeters,
    samplingProtocol
  )

# Creating a single Occurrence dataframe ---------------------------------------
Occurrence_subplots <- subplots %>%
  select(uid, endDate, plotID, subplotID, boutNumber, release, publicationDate, 
         samplingImpractical, samplingImpracticalRemarks, samplingProtocolVersion,
         otherVariablesPresent, divDataType, otherVariables, percentCover,  #divDataType is plantSpecies, otherVariables
         heightPlantOver300cm, morphospeciesID, morphospeciesIDRemarks, 
         taxonID, taxonIDRemarks, remarks, nativeStatusCode, 
         taxonRank, family, scientificName, targetTaxaPresent, heightPlantOver300cm,
         heightPlantSpecies, recordedBy, measuredBy, identificationHistoryID, identificationReferences, identificationQualifier)
# ^ this gets rbinded in the next section
Occurrence <- ARID_veg %>%
  mutate(divDataType = NA,
         otherVariablesPresent = NA,
         otherVariables = NA,
         percentCover = NA,
         heightPlantOver300cm = NA, 
         heightPlantSpecies = NA) %>%
  select(names(Occurrence_subplots)) %>%
  rbind(Occurrence_subplots) %>% 
  # dplyr::filter(
  #   targetTaxaPresent == "Y" &
  #     !is.na(scientificName)
  # ) %>%
  dplyr::rename(
    occurrenceID = uid #debatable
  ) %>%
  mutate(
    eventDate = lubridate::date(endDate),
    year = format(as.Date(endDate, format="%Y-%m-%d"),"%Y"),
    #event ID will indicate both subplot and bout; if there are two bouts, can take max.
    #eventID should align with the actual event
    eventID = paste0(plotID, '_', subplotID, '_', year, '_', boutNumber),
    genus = ifelse ((taxonRank == "genus" | 
                       taxonRank == "species" |
                       taxonRank == "variety" | 
                       taxonRank == "subspecies" | 
                       taxonRank == 'speciesGroup'),
                    word(scientificName, start = 1), NA),
    speciesEpithet = ifelse ((taxonRank == "species" |
                                taxonRank == "variety" | 
                                taxonRank == "subspecies"| 
                                taxonRank == 'speciesGroup'),
                             word(scientificName, start = 2), NA),
    establishmentMeans = nativeStatusCode %>%
      recode(
        N = "native",
        I = "introduced",
        UNK = "uncertain", #I replaced NA here
        NI = "native & introduced", #some infrataxa are native and others are introduced
        A = "", #Presumed absent... shouldn't be any cases of this?
        'I?' = "uncertain, likely introduced",
        'N?' = "uncertain, likely native",
        'NI?' = "uncertain, likely native & introduced"),
    taxonRank = case_when(
      taxonRank == 'species' ~ 'species',
      taxonRank == 'speciesGroup' ~ 'species', #added with 2025
      taxonRank == 'variety' ~ 'species',
      taxonRank == 'subspecies' ~ 'species', 
      taxonRank == 'form' ~ 'species',
      taxonRank == 'genus' ~ 'genus',
      taxonRank == 'family' ~ 'family', 
      taxonRank == 'order' ~ 'order',
      taxonRank == 'class' ~ 'class',
      taxonRank == 'unranked' ~ 'unranked',
      taxonRank == 'phylum' ~ 'phylum',
      taxonRank == 'kingdom' ~ 'kingdom'),
    scientificName = case_when(     # scientific name to species level (i.e. simplified if sub-species)
      #if taxonomiccertainty == species, take first two words 
      taxonRank == 'species' | 
        taxonRank == "variety" | 
        taxonRank == "subspecies"| 
        taxonRank == 'speciesGroup' ~ stringr::word(scientificName, 1, 2),
      #if taxonomiccertainty == genus, take the first word and add sp.
      taxonRank == 'genus' ~ paste(stringr::word(scientificName, 1), 'sp.')),
    verbatimIdentification = taxonIDRemarks,
    fieldNotes = remarks,
    basisOfRecord = 'HumanObservation',
    
    #plant % covers
    organismQuantity = ifelse(divDataType == 'plantSpecies', percentCover, NA),
    organismQuantityType = ifelse(divDataType == 'plantSpecies', '% cover', NA),
    #measurementType will be either otherVariables (e.g., litter, etc., or understory height)
    measurementType = case_when(
      # non vascular plant % cover
      # although "overstory" is a category, which is plant species.
      divDataType == 'otherVariables' ~ '% cover',
      divDataType == 'plantSpecies' & heightPlantOver300cm == 'N' ~ 'understory height',
      divDataType == 'plantSpecies' & heightPlantOver300cm == 'Y' ~ 'overstory % cover'),
    verbatimMeasurementType = ifelse(divDataType == 'otherVariables', 'percentCover', 
                                     ifelse(divDataType == 'plantSpecies', 'heightPlantSpecies', NA)),
    measurementValue = ifelse(divDataType == 'otherVariables', percentCover, 
                              ifelse(divDataType == 'plantSpecies', heightPlantSpecies, percentCover)),
    measurementUnit = ifelse(divDataType == 'otherVariables', '%',
                             ifelse(divDataType == 'plantSpecies', 'cm', '%')),
    #detailed version
    coverType = case_when(divDataType == 'plantSpecies' & taxonRank == 'species' ~ scientificName,
                          divDataType == 'plantSpecies' & taxonRank == 'genus' ~ scientificName,
                          divDataType == 'plantSpecies' & taxonRank != 'species' & taxonRank != 'genus' ~ 'plants',
                          divDataType == 'otherVariables' & otherVariables != 'overstory' ~ otherVariables,
                          divDataType == 'otherVariables' & otherVariables == 'overstory' ~ NA),
    coverLocation = case_when(divDataType == 'plantSpecies' & heightPlantOver300cm == 'Y' ~ 'overstory',
                              divDataType == 'plantSpecies' & heightPlantOver300cm == 'N' ~ 'understory',
                              divDataType == 'plantSpecies' & is.na(measurementValue) ~ 'ground',
                              divDataType == 'otherVariables' & coverType != 'overstory' ~ 'ground'),
    #relabeling standing dead and biocrust specific categories into just standing dead and biocrust.
    coverTypeGeneral = case_when(divDataType == 'plantSpecies' ~ 'plants', 
                                 divDataType == 'otherVariables' & 
                                   !str_detect(otherVariables, 'standingDead') & 
                                   !str_detect(otherVariables, 'biocrust') ~ otherVariables,
                                 divDataType == 'otherVariables' & str_detect(otherVariables, 'standingDead') ~ 'standingDead',
                                 divDataType == 'otherVariables' & str_detect(otherVariables, 'biocrust') ~ 'biocrust')) %>%
  print()

# Combining Event and Occurrence ------------------------------------------
#All observations:
ARID_Obs <- Occurrence %>%
  filter(!is.na(coverType)) %>% 
  left_join(Event, 
            by = c('eventID', 'eventDate', 'year'),
            relationship = "many-to-many") %>%
  distinct()

#summary of % cover at each coverlocation (ground, understory, overstory)
tot_cov_by_group <- ARID_Obs %>% 
  group_by(subplotID, plotID, year, boutNumber, coverLocation) %>% 
  summarise(totcov = sum(percentCover)) 
tot_cov_by_subplot <- tot_cov_by_group %>%
  group_by(subplotID, plotID, boutNumber, year) %>%
  summarise(totcov2 = sum(totcov))

#Looking at bouts within a year:
ARID_Obs %>%
  group_by(year, siteNumber) %>%
  summarise(n_bouts = n_distinct(boutNumber)) %>%
  filter(n_bouts > 1) %>% 
  distinct(siteNumber, year, n_bouts)

#ANNUAL Summary: simplifying multiple bouts per year by taking max val, just to characterize overall
ARID_Annual_Obs <- ARID_Obs %>%
  #group by everything but eventID and eventRemarks to get annual obs for subplot
  group_by(locationID, 
           year, 
           parentEventID, 
           sampleSizeValue, 
           #organismQuantityType, 
           coverType,
           coverTypeGeneral,
           measurementType, 
           coverLocation,
           taxonRank,
           scientificName, 
           family, 
           genus, 
           speciesEpithet,
           establishmentMeans,
           contained_within_10, 
           contained_within_100, 
           habitat, 
           decimalLatitude, 
           decimalLongitude, 
           samplingImpractical, samplingImpracticalRemarks,
           identificationHistoryID, verbatimIdentification) %>%
  summarise(n_bouts = n_distinct(eventID),
            n_remarks = n_distinct(na.omit(eventRemarks)),
            remarks = ifelse(n_remarks >= 1, 
                             unique(na.omit(eventRemarks)),
                             NA),
            perccovs = n_distinct(percentCover, na.rm = TRUE),
            mvals = n_distinct(measurementValue, na.rm = TRUE),
            mean_bout = mean(percentCover),
            max_bout = max(percentCover),
            mean_val = mean(measurementValue),
            max_val = max(measurementValue)) %>%
  mutate(percentCover = case_when(
    #one bout value: just take the quantity
    perccovs == 1 ~ mean_bout, #calculated above
    
    #two bout values, take the max:
    perccovs == 2 ~ max_bout),
    measurementValue = case_when(
      mvals == 1 ~ mean_val,
      mvals == 2 ~ max_val
    ),
    #now that that's dealt with, set the parent event back to the full plot and set event to annual survey
    eventID = parentEventID,
    parentEventID = paste0(stringr::str_sub(parentEventID, 1, 8), '_', year)) %>%
  ungroup() %>%
  select(-mvals, -perccovs, -mean_val, -max_val, -max_bout, -mean_bout, -n_bouts, -n_remarks) %>%
  print()

# Visualize and understand the distributions ?

# Summarize: 
# Whole plot: P/A
# whole plot: size-specific presences and % cover (for 1m)
# subplot specific: P/A (each subplot)
# vegetation vibe: PA, % cover

#whole plot PA

#whole plot summary of subplots

#subplot specific pairings

# capturing vegetation "vibe" using PA

# capturing vegetation “vibe” using % cover
# Analysis of 1m plant cover, CPER 2024 ------------------------------------

#this is every single 1m subplot from CPER 2024
cover_1m <- ARID_Annual_Obs %>% #or can do ARID_Obs for multiple bouts!
  filter(sampleSizeValue == 1) %>%
  filter(locationID == 'CPER') %>%
  filter(year == 2024) %>%
  distinct(eventID,
           coverType, coverTypeGeneral, percentCover, 
           contained_within_10, contained_within_100)

#two options: 

#one, group all plants (coverTypeGeneral)
percentCoverGeneral <- cover_1m %>%
  group_by(eventID, coverTypeGeneral) %>%
  summarise(cov = sum(percentCover))
percentCoverGeneral_wide <- percentCoverGeneral %>%
  distinct(eventID, coverTypeGeneral, cov) %>%
  group_by(eventID) %>%
  pivot_wider(names_from = coverTypeGeneral,
              values_from = cov,
              values_fill = 0) %>%
  ungroup() 
mat_gen <- percentCoverGeneral_wide %>%
  select(-eventID) %>%
  as.matrix()
row.names(mat_gen) <- percentCoverGeneral_wide$eventID

#two: don't group all plants and use scientific name
percentCoverSpecific <- cover_1m %>%
  group_by(eventID, coverType) %>%
  summarise(cov = sum(percentCover))
percentCoverSpecific_wide <- percentCoverSpecific %>%
  distinct(eventID, coverType, cov) %>%
  group_by(eventID) %>%
  pivot_wider(names_from = coverType,
              values_from = cov,
              values_fill = 0) %>%
  ungroup()
mat_spec <- percentCoverSpecific_wide %>%
  select(-eventID) %>%
  as.matrix()
row.names(mat_spec) <- percentCoverSpecific_wide$eventID


# PCA: specific and general --------------------------------------------
# Specific cover PCA ------------------------------------------------------
spec_pca <- pca(mat_spec, scale = TRUE)

#each subplot in component space:
subplot_locs <- spec_pca$CA$u[,1:3] %>% #198 x 22
  as.data.frame() %>%
  mutate(subplot = rownames(spec_pca$CA$u)) %>%
  pivot_longer(PC1:PC3, names_to = 'PC', values_to = 'subplot_loc')
#each species
sp_load <- spec_pca$CA$v[,1:3] %>% #22 x 22
  as.data.frame() %>%
  mutate(covtype = rownames(spec_pca$CA$v)) %>%
  pivot_longer(PC1:PC3, names_to = 'PC', values_to = 'PC_val')
loads_locs <- sp_load %>%
  left_join(subplot_locs) %>%
  mutate(subplot_sp_relationship = PC_val*subplot_loc) %>%
  mutate(rel_cat = case_when(subplot_sp_relationship < -0.05 ~ 'strong negative',
                             subplot_sp_relationship > -0.05 & subplot_sp_relationship < 0.05 ~ 'neutral',
                             subplot_sp_relationship > 0.05 ~ 'strong positive'))
#how each covertype loads in the first 3 components of PCA:
loads_locs %>%
  filter(PC_val > 0.2 | PC_val < -0.2) %>%
  ggplot() +
  geom_point(aes(x = covtype, y = PC_val, colour = PC)) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  facet_wrap(~PC, scales = 'free_x') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#PC1: + standing dead herbaceous; - soil; - Salsola, Plantago, Amaranthera?
#PC2: only positive: litter, Sisymbrium altissium; less so vulpia octoflora, Elymus elymoides, Lappula occidentalis, Lepidium densiflorum.
#PC3: most negative: lappula occidentalis, Sisymbrium altissium. lepidium densiflorum

# + PC1: hella standing dead herbaceous, little soil, Salsola, Plantago, Amaranthera.
# - PC1: hella soil and little standing dead herbaceous. Maybe also hella Salsola, Plantago, and Amaranthera?

# + PC2: lots of litter and sisymbrium
# - PC2: little litter and sisymbrium

# + PC3: little lap occ, sisymbrium altissium, and lepidium densiflorum.
# - PC3: lots of lap occ, sisymbrium altissium, and lepidium densiflorum.

# General cover PCA -------------------------------------------------------
gen_pca <- pca(mat_gen, scale = TRUE)

#each subplot in component space:
subplot_locs <- gen_pca$CA$u[,1:3] %>% #198 x 22
  as.data.frame() %>%
  mutate(subplot = rownames(gen_pca$CA$u)) %>%
  pivot_longer(PC1:PC3, names_to = 'PC', values_to = 'subplot_loc')
#each species
sp_load <- gen_pca$CA$v[,1:3] %>% #22 x 22
  as.data.frame() %>%
  mutate(covtype = rownames(gen_pca$CA$v)) %>%
  pivot_longer(PC1:PC3, names_to = 'PC', values_to = 'PC_val')
loads_locs <- sp_load %>%
  left_join(subplot_locs) %>%
  mutate(subplot_sp_relationship = PC_val*subplot_loc) %>%
  mutate(rel_cat = case_when(subplot_sp_relationship < -0.05 ~ 'strong negative',
                             subplot_sp_relationship > -0.05 & subplot_sp_relationship < 0.05 ~ 'neutral',
                             subplot_sp_relationship > 0.05 ~ 'strong positive'))
#how each covertype loads in the first 3 components of PCA:
loads_locs %>%
  ggplot() +
  geom_point(aes(x = covtype, y = PC_val, colour = PC)) +
  geom_hline(yintercept = 0) +
  theme_classic()

# shows us that PC1 is strongly affected by... 1) soil (positively), 2) standing dead (negatively), plants (negatively), and 3) rock (positively)
# PC2 is strongly affected by litter, wood (positively) and negatively with lichen and moss
# PC3 affected by fungi (negatively), plants (positively)

# + PC1 = hella soil, some rock; little standing dead or plants
# - PC1 = lots of standing dead or plants

# + PC2 = hella litter and wood, little lichen and moss
# - PC2 = hella lichen and moss, little litter and wood

# + PC3 = plants
# - PC3 = fungi

#relationship between each covertype and plot
loads_locs %>%
  ggplot() +
  geom_point(aes(x = subplot, y = subplot_sp_relationship, colour = covtype)) +
  #geom_label(aes(label = covtype)) +
  theme_classic()

# Fractional Cover: Soil --------------------------------------------------
#characterizing dominant cover:
soil_dom <- cover_1m %>%
  mutate(coverTypeRecode = ifelse(coverType == 'soil', 'soil', 'other')) %>%
  group_by(eventID, coverTypeRecode) %>%
  summarise(totcov = sum(percentCover))

soil_dom_summ <- soil_dom %>%
  group_by(eventID) %>%
  summarise(total = sum(totcov)) %>%
  left_join(soil_dom) %>%
  mutate(rel_cov = 100*totcov/total) %>%
  select(eventID, coverTypeRecode, rel_cov) %>%
  filter(coverTypeRecode == 'soil') %>%
  arrange(desc(rel_cov)) %>%
  print(n = 200)

#use their soil number:
actual_soil_cover <- cover_1m %>%
  filter(coverType == 'soil') %>%
  distinct(eventID, percentCover)
actual_soil_cover

# Visualize and understand the distributions 

# Summarize: 
# Whole plot, subplot size-specific analyses of % cover
# relative frequency for plot and larger subplots

# among (matrix rows) in terms of cover and composition, 
#which is then used to characterize each organizational level: 
 # [a diagram of subplot, plot, site, domain, all NEON sites]
#and what we're dealing w here...

#a single site for now

#output:
#loadings of each species in each PC; score of each component in each subplot /. plot

#show clustering by species groups hopefully

#consider how this might translate to a map, spatial / grid patterns added to what we already see
