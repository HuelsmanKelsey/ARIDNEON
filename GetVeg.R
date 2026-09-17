
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
# Load vegetation dataset and choose a site to filter:
repo_dir <- '/Users/khuelsma/Desktop/ARIDNEON/'
setwd(repo_dir)
#plots and subplots from ARID sites: CPER, RMNP, STER, JORN, SRER
ARID_veg <- read.csv('NEON_ARID_veg.csv')
subplots <- read.csv('NEON_ARID_subs.csv')

#subplots from which_site only
which_site <- 'CPER'
site_veg_subs <- subplots %>%
  filter(siteID == which_site)
site_veg <- ARID_veg %>%
  filter(siteID == which_site)

# Crosswalk / interpret the data
# Creating a single Event dataframe ---------------------------------------
Event_subplots <- site_veg_subs %>%
  dplyr::select(domainID, 
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
Event <- site_veg %>%
  #create variables that are in the subplot dataset so we can combine the two
  mutate(divDataType = NA,
         otherVariablesPresent = NA,
         otherVariables = NA,
         percentCover = NA) %>%
  dplyr::select(names(Event_subplots)) %>%
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
      sampleSizeValue == 100 ~ paste0(plotID, '_', subplotID, '_', year, '_', boutNumber),
      sampleSizeValue == 10 ~ paste0(plotID, '_', as.character(strtrim(subplotID, 5)), '0', '_', year, , '_', boutNumber),
      sampleSizeValue == 1 ~ paste0(plotID, '_', strtrim(subplotID, 4), '00', '_', year, '_', boutNumber)),
    contained_within_10 = case_when(
      sampleSizeValue == 1 ~ paste0(plotID, '_', str_replace(subplotID, '_1_', '_10_'), '_', year, '_', boutNumber),
      sampleSizeValue == 10 ~ paste0(plotID, '_', subplotID, '_', year, '_', boutNumber))
  ) %>% 
  dplyr::select(
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
Occurrence_subplots <- site_veg_subs %>%
  dplyr::select(uid, endDate, plotID, subplotID, boutNumber, release, publicationDate, 
         samplingImpractical, samplingImpracticalRemarks, samplingProtocolVersion,
         otherVariablesPresent, divDataType, otherVariables, percentCover,  #divDataType is plantSpecies, otherVariables
         heightPlantOver300cm, morphospeciesID, morphospeciesIDRemarks, 
         taxonID, taxonIDRemarks, remarks, nativeStatusCode, 
         taxonRank, family, scientificName, targetTaxaPresent, heightPlantOver300cm,
         heightPlantSpecies, recordedBy, measuredBy, identificationHistoryID, identificationReferences, identificationQualifier)
# ^ this gets rbinded in the next section
Occurrence <- site_veg %>%
  mutate(divDataType = NA,
         otherVariablesPresent = NA,
         otherVariables = NA,
         percentCover = NA,
         heightPlantOver300cm = NA, 
         heightPlantSpecies = NA) %>%
  dplyr::select(names(Occurrence_subplots)) %>%
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
 # filter(!is.na(coverType)) %>% 
  left_join(Event, 
            by = c('eventID', 'eventDate', 'year'),
            relationship = "many-to-many") %>%
  distinct()

#ARID_Obs is our input!

# Make PA dataframe:
#if which_size is null it will do whole plot; if which size isn't, it will do the size provided (1, 10, or 100)
make_PA_df <- function(which_site, 
                       which_year, 
                       which_species, 
                       which_size = NULL,
                       matrix = FALSE) #default
  {
  
  if (is.null(which_size)) { #default for whole plot
    wide_df <- ARID_Obs %>%
      filter(year == which_year) %>%
      distinct(plotID, boutNumber, year, scientificName) %>%
      mutate(recode = case_when(scientificName == which_species ~ 'SOI',
                                scientificName != which_species ~ 'other'),
             eventID = paste0(plotID, '_', year, '_', boutNumber))
  }

  if (!is.null(which_size)) { 
    if (which_size == 1) {
      ARID_Obs <- ARID_Obs %>%
        filter(sampleSizeValue == which_size)
    }
    
    wide_df <- ARID_Obs %>%
      filter(year == which_year) %>%
      distinct(plotID, subplotID, sampleSizeValue, contained_within_10, contained_within_100, boutNumber, year, scientificName) %>%
      mutate(recode = case_when(scientificName == which_species ~ 'SOI',
                                scientificName != which_species ~ 'other'),
             eventID = case_when(which_size == 1 ~ paste0(plotID, '_', subplotID, '_', year, '_', boutNumber), #if 1m use subplot ID; 
                                 which_size == 10 ~ paste0(contained_within_10),
                                 which_size == 100 ~ paste0(contained_within_100)))
  }
  
  wide_df <- wide_df %>%
    group_by(eventID) %>%
    count(recode) %>%
    na.omit() %>%
    pivot_wider(names_from = recode, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(eventID, SOI, other)
  
  #if want a matrix:
  if (matrix == TRUE) {
    PA_mat <- wide_df %>%
      dplyr::select(SOI, other) %>%
      as.matrix()
    rownames(PA_mat) <- wide_df$eventID
    
    return(PA_mat)
  }
  return(wide_df)
}


make_PA_df('CPER', 2021, 'soil', which_size = 1, matrix = FALSE)

#Can use make_PA_df to also characterize plot or subplot, but will need to make it so there is a which_species == NULL option
#because we want to characterize by everything that is there, not just whether one species is present or absent

# PA_sp_characterization <- function(which_site, which_year, which_species, which_size, matrix) {
#   which_plant_deets <- make_PA_df(which_site, 
#                                   which_year, 
#                                   which_species = NULL, 
#                                   which_size, 
#                                   matrix = TRUE)
#   vegdeets_pca <- vegan::pca(which_plant_deets, scale = TRUE)
#   subplot_locs <- vegdeets_pca$CA$u[,1:3] %>%
#     as.data.frame() %>%
#     mutate(subplot = rownames(vegdeets_pca$CA$u)) %>%
#     pivot_longer(PC1:PC3, names_to = 'PC', values_to = 'subplot_loc')
#   #each species
#   plant_deets_load <- vegdeets_pca$CA$v[,1:3] %>% #22 x 22
#     as.data.frame() %>%
#     mutate(covtype = rownames(vegdeets_pca$CA$v)) %>%
#     pivot_longer(PC1:PC3, names_to = 'PC', values_to = 'PC_val')
#   loads_locs <- plant_deets_load %>%
#     left_join(subplot_locs, relationship = 'many-to-many') %>%
#     mutate(subplot_veg_relationship = PC_val*subplot_loc) %>%
#     mutate(rel_cat = case_when(subplot_veg_relationship < -0.05 ~ 'strong negative',
#                                subplot_veg_relationship > -0.05 & subplot_veg_relationship < 0.05 ~ 'neutral',
#                                subplot_veg_relationship > 0.05 ~ 'strong positive'))
#   #how each covertype loads in the first 3 components of PCA:
#   loads_locs %>%
#     filter(PC_val > 0.2 | PC_val < -0.2) %>%
#     ggplot() +
#     geom_point(aes(x = covtype, y = PC_val, colour = PC)) +
#     geom_hline(yintercept = 0) +
#     theme_classic() +
#     facet_wrap(~PC, scales = 'free_x') +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))
#   
#   #relationship between each covertype and plot
#   loads_locs %>%
#     ggplot() +
#     geom_point(aes(x = subplot, y = subplot_veg_relationship, colour = covtype)) +
#     #geom_label(aes(label = covtype)) +
#     theme_classic()
#   
#   return(loads_locs)
# }

  
# explore and make % cover dataframe:

# summary of % cover at each coverlocation (ground, understory, overstory)
tot_cov_by_group <- ARID_Obs %>% 
  filter(year == which_year) %>%
  filter(sampleSizeValue == 1) %>%
  group_by(subplotID, plotID, year, boutNumber, coverLocation) %>% 
  summarise(totcov = sum(percentCover)) 

#caveat #1: how % cover is reported:
tot_cov_by_group

tot_cov_by_subplot <- tot_cov_by_group %>%
  group_by(subplotID, plotID, boutNumber, year) %>%
  summarise(totcov2 = sum(totcov))

tot_cov_by_subplot

#percent cover of species or cover of interest
make_perccov_df <- function(which_site, 
                            which_year, 
                            which_cover,
                            which_size = NULL) {
  
  #this is every single 1m subplot from CPER 2024
  cover_1m <- ARID_Obs %>% #or can do ARID_Obs for multiple bouts!
    filter(sampleSizeValue == 1) %>%
    filter(locationID == which_site) %>%
    filter(year == which_year) %>%
    distinct(eventID, plotID, subplotID, boutNumber, year,
             coverType, coverTypeGeneral, percentCover, 
             contained_within_10, contained_within_100)
  
  if (is.null(which_size)) { #default for whole plot
    cover_df <- cover_1m %>%
      mutate(recode = case_when(coverType == which_cover ~ 'SOI',
                                coverType != which_cover ~ 'other'),
             eventID = paste0(plotID, '_', year, '_', boutNumber))
  }
  
  if (!is.null(which_size)) { 
    cover_df <- cover_1m %>%
        mutate(recode = case_when(coverType == which_cover ~ 'SOI',
                                  coverType != which_cover ~ 'other'),
             
               eventID = case_when(which_size == 1 ~ paste0(plotID, '_', subplotID, '_', year, '_', boutNumber), #if 1m use subplot ID; 
                                 which_size == 10 ~ paste0(contained_within_10),
                                 which_size == 100 ~ paste0(contained_within_100)))
  }
  
  wide_df <- cover_df %>%
    group_by(eventID, recode) %>%
    summarise(mean_cov = mean(percentCover, na.rm = TRUE)) %>%
    pivot_wider(names_from = recode, values_from = mean_cov, values_fill = 0) %>%
    dplyr::select(eventID, SOI) %>%
    mutate(which_cover = paste(which_cover)) %>%
    print()
}


#3 options to characterize the vegetation at the plot using 1m percent cover data:
#1) group all plants (coverTypeGeneral) (general = TRUE, species_only = FALSE)
#2) use scientific name to separate all species but include % cover (general = FALSE, species_only = FALSE)
#3) percent cover of plants only (general = FALSE, species_only = TRUE)

#make veg matrix, then characterize veg, which returns "loading_locs"
#loadings of each species in each PC; score of each component in each subplot /plot

make_perccov_matrix <- function(which_site, which_year, general, species_only = FALSE #by default
                            ) {
  
  cover_1m <- ARID_Obs %>% #or can do ARID_Obs for multiple bouts!
    filter(sampleSizeValue == 1) %>%
    filter(locationID == which_site) %>%
    filter(year == which_year) %>%
    distinct(eventID, plotID, subplotID, boutNumber, year,
             coverType, coverTypeGeneral, percentCover, 
             contained_within_10, contained_within_100)
  
  if (general == TRUE) {
    group_summary <- cover_1m %>%
      group_by(eventID, coverTypeGeneral) %>%
      summarise(cov = sum(percentCover))
    
    group_summary_wide <- group_summary %>%
      distinct(eventID, coverTypeGeneral, cov) %>%
      group_by(eventID) %>%
      pivot_wider(names_from = coverTypeGeneral,
                  values_from = cov,
                  values_fill = 0) %>%
      ungroup() 
    mat_gen <- group_summary_wide %>%
      dplyr::select(-eventID) %>%
      as.matrix()
    row.names(mat_gen) <- group_summary_wide$eventID
    
    return(mat_gen)
  }
  
  if (general == FALSE) {
   
    group_summary <- cover_1m %>%
      group_by(eventID, coverType) %>%
      summarise(cov = sum(percentCover))
    
    if (species_only == TRUE) {
      group_summary <- cover_1m %>%
        filter(coverTypeGeneral == 'plants') %>%
        group_by(eventID, contained_within_10, contained_within_100, coverType) %>%
        summarise(cov = sum(percentCover))
    }
    group_summary_wide <- group_summary %>%
      ungroup() %>%
      distinct(eventID, coverType, cov) %>%
      pivot_wider(names_from = coverType,
                  values_from = cov,
                  values_fill = 0) %>%
      ungroup()
    
    mat_spec <- group_summary_wide %>%
      dplyr::select(-eventID) %>%
      as.matrix()
    row.names(mat_spec) <- group_summary_wide$eventID
  
    return(mat_spec)
    }
}


characterize_veg <- function(which_site, which_year, general, species_only) {
  which_plant_deets <- make_perccov_matrix(which_site, which_year, general, species_only)
  vegdeets_pca <- vegan::pca(which_plant_deets, scale = TRUE)
  subplot_locs <- vegdeets_pca$CA$u[,1:3] %>%
    as.data.frame() %>%
    mutate(subplot = rownames(vegdeets_pca$CA$u)) %>%
    pivot_longer(PC1:PC3, names_to = 'PC', values_to = 'subplot_loc')
  #each species
  plant_deets_load <- vegdeets_pca$CA$v[,1:3] %>% #22 x 22
    as.data.frame() %>%
    mutate(covtype = rownames(vegdeets_pca$CA$v)) %>%
    pivot_longer(PC1:PC3, names_to = 'PC', values_to = 'PC_val')
  loads_locs <- plant_deets_load %>%
    left_join(subplot_locs, relationship = 'many-to-many') %>%
    mutate(subplot_veg_relationship = PC_val*subplot_loc) %>%
    mutate(rel_cat = case_when(subplot_veg_relationship < -0.05 ~ 'strong negative',
                               subplot_veg_relationship > -0.05 & subplot_veg_relationship < 0.05 ~ 'neutral',
                               subplot_veg_relationship > 0.05 ~ 'strong positive'))
  #how each covertype loads in the first 3 components of PCA:
  loads_locs %>%
    filter(PC_val > 0.2 | PC_val < -0.2) %>%
    ggplot() +
    geom_point(aes(x = covtype, y = PC_val, colour = PC)) +
    geom_hline(yintercept = 0) +
    theme_classic() +
    facet_wrap(~PC, scales = 'free_x') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  #relationship between each covertype and plot
  loads_locs %>%
    ggplot() +
    geom_point(aes(x = subplot, y = subplot_veg_relationship, colour = covtype)) +
    #geom_label(aes(label = covtype)) +
    theme_classic()
  
  return(loads_locs)
}

veg_char_2021 <- characterize_veg('CPER', 2021, general = TRUE, species_only = FALSE)
veg_char_2021 %>% group_by(PC, subplot, covtype, rel_cat) %>% summarise(n = n()) %>%
  filter(rel_cat != 'neutral') %>%
  group_by(PC, subplot, rel_cat) %>%
  summarise(covs = list(unique(covtype))) %>%
  ungroup() %>%
  mutate(rel_cat_simplified = ifelse(rel_cat == 'strong negative', 'neg', 'pos'),
         summary = paste(rel_cat_simplified, covs)) %>%
  select(PC, subplot, summary)

#plot: CPER 1
#3 subplots: 40 1 1, 40 1 3, 41 4 1, -litter, plants, + standing dead; generally the vibe across the plot

#CPER3, 1 subplot: neg standing dead, positive plants. 

#CPER 5: positive litter, plants, wood

veg_char_2024 <- characterize_veg('CPER', 2024, general = TRUE, species_only = FALSE)


#Looking at bouts within a year:
ARID_Obs %>%
  group_by(year, siteNumber) %>%
  summarise(n_bouts = n_distinct(boutNumber)) %>%
  filter(n_bouts > 1) %>% 
  distinct(siteNumber, year, n_bouts)

#only one year with two veg bouts??

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
  dplyr::select(-mvals, -perccovs, -mean_val, -max_val, -max_bout, -mean_bout, -n_bouts, -n_remarks) %>%
  print()