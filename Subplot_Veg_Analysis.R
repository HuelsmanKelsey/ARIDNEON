library(mixOmics)
library(pls) 

#how to combine vegetation details and spectra:

# workflow: ---------------------------------------------------------------
which_site <- 'CPER' #4 letter NEON code
which_rast <- 'hsi' #or rgb
which_size <- 1

# workflow: 2020 ----------------------------------------------------------
#TBD because of multiple bouts

# workflow: 2021 ----------------------------------------------------------
which_year <- 2021 
CPER_spectra_2021 <- extract_subplot_spectra(which_rast, which_site, which_year, which_size)
CPER_gen_cover_2021 <- make_perccov_matrix(which_site, which_year, general = TRUE, species_only = FALSE)
CPER_spec_cover_2021 <- make_perccov_matrix(which_site, which_year, general = FALSE, species_only = FALSE)
CPER_plants_only_2021 <- make_perccov_matrix(which_site, which_year, general = FALSE, species_only = TRUE)
# workflow: 2024 ----------------------------------------------------------
which_year <- 2024 
CPER_spectra_2024 <- extract_subplot_spectra(which_rast, which_site, which_year, which_size)
CPER_gen_cover_2024 <- make_perccov_matrix(which_site, which_year, general = TRUE, species_only = FALSE)
CPER_spec_cover_2024 <- make_perccov_matrix(which_site, which_year, general = FALSE, species_only = FALSE)
CPER_plants_only_2024 <- make_perccov_matrix(which_site, which_year, general = FALSE, species_only = TRUE)


#choose a SPECIFIC cover type, choose a size to focus on subplots or whole plot if a size isn't given.

#do for all the covers, all the years, all the buffers
covers <- c('Opuntia polyacantha', 
            'Bouteloua gracilis',
            'Bouteloua dactyloides',
            'soil')
which_buff <- c(0, 2)

get_perccov_spectra <- function(which_site, which_year, which_cover, which_buff) {
  extracted <- extract_subplot_spectra('hsi', 
                                       which_site, 
                                       which_year,
                                       1,  #1 = which_size in extract subplot spectra
                                       which_buff) 
  perccov <- make_perccov_df(which_site, 
                             which_year, 
                             which_cover,
                             1)
  extracted_perccov <- extracted %>%
    left_join(perccov)
  kept_bands <- get_wvs(which_site, which_year)
  input <- extracted_perccov %>%
    rename(perc_cov = SOI,
           cover_type = which_cover) %>%
    pivot_longer(cols = any_of(unique(kept_bands$band)),
                 names_to = 'band', values_to = 'refl') %>%
    left_join(kept_bands) %>%
    group_by(eventID, plotID, subplotID, nlcdClass, easting, northing, latitude, longitude, horzUncert,
             subpltSize, refPoint, band, year, wv, cover_type, perc_cov) %>%
    summarise(mean_refl = mean(refl),
              sd_refl = sd(refl),
              se_refl = sd(refl)/sqrt(n()),
              CV_refl = sd_refl/mean_refl) %>%
    select(eventID, plotID, subplotID, nlcdClass, band, wv, year, cover_type, perc_cov,
           mean_refl, se_refl, CV_refl)
}
#if we are extracting rgb it's really just to see what subplot heterogeneity looks like...

perccov_soil_2021_buffed <- get_perccov_spectra('CPER', 2021, 'soil', which_buff = 2)
perccov_soil_2021 <- get_perccov_spectra('CPER', 2021, 'soil', which_buff = 0)
perccov_soil_2024_buffed <- get_perccov_spectra('CPER', 2024, 'soil', which_buff = 2)
perccov_soil_2024 <- get_perccov_spectra('CPER', 2024, 'soil', which_buff = 0)



perccov_soil_2021 %>%
  filter(!is.na(cov_cat)) %>%
  group_by(cov_cat, wv) %>%
  summarise(mean_refl_by_cover = mean(mean_refl),
            se_refl_by_cover = sd(mean_refl)/sqrt(n())) %>%
  ggplot(aes(x = wv, y = mean_refl_by_cover)) +
  geom_point(aes(colour = cov_cat)) +
  geom_errorbar(aes(x = wv, 
                    ymin = mean_refl_by_cover + 2*se_refl_by_cover,
                    ymax = mean_refl_by_cover + 2*se_refl_by_cover,
                    colour = cov_cat))+
  theme_classic()



# Each Year ---------------------------------------------------------------

input <- perccov_soil_2021
bands <- paste0('X', unique(input$wv))
wvs <- unique(input$wv)

input <- perccov_soil_2021 %>%
  ungroup() %>%
  select(eventID, plotID, subplotID, year, 
         nlcdClass, 
         cover_type, perc_cov, wv, mean_refl) %>%
  pivot_wider(names_from = wv, names_prefix = 'X',
              values_from = mean_refl) %>%
  filter(!is.na(cover_type))

set.seed(123)
#sample 70% of the data; choosing what rows to sample for consistency
sampled_indices <- sample(nrow(input), size = .7*nrow(input))

train <- input[sampled_indices, ]
test <- input[-sampled_indices, ]

spectra_train <- train %>%
  dplyr::select(all_of(bands))

spectra_test <- test %>%
  dplyr::select(all_of(bands))

info_train <- train %>%
  dplyr::select(-names(spectra_train))

info_test <- test %>%
  dplyr::select(-names(spectra_test))

train_df <- data.frame(info_train, 
                       Spectra = I(as.matrix(spectra_train)))

test_df <- data.frame(info_test, 
                      Spectra = I(as.matrix(spectra_test)))
which_var <- 'perc_cov'
#pls::plsr(IV, DV, ncomp = nComps, scale = TRUE)
plsr.out <- pls::plsr(
  as.formula(paste(which_var,"~",'Spectra')),
  scale = TRUE, 
  ncomp = nComps, 
  validation = "LOO", #leave one out
  method = "oscorespls", #kernelpls
  segment.type = "interleaved", 
  seg = 10,
  trace = FALSE, 
  center = TRUE, #default
  data = train_df) #cal data

plsr.out$model
plsr.out$coefficients
plsr.out$scores
plsr.out$fitted.values
plsr.out$residuals


#get R2 and RMSE(P)
R2 = round(
  pls::R2(plsr.out, intercept = F) [[1]] [nComps], 
  2)
RMSEP = pls::RMSEP(
  plsr.out,
  estimate = c("test"),
  newdata = test_df)
cal_df <- train_df
val_df <- test_df

#need to compare each original dataset to model:
cal_df <- cal_df %>%
  dplyr::mutate(
    dataset = 'cal',
    pred = as.vector(
      plsr.out$validation$pred[,,nComps]),
  )
cal_df

val_df <- val_df %>%
  dplyr::mutate(
    dataset = 'val',
    pred = as.vector(
      predict(
        plsr.out, 
        newdata = val_df, 
        ncomp = nComps, 
        type = 'response')
    ))


cal_val_df <- rbind(cal_df, val_df) %>%
  rename(var = perc_cov) %>%
  dplyr::mutate(resid = pred - var,
                r2 = case_when(
                  dataset == 'cal' ~ R2, #calc above?
                  dataset == 'val' ~ resid^2),
                min = quantile(x = var, probs = 0.001),
                max = quantile(x = var, probs = 0.999))
scatterplot

scatterplot <- cal_val_df %>%
  ggplot(aes(x = pred, y = var)) + 
  theme_bw() + 
  geom_point() +
  geom_abline(intercept = 0,
              slope = 1,
              color = "dark grey",
              linetype = "dashed",
              linewidth = 1.5) +
  #xlim using min and max??
  labs(x = paste0("Predicted ", paste(which_var), " (units)"),
       y = paste0("Observed ", paste(which_var), " (units)")) +
  theme(axis.text = element_text(size = 18), 
        legend.position = "none",
        axis.title=element_text(size=20, face="bold"),
        axis.text.x = element_text(angle = 0,vjust = 0.5)) +
  facet_wrap(~dataset)

resid_histogram <- cal_val_df %>%
  ggplot(aes(x = resid)) +
  geom_histogram(alpha = .5, position = "identity") + 
  geom_vline(xintercept = 0, color="black", 
             linetype = "dashed", linewidth = 1) + 
  theme_bw() + 
  theme(axis.text = element_text(size=18), 
        legend.position = "none", #no legend
        axis.title = element_text(size=20, face="bold"), 
        axis.text.x = element_text(angle = 0, vjust = 0.5)) +
  facet_wrap(~dataset)


this_pls <- mixOmics::pls(IV, DV, ncomp = 2, scale = TRUE, mode = "regression")

this_pls$prop_expl_var

#how X (IV) and Y (DV) are positioned in the components:
IV_position <- this_pls$variates$X
DV_position <- this_pls$variates$Y

#how the wavelengths load
wv_loads <- as.data.frame(this_pls$loadings$X) %>%
  mutate(wv = wvs)
wv_loads %>%
  ggplot(aes(x = wv, y = comp1)) +
  geom_point() +
  theme_classic()

predict(test_df, this_pls)

this_pls

test_df <- data.frame(info_test, 
                     Spectra = I(as.matrix(spectra_test))) %>%
  print()

perccov_soil_2024


# Both years --------------------------------------------------------------

all_soil <- soil_2024 %>%
  rbind(soil_2021) %>%
  filter(mean_cov >= 1) %>%
  pivot_wider(names_from = year, names_prefix = 'Y', values_from = mean_cov) %>%
  mutate(change = 100*(Y2024 - Y2021)/(Y2021+0.01)) %>%
  filter(!is.na(change))

all_soil %>%
  arrange(change) %>%
  filter(!is.na(change)) %>%
  print(n = 600)

all_soil %>%
  group_by(plotID) %>%
  summarise(plot_mean_change = mean(change, na.rm = TRUE),
            plot_sd_change = sd(change, na.rm = TRUE),
            plot_se_change = plot_sd_change/sqrt(n())) %>%
  ggplot(aes(x = plotID, y = plot_mean_change)) +
  geom_point() +
  geom_errorbar(aes(x = plotID, 
                    ymin = plot_mean_change - plot_se_change,
                    ymax = plot_mean_change + plot_se_change))
#let's look at what happened to plot 25, 16, 14, 20, 22, 2, and 28... 


# 2) get_PA_spectra:
# req raster inputs; choose a cover type, choose a size to focus on subplots or whole plot if a size isn't given.

#presence / absence in each subplot or whole plot
get_PA_spectra <- function(which_site, which_year, which_rast, which_size, which_species) {
  #use extract_subplot_spectra with hsilist + which_size;
  extracted <- extract_subplot_spectra(which_rast, which_site, which_year, which_size)
  #make a P/A dataframe for the species of interest + which_size for every subplot
  PA <- make_PA_df(which_site, which_year, which_species, which_size, matrix = FALSE)
  
  extracted_PA <- extracted %>%
    left_join(PA) %>%
    na.omit() %>%
    distinct() %>%
    print()
}



which_size = 1

ggplot(soil_cov_1) +
  #geom_density(aes(mean_cov), lwd = 2) +
  geom_histogram(aes(), bins = 20, alpha = 0.5) +
  theme_classic()
#shows us that most 1m2 plots have 10-15% soil cover (nearly 1000?)
#second most: 5-10, then 20-25


#train on the 8 x 1m2 subplots in each plot: 
plot_summary_1m <- soil_spec_1m %>%
  group_by(plotID, wv) %>%
  summarise(n = n_distinct(subplotID))

  
#can add in details from RGB imagery: 


#generate a field-based map of something:

#generate a map from RGB imagery of site or subplot:





  # Use permutation to determine the optimal number of components
  if(grepl("Windows", sessionInfo()$running)){
    pls.options(parallel = NULL)
  } else {
    pls.options(parallel = parallel::detectCores()-1)
  }
  
  seg <- 10
  maxComps <- 10
  iterations <- 10 #gets used later; was 1000 previously but keeping it simple
  inVar <- 'var' 
  
  nComps <- spectratrait::find_optimal_components(
    dataset = cal_df,
    targetVariable = 'var',
    method = 'pls',
    seg = seg,
    maxComps = maxComps, 
    prop = 0.70,
    random_seed = 17)
  
  print(paste0("*** Optimal number of components: ", nComps))
  
  #preliminary assessment
  plsr.out <- plsr(
    as.formula(paste(inVar,"~","Spectra")),
    scale = TRUE, 
    ncomp = nComps, #calc above
    validation = "LOO", #leave one out
    method = "oscorespls", #kernelpls
    segment.type = "interleaved", 
    seg = seg,
    trace = FALSE, 
    center = TRUE, #default
    data = cal_df) #cal data
  
  #get R2 and RMSE(P)
  R2 = round(
    pls::R2(plsr.out, intercept = F) [[1]] [nComps], 
    2)
  RMSEP = pls::RMSEP(
    plsr.out,
    estimate = c("test"),
    newdata = val_df)
  
  #need to compare each original dataset to model:
  cal_df <- cal_df %>%
    dplyr::mutate(
      dataset = 'cal',
      pred = as.vector(
        plsr.out$validation$pred[,,nComps]),
    )
  
  val_df <- val_df %>%
    dplyr::mutate(
      dataset = 'val',
      pred = as.vector(
        predict(
          plsr.out, 
          newdata = val_df, 
          ncomp = nComps, 
          type = 'response')
      ))
  
  cal_val_df <- rbind(cal_df, val_df) %>%
    dplyr::mutate(resid = pred - var,
                  r2 = case_when(
                    dataset == 'cal' ~ R2, #calc above?
                    dataset == 'val' ~ resid^2),
                  min = quantile(x = var, probs = 0.001),
                  max = quantile(x = var, probs = 0.999))
  
  scatterplot <- cal_val_df %>%
    ggplot(aes(x = pred, y = var)) + 
    theme_bw() + 
    geom_point() +
    geom_abline(intercept = 0,
                slope = 1,
                color = "dark grey",
                linetype = "dashed",
                linewidth = 1.5) +
    #xlim using min and max??
    labs(x = paste0("Predicted ", paste(which_var), " (units)"),
         y = paste0("Observed ", paste(which_var), " (units)")) +
    theme(axis.text = element_text(size = 18), 
          legend.position = "none",
          axis.title=element_text(size=20, face="bold"),
          axis.text.x = element_text(angle = 0,vjust = 0.5),
          panel.border = element_rect(linetype = "solid", fill = NA, size=1.5)) +
    facet_wrap(~dataset)
  
  resid_histogram <- cal_val_df %>%
    ggplot(aes(x = resid)) +
    geom_histogram(alpha = .5, position = "identity") + 
    geom_vline(xintercept = 0, color="black", 
               linetype = "dashed", size = 1) + 
    theme_bw() + 
    theme(axis.text = element_text(size=18), 
          legend.position = "none", #no legend
          axis.title = element_text(size=20, face="bold"), 
          axis.text.x = element_text(angle = 0, vjust = 0.5),
          panel.border = element_rect(linetype = "solid", fill = NA, size=1.5)) +
    facet_wrap(~dataset)
  
  norm_check <- grid.arrange(scatterplot,  
                             resid_histogram,
                             nrow = 2, 
                             ncol = 1)
  print(norm_check)
}

norm_check_plsr('all', 'all', 'N_area', buck_only = TRUE)
norm_check_plsr('all', 'all', 'N_area')

norm_check_plsr('all', 'all', 'CN', buck_only = TRUE)
norm_check_plsr('all', 'all', 'LMA', buck_only = TRUE)

norm_check_plsr('all', 'all', 'N_area')
norm_check_plsr('all', 'all', 'CN')
norm_check_plsr('all', 'all', 'LMA')


spectratrait::f.coef.valid()
spectratrait::get_ecosis_data()
spectratrait::percent_rmse()
spectratrait::VIP() #object; where object is the fitted plsr object
#spectratrait::VIPjh() #variable j with h components

#plsr.out = model from jackknife procedure,

#data_plsr = spectra of matrix for model,

#nComp and inVar
norm_check_plsr(cal_data = 'all', val_data = 'all', which_var = 'LMA')


#final assessment: jackknife data 
inVar <- 'LMA'
nComps <- 2
segments = 1000

jk.plsr.out <- pls::plsr(
  as.formula(paste(inVar,"~","Spectra")),
  scale = FALSE, 
  center = TRUE, 
  ncomp = nComps, 
  validation = "CV", 
  method = "oscorespls", #kernelpls
  segments = 100, 
  segment.type = "interleaved",
  trace = FALSE, 
  jackknife = TRUE, #this does it all automatically
  data = cal_df)

Jackknife_coef <- spectratrait::f.coef.valid(
  plsr.out = jk.plsr.out, 
  data_plsr = cal_df, 
  ncomp = nComps, 
  inVar = inVar)

Jackknife_intercept <- Jackknife_coef[1,,,]
Jackknife_coef <- Jackknife_coef[2:dim(Jackknife_coef)[1],,,]
interval <- c(0.025,0.975)

Jackknife_intercept <- matrix( 
  rep(Jackknife_intercept, length(val_df$var), by_row = TRUE),
  ncol = length(Jackknife_intercept))

Jackknife_Pred <- val_df$Spectra %*% Jackknife_coef + Jackknife_intercept

Interval_Conf <- apply(X = Jackknife_Pred, 
                       MARGIN = 1, 
                       FUN = quantile, 
                       probs = c(interval[1], interval[2])) #gives us two values for 95% CI 
#Interval_Conf[1,] and [2,] are LCI and UCI, respectively

sd_mean <- apply(X = Jackknife_Pred, 
                 MARGIN = 1, 
                 FUN = sd)
#standard deviation of the residuals
sd_res <- sd(val_df$resid) #r2 is other value
#total standard deviation
sd_tot <- sqrt(sd_mean^2+sd_res^2) #total standard deviation

Jackknife_Pred #spectra*coeffs + intercept
Interval_Conf #95% CI

#vips
vips <- spectratrait::VIP(jk.plsr.out)

#predicted values
#LCI -- confidence intervals
#LPI - lowest predicted interval (predicted value - 1.96*st. error of mean)

#possible response variables: LMA, N_area, C_area, CN

spectratrait::f.plot.coef(
  Z = t(Jackknife_coef), wv = wvs, 
  plot_label="Jackknife regression coefficients",
  position = 'topleft')
abline(h=0,lty=2,col="grey50")
box(lwd=2.2)

# JK validation plot
rmsep_percrmsep <- spectratrait::percent_rmse(
  plsr_dataset = val.plsr.output, 
  inVar = inVar, 
  residuals = val.plsr.output$PLSR_Residuals, 
  range="full")

RMSEP <- rmsep_percrmsep$rmse

perc_RMSEP <- rmsep_percrmsep$perc_rmse

r2 <- round(pls::R2(plsr.out, newdata = val.plsr.data, intercept=F)$val[nComps],2)
}

