library(tidyverse)
library(sf)
library(terra)
library(mapview)
options(mapview.maxpixels = 9000000)
library(ggrepel)
theme_set(theme_minimal())

# TODO: GET RASTER DATED FOR 2020!

## Helper functions

# Load a raster, cropped, projected, and masked to the study area
load_raster = function(path_rast) {
  r = rast(path_rast)
  r_crop = crop(r, vect(st_transform(study_area, crs(r))))
  r_proj = project(r_crop, crs(study_area))
  r_mask = mask(r_proj, vect(study_area))
  return(r_mask)
}

# Computer raster function value for a buffer
computer_raster_buffer_value_func = function(raster, points, buffer_distance, func) {
  r = project(raster, crs(points))
  site_buffers = st_buffer(points, dist = buffer_distance)
  buffer_values = terra::extract(r, vect(site_buffers), fun = func, na.rm = TRUE)
}

### Study area

crs_m = 32610 # EPSG:32610, UTM Zone 10N, (meters)

# ARU locations
aru_sites = st_read('data/environment/GIS Data/AcousticStations.shp') %>%
  st_drop_geometry() %>% janitor::clean_names() %>% as.data.frame() %>%
  st_as_sf(coords = c("utm_e", "utm_n"), crs = crs_m) %>%
  select(name, ces, treatment, geometry) %>% rename(site = name)
mapview(aru_sites, label = aru_sites$name)

site_data = aru_sites

# site_data = gpx::read_gpx('data/sites/PAMlocations20230510.gpx')$waypoints %>%
#   janitor::clean_names() %>% rename(site = name) %>% select(site, latitude, longitude)
# site_data = st_transform(st_as_sf(site_data, coords = c('longitude', 'latitude'), crs = crs_data), crs = crs_projected)

# Study area bounding buffer
study_area = st_buffer(st_as_sfc(st_bbox(aru_sites)), dist = 100)

# Roads
roads = st_read('data/environment/GIS Data/T3roads.shp') %>%
  st_zm(drop = TRUE) %>% janitor::clean_names() %>% st_transform(crs = crs_m)

# Strata and forest inventory units

  # Seral stage class
  seral_classes = c("initiation", "canclose", "stemex", "mature")
  class_colors = viridis::viridis(4, option = 'D', direction = -1)
  
  poly_seral_class = st_read('data/environment/GIS Data/Forest Development Strata/AgeStrataFixedDIS_RSFRIS20200130.shp') %>%
    st_transform(crs_m) %>% janitor::clean_names() %>% mutate(class_seral = stratum %>% str_to_lower()) %>%
    mutate(stratum = factor(stratum, levels = seral_classes)) %>% select(-stratum)
  # rsfris_units = st_set_crs(st_read('data/environment/GIS Data/Forest Development Strata/AgeStrataFixedRSFRIS20200130.shp'), 'NAD83')
  
  # Thinning treatment
  poly_thinning_treatment = st_read('data/environment/GIS Data/Forest Development Strata/ThinAfter94NoHarvSinceClipByInitBuf3.shp') %>% 
    st_transform(crs_m) %>% select(TECHNIQUE_, FMA_DT, FMA_STATUS) %>% janitor::clean_names() %>%
    mutate(technique = technique %>% str_to_lower(), fma_status = fma_status %>% str_to_lower()) %>%
    rename(thinning_treatment = technique, thinning_status = fma_status, thinning_date = fma_dt)

# Classify sites based on the seral class and thinning prescription polygons in which they reside
site_data = st_join(site_data, poly_seral_class)
table(site_data$class_seral, useNA = 'ifany')
site_data = st_join(site_data, poly_thinning_treatment)
table(site_data$thinning_treatment, useNA = 'ifany')

mapview(study_area, col.regions = "transparent", color = "black", lwd = 2) +
  mapview(poly_seral_class, col.regions = class_colors) +
  mapview(site_data, zcol = 'class_seral', col.regions = class_colors, legend = T)

# Forest origin (age)
rast_origin = load_raster('data/environment/rs_fris/rs_fris_Origin_Year/RS_FRIS_ORIGIN_YEAR.tif')
sort(unique(terra::values(rast_origin)))

rast_age = 2024 - rast_origin

# Classify origin raster based on number of years since 2020
age_classes = c("early", "mid", "late", "old-growth")
rast_age_class = classify(rast_age, rcl = matrix(c(
  0,   24,  1,   # "early"
  25,  79,  2,   # "mid"
  80, 199,  3,   # "late"
  200, Inf, 4    # "old-growth"
), ncol = 3, byrow = TRUE))
levels(rast_age_class) = data.frame(value = 1:4, class_age = age_classes)
rast_age_class = ifel(rast_age_class %in% 1:4, rast_age_class, NA) # discard erroneous / NA values
unique(terra::values(rast_age_class))
levels(rast_age_class) = data.frame(value = 1:4, class_age = age_classes)

# Classify sites based on most abundant origin class within a 100m buffer
site_data$class_age = as.factor(computer_raster_buffer_value_func(rast_age_class, site_data, 100, modal)[,2])
table(as.numeric(site_data$class_age))

# Get mean age
site_data$age = as.numeric(computer_raster_buffer_value_func(rast_age, site_data, 100, modal)[,2])
table(site_data$age)

mapview(study_area, col.regions = "transparent", color = "black", lwd = 2) +
  mapview(as.factor(rast_age_class), col.regions = class_colors) +
  mapview(site_data, zcol = 'class_age', col.regions = class_colors, legend = T)

### Plot scale variables

# Elevation (m)

  # Local: TODO

  # Regional: TODO

# Total basal area (m2/ha)

  # Local: TODO

  # Regional
  rast_ba = load_raster('data/environment/rs_fris/rs_fris_BA/RS_FRIS_BA.tif')
  site_data$ba_mean_100m = as.numeric(computer_raster_buffer_value_func(rast_ba, site_data, 100, mean)[,2])
  site_data$ba_sd_100m = as.numeric(computer_raster_buffer_value_func(rast_ba, site_data, 100, sd)[,2])
  
  ggplot(site_data, aes(x = age, y = ba_mean_100m, label = site)) +
    geom_point() + geom_text_repel(size = 2)
  ggplot(site_data, aes(x = age, y = ba_sd_100m, label = site)) +
    geom_point() + geom_text_repel(size = 2)
  
  mapview(rast_ba) +
    mapview(site_data, zcol = 'ba_mean_100m', legend = T)

# Tree height distribution (e.g. SD of height) (m)
  
  # Regional
  rast_ht_lorey = load_raster('data/environment/rs_fris/rs_fris_HT_LOREY/RS_FRIS_HT_LOREY.tif')
  site_data$ht_lorey_mean_100m = as.numeric(computer_raster_buffer_value_func(rast_ht_lorey, site_data, 100, mean)[,2])
  site_data$ht_lorey_sd_100m = as.numeric(computer_raster_buffer_value_func(rast_ht_lorey, site_data, 100, sd)[,2])
  
  ggplot(site_data, aes(x = age, y = ht_lorey_mean_100m, label = site)) +
    geom_point() + geom_text_repel(size = 2)
  ggplot(site_data, aes(x = age, y = ht_lorey_sd_100m, label = site)) +
    geom_point() + geom_text_repel(size = 2)
  
  mapview(rast_ht_lorey) +
    mapview(site_data, zcol = 'ht_lorey_mean_100m', legend = T)
  
  rast_ht_max = load_raster('data/environment/rs_fris/rs_fris_HTMAX/RS_FRIS_HTMAX.tif')
  site_data$ht_max_mean_100m = as.numeric(computer_raster_buffer_value_func(rast_ht_max, site_data, 100, mean)[,2])
  site_data$ht_max_sd_100m = as.numeric(computer_raster_buffer_value_func(rast_ht_max, site_data, 100, sd)[,2])
  
  ggplot(site_data, aes(x = age, y = ht_max_mean_100m, label = site)) +
    geom_point() + geom_text_repel(size = 2)
  ggplot(site_data, aes(x = age, y = ht_max_sd_100m, label = site)) +
    geom_point() + geom_text_repel(size = 2)
  
  mapview(rast_ht_lorey) +
    mapview(site_data, zcol = 'ht_max_mean_100m', legend = T)

# Canopy layers (#)
  
  # Local: TODO
  
  # Regional
  rast_canopy_layers = load_raster('data/environment/rs_fris/rs_fris_CANOPY_LAYERS/RS_FRIS_CANOPY_LAYERS.tif')
  site_data$canopy_layers_mean_100m = as.numeric(computer_raster_buffer_value_func(rast_canopy_layers, site_data, 100, mean)[,2])
  site_data$canopy_layers_sd_100m = as.numeric(computer_raster_buffer_value_func(rast_canopy_layers, site_data, 100, sd)[,2])
  
  ggplot(site_data, aes(x = age, y = canopy_layers_mean_100m, label = site)) +
    geom_point() + geom_text_repel(size = 2)
  ggplot(site_data, aes(x = age, y = canopy_layers_sd_100m, label = site)) +
    geom_point() + geom_text_repel(size = 2)
  
  mapview(rast_canopy_layers) +
    mapview(site_data, zcol = 'canopy_layers_mean_100m', legend = T)

# Canopy cover
  
  # Local: TODO
  
  # Regional
  rast_canopy_cover = load_raster('data/environment/rs_fris/rs_fris_COVER/RS_FRIS_COVER.tif')
  site_data$canopy_cover_mean_100m = as.numeric(computer_raster_buffer_value_func(rast_canopy_cover, site_data, 100, mean)[,2])
  site_data$canopy_cover_sd_100m = as.numeric(computer_raster_buffer_value_func(rast_canopy_cover, site_data, 100, sd)[,2])
  
  ggplot(site_data, aes(x = age, y = canopy_cover_mean_100m, label = site)) +
    geom_point() + geom_text_repel(size = 2)
  ggplot(site_data, aes(x = age, y = canopy_cover_sd_100m, label = site)) +
    geom_point() + geom_text_repel(size = 2)
  
  mapview(rast_canopy_cover) +
    mapview(site_data, zcol = 'canopy_cover_mean_100m', legend = T)
  
# Canopy closure
  
  # Regional
  rast_canopy_closure = load_raster('data/environment/rs_fris/rs_fris_CLOSURE/RS_FRIS_CLOSURE.tif')
  site_data$canopy_closure_mean_100m = as.numeric(computer_raster_buffer_value_func(rast_canopy_closure, site_data, 100, mean)[,2])
  site_data$canopy_closure_sd_100m = as.numeric(computer_raster_buffer_value_func(rast_canopy_closure, site_data, 100, sd)[,2])
  
  ggplot(site_data, aes(x = age, y = canopy_closure_mean_100m, label = site)) +
    geom_point() + geom_text_repel(size = 2)
  ggplot(site_data, aes(x = age, y = canopy_closure_sd_100m, label = site)) +
    geom_point() + geom_text_repel(size = 2)
  
  mapview(rast_canopy_closure) +
    mapview(site_data, zcol = 'canopy_closure_mean_100m', legend = T)

# Number/density of all snags
  
  # Local: TODO
  
  # Regional
  rast_snags = load_raster('data/environment/rs_fris/rs_fris_SNAG_ACRE_15/RS_FRIS_SNAG_ACRE_15.tif')
  site_data$snags_mean_100m = as.numeric(computer_raster_buffer_value_func(rast_snags, site_data, 100, mean)[,2])
  site_data$snags_sd_100m = as.numeric(computer_raster_buffer_value_func(rast_snags, site_data, 100, sd)[,2])
  
  ggplot(site_data, aes(x = age, y = snags_mean_100m, label = site)) +
    geom_point() + geom_text_repel(size = 2)
  ggplot(site_data, aes(x = age, y = snags_sd_100m, label = site)) +
    geom_point() + geom_text_repel(size = 2)
  
  mapview(rast_snags) +
    mapview(site_data, zcol = 'snags_mean_100m', legend = T)

# Understory vegetation density/cover/volume
  
  # Local: TODO
  
  # Regional: TODO
  
# Dead and downed woody material
  
  # Local: TODO
  
  # Regional
  rast_ddwm = load_raster('data/environment/rs_fris/rs_fris_CFVOL_DDWM/RS_FRIS_CFVOL_DDWM.tif')
  site_data$ddwm_mean_100m = as.numeric(computer_raster_buffer_value_func(rast_ddwm, site_data, 100, mean)[,2])
  site_data$ddwm_sd_100m = as.numeric(computer_raster_buffer_value_func(rast_ddwm, site_data, 100, sd)[,2])
  
  ggplot(site_data, aes(x = age, y = ddwm_mean_100m, label = site)) +
    geom_point() + geom_text_repel(size = 2)
  ggplot(site_data, aes(x = age, y = ddwm_sd_100m, label = site)) +
    geom_point() + geom_text_repel(size = 2)
  
  mapview(rast_ddwm) +
    mapview(site_data, zcol = 'ddwm_mean_100m', legend = T)

# Tree species composition/richness/diversity
  
  # Local: TODO
  
  # Regional
  rast_bap_hwd = load_raster('data/environment/rs_fris/rs_fris_BAP_HWD/RS_FRIS_BAP_HWD.tif')
  site_data$bap_hwd_mean_100m = as.numeric(computer_raster_buffer_value_func(rast_bap_hwd, site_data, 100, mean)[,2])
  site_data$bap_hwd_sd_100m = as.numeric(computer_raster_buffer_value_func(rast_bap_hwd, site_data, 100, sd)[,2])
  
  ggplot(site_data, aes(x = age, y = bap_hwd_mean_100m, label = site)) +
    geom_point() + geom_text_repel(size = 2)
  ggplot(site_data, aes(x = age, y = bap_hwd_sd_100m, label = site)) +
    geom_point() + geom_text_repel(size = 2)
  
  mapview(rast_bap_hwd) +
    mapview(site_data, zcol = 'bap_hwd_mean_100m', legend = T)

### Patch scale variables

# Area (ha)

# Core area index (%)

# Distance to edge (m)

# Edge contrast index (%)

### Landscape scale variables (in a buffered radius)

# Proportional abundance of each patch class (%)

# Patch diversity/evenness/richness (#)

# Patch similarity/isolation index (#)

# Edge density

# Density of roads (km/km2)

# Density of streams (km/km2)

####################################################################################################
# Local habitat data

plot_data_path = 'data/environment/PAM_PreHarvest_Habitat_results_DD_WD_TM.xlsx'

# Load and clean plot-level data
plot_data = list(
  readxl::read_xlsx(plot_data_path, sheet = 2, skip = 1) %>% janitor::clean_names(),
  readxl::read_xlsx(plot_data_path, sheet = 4, skip = 1) %>% janitor::clean_names(),
  readxl::read_xlsx(plot_data_path, sheet = 5, skip = 1) %>% janitor::clean_names(),
  readxl::read_xlsx(plot_data_path, sheet = 6, skip = 1) %>% janitor::clean_names()
) %>%
  reduce(full_join, by = c("station", "strata")) %>%
  rename(site = station, stratum = strata) %>%
  mutate(
    tag = str_extract(site, "_.*$") %>% str_remove("^_"),
    site = str_remove(site, "_.*$")
  )

# Coalesce duplicate site entries
plot_data = plot_data %>% group_by(site) %>% summarise(across(everything(), ~ coalesce(.[!is.na(.)][1], NA)))

# NA entries by column
na_values =  t(as.data.frame(plot_data %>% summarise(across(everything(), ~ sum(is.na(.))))))

# TODO: Override NAs with 0 for select columns?
# selected_cols = c(avg_dbh_cm, avg_height_cm)
# plot_data %>% mutate(across(all_of(selected_cols), ~replace(., is.na(.), 0)))

cor_data = plot_data %>% select(-c(site, stratum, tag))

# Drop "unknown" variables
cor_data = cor_data %>% select(-ends_with("_un"))
# Drop variables with sparse observations and low variance
cor_data = cor_data %>% select(-c(avg_dbh_cm_abam, avg_dbh_cm_pisi, avg_height_m_abam, avg_height_m_alru, avg_height_m_pisi, avg_hlc_abam, avg_hlc_alru, avg_hlc_pisi, avg_llc_abam, avg_llc_alru, avg_llc_pisi, avg_lcr_abam, avg_lcr_alru, avg_lcr_pisi))
# i = 75
# cor(cor_data[1:i], use = "pairwise.complete.obs", method = "pearson")
# colnames(cor_data)[i]

cor_matrix = cor(cor_data, use = "pairwise.complete.obs", method = "pearson")

corrplot::corrplot(cor_matrix, method = "color", type = "upper",
                   tl.cex = 0.2, tl.col = "black", tl.srt = 45, # text color and rotation
                   addCoef.col = "black", # add correlation coefficients
                   number.cex = 0.2, # size of the numbers
                   diag = FALSE) # hide diagonal

correlation_threshold = function(cm, threshold) {
  highly_correlated_idx = which(abs(cm) >= threshold & abs(cm) < 1, arr.ind = TRUE)
  highly_correlated_idx = highly_correlated_idx[highly_correlated_idx[,1] < highly_correlated_idx[,2], ]
  data.frame(
    rownames(cm)[highly_correlated_idx[,1]],
    colnames(cm)[highly_correlated_idx[,2]],
    correlation = cm[highly_correlated_idx]
  )
}

correlation_threshold(cor_matrix, 0.8)

# Remove the following variables to reduce multicollinearity
cor_data_reduced = cor_data %>% select(-c(
  # perc_deciduous_shrub
  cv_deciduous_shrub,
  # per_cover_all_shrubs
  per_cover_total_understory, per_cover_medium_shrub, shrub_layer_vol,
  # cv_all_shrubs
  cv_shrub_layer_vol,
  # avg_understory_heights
  max_understory_heights,
  # cv_avg_understory_heights
  cv_max_understory_heights,
  # avg_dbh_cm_psme, avg_dbh_cm_tshe
  avg_dbh_cm_all, avg_height_m_psme,
  # ba_ha_psme
  large_per_hectare_psme,
  # ba_ha_thpl
  avg_dbh_cm_thpl,
  # ba_ha_abam
  large_per_hectare_abam,
  # snags_ha
  decay_1_snags_ha,
  # vol_alldown_m3
  vol_rottendown_m3
))
# Discard all HLC and LCR measurements, but retain LLC
# These measurements are highly correlated and some are interdependent.
# LLC is a directly interpretable measurement of the amount of crown habitat
cor_data_reduced = cor_data_reduced %>% select(-matches("_hlc_|_lcr_"))
# Drop all species-specific large tree average heights and LLC measurements
cor_data_reduced = cor_data_reduced %>% select(-starts_with("avg_height_m_"), avg_height_m_all)
cor_data_reduced = cor_data_reduced %>% select(-matches("_llc_"), avg_llc_all)

cor_matrix_reduced = cor(cor_data_reduced, use = "pairwise.complete.obs", method = "pearson")
corrplot::corrplot(cor_matrix_reduced, method = "color", type = "upper",
                   tl.cex = 0.2, tl.col = "black", tl.srt = 45, # text color and rotation
                   addCoef.col = "black", # add correlation coefficients
                   number.cex = 0.2, # size of the numbers
                   diag = FALSE) # hide diagonal

correlation_threshold(cor_matrix_reduced, 0.7)

####################################################################################################
# Remote sensing data



