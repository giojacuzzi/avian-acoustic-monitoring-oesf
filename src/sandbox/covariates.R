library(tidyverse)
library(sf)
library(terra)
library(viridis)
library(mapview)
options(mapview.maxpixels = 9000000)

### Study area

crs_data = 4326 # EPSG:4326
crs_projected = 32610 # EPSG:32610 (for buffering in meter units)

# ARU locations
site_data = gpx::read_gpx('data/sites/PAMlocations20230510.gpx')$waypoints %>%
  janitor::clean_names() %>% rename(site = name) %>% select(site, latitude, longitude)
site_data = st_transform(st_as_sf(site_data, coords = c('longitude', 'latitude'), crs = crs_data), crs = crs_projected)

# Study area bounding buffer
study_area = st_buffer(st_as_sfc(st_bbox(site_data)), dist = 100)

# Roads
roads = st_zm(st_read('data/environment/GIS Data/T3roads.shp'), drop = TRUE)

# Strata and forest inventory units

  # Seral stage class
  seral_classes = c("initiation", "canclose", "stemex", "mature")
  poly_seral_class = st_transform(st_read('data/environment/GIS Data/Forest Development Strata/AgeStrataFixedDIS_RSFRIS20200130.shp'), crs = crs_projected) %>% janitor::clean_names() %>%
    mutate(class_seral = stratum %>% str_to_lower())
  poly_seral_class = poly_seral_class %>%
    mutate(stratum = factor(stratum, levels = seral_classes)) %>% select(-stratum)
  # rsfris_units = st_set_crs(st_read('data/environment/GIS Data/Forest Development Strata/AgeStrataFixedRSFRIS20200130.shp'), 'NAD83')
  
  # Thinning treatment
  poly_thinning_treatment = st_transform(st_read('data/environment/GIS Data/Forest Development Strata/ThinAfter94NoHarvSinceClipByInitBuf3.shp'), crs = crs_projected) %>%
    select(TECHNIQUE_, FMA_DT, FMA_STATUS) %>%
    janitor::clean_names() %>%
    mutate(technique = technique %>% str_to_lower(), fma_status = fma_status %>% str_to_lower()) %>%
    rename(thinning_treatment = technique, thinning_status = fma_status, thinning_date = fma_dt)

# Classify sites based on the seral class and thinning prescription polygons in which they reside
site_data = st_join(site_data, poly_seral_class)
table(site_data$class_seral, useNA = 'ifany')
site_data = st_join(site_data, poly_thinning_treatment)
table(site_data$thinning_treatment, useNA = 'ifany')

# Forest origin (age)
rast_origin = rast('data/environment/rs_fris/rs_fris_Origin_Year/RS_FRIS_ORIGIN_YEAR.tif') # TODO: GET RASTER DATED FOR 2020!
crop_vector = vect(st_transform(study_area, crs(rast_origin)))
rast_origin = mask(crop(rast_origin, crop_vector), crop_vector)
sort(unique(terra::values(rast_origin)))

rast_age = 2020 - rast_origin

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
site_buffers = st_buffer(st_transform(site_data, crs(rast_age_class)), dist = 100)
surrounding_classes = terra::extract(rast_age_class, vect(site_buffers), fun = modal, na.rm = TRUE)
site_data$class_age = as.factor(surrounding_classes[, 2])
table(as.numeric(site_data$class_age))

# Get mean age
site_buffers = st_buffer(st_transform(site_data, crs(rast_age)), dist = 100)
surrounding_classes = terra::extract(rast_age, vect(site_buffers), fun = modal, na.rm = TRUE)
site_data$age = surrounding_classes[, 2]
table(as.numeric(site_data$age))

class_colors = viridis(4, option = 'D', direction = -1)

mapview(study_area, col.regions = "transparent", color = "black", lwd = 2) +
  mapview(as.factor(rast_age_class), col.regions = class_colors) +
  mapview(poly_seral_class, col.regions = class_colors) +
  mapview(site_data, zcol = 'class_age', col.regions = class_colors, legend = T) +
  mapview(site_data, zcol = 'class_seral', col.regions = class_colors, legend = T) +
  mapview(site_data, zcol = 'thinning_treatment', col.regions = class_colors, legend = T)

### Plot scale variables

# Elevation (m)

  # TODO

# Total basal area (m2/ha)

  # Local

  # Regional
  rast_ba = rast('data/environment/rs_fris/rs_fris_BA/RS_FRIS_BA.tif')
  crop_vector = vect(st_transform(study_area, crs(rast_ba)))
  rast_ba = mask(crop(rast_ba, crop_vector), crop_vector)
  
  mapview(rast_ba) +
    mapview(site_data, zcol = 'forest_class', col.regions = class_colors, legend = T)

# Tree height distribution (e.g. SD of height) (m)

# Canopy layers (#)

# Canopy cover

# Number/density of all snags

# Understory vegetation density/cover/volume

# Tree species composition/richness/diversity

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



