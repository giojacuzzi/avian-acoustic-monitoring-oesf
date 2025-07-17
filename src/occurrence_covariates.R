##############################################################################
# Quantify covariates on occurrence
#
# Note "hs" refers to data collected from in-person habitat surveys, while
# "rs" refers to data derived via remote-sensing imagery.
##############################################################################

overwrite_rast_cover_cache = FALSE
path_rast_cover_clean_out = "data/cache/occurrence_covariates/rast_cover_clean.tif"

library(tidyverse)
library(sf)
library(terra)
library(mapview)
options(mapview.maxpixels = 2117676)
library(ggrepel)
theme_set(theme_minimal())
library(landscapemetrics)
library(lwgeom)
library(units)
library(viridis)
library(vegan)

##############################################################################
# Study area, sites, boundaries, and helper functions

crs_m = 32610 # coordinate reference system EPSG:32610, UTM Zone 10N, (meters)
crs_m_rast = "EPSG:32610"

# ARU locations (sites)
aru_sites = st_read('data/environment/GIS Data/AcousticStations.shp') %>%
  st_drop_geometry() %>% janitor::clean_names() %>% as.data.frame() %>%
  st_as_sf(coords = c("utm_e", "utm_n"), crs = crs_m) %>%
  select(name, ces, treatment, geometry) %>% rename(site = name)
mapview(aru_sites, label = aru_sites$name)

# Study area boundary buffer
study_area = st_buffer(st_as_sfc(st_bbox(aru_sites)), dist = 100)

# Load a raster, cropped, projected, and masked to the study area
load_raster = function(path_rast) {
  r = rast(path_rast)
  sa_r = vect(st_transform(study_area, crs(r)))
  r = crop(r, sa_r)
  r = mask(r, sa_r)
  r = project(r, crs_m_rast)
  # r_proj = project(r, crs_m_rast)
  # r_crop = crop(r_proj, vect(study_area))
  # r_mask = mask(r_crop, vect(study_area))
  return(r)
}

# Compute raster function value for a buffer
compute_raster_buffer_value_func = function(raster, points, buffer_distance, func) {
  r = project(raster, crs(points))
  site_buffers = st_buffer(points, dist = buffer_distance)
  buffer_values = terra::extract(r, vect(site_buffers), fun = func, na.rm = TRUE)
}

compute_cv = function(x, na.rm = TRUE) {
  if (na.rm) x = x[!is.na(x)]
  if (length(x) == 0 || mean(x) == 0) {
    return(NA)
  } else {
    return(sd(x) / mean(x))
  }
}

# Conversion factors
conv_sqftAcre_m2Ha = 0.229568 # ft2/acre to m2/ha
conv_acre_hectare = 2.47105 # acre to hectare
conv_hectare_acre = 0.404686 # hectare to acre
conv_in_cm = 2.54 # inches to centimeters
conv_ft_m = 0.3048 # feet to meters
conv_ft3_m3 = 0.0283168 # cubic feet to cubic meters
conv_m2_ha = 0.0001 # square meter to hectare

# Base RS-FRIS data (0.1 acre resolution, i.e. ~404m2 or 20.10836 * 20.10836 sides, roughly 1% of the area of a 100m radius circle)
# RS-FRIS 4.0 uses a combination of 2019 and 2020 photogrammetry.
# RS-FRIS 5.0 uses a combination of 2021 and 2022 photogrammetry. 
# TODO: Calculate on a yearly basis!
dir_rsfris_version = 'data/environment/rsfris_v4.0' # Only use 2020 for now

##############################################################################
# Patch delineation
#
# Discrete habitat patches are delineated by forest stand developmental stage classes according to O'Hara et al. 1996, Oliver and Larson 1996, and Spies 1997, namely: stand initiation, stem exclusion, understory reinitiation, and old forest
# Patch boundaries are further delineated by paved roads and watercourses/waterbodies.

# Store site-specific covariate data
site_data = aru_sites
site_data$hs = FALSE # flag for habitat survey data availability

seral_classes = c("stand initiation", "canopy closure", "stem exclusion", "mature/old forest")

# WADNR delineated patch vectors. Manually determined from a combination of aerial imagery and remote sensing.
vect_patches = st_read('data/environment/GIS Data/Forest Development Strata/AgeStrataFixedDIS_RSFRIS20200130.shp') %>%
  st_transform(crs_m) %>% janitor::clean_names() %>% mutate(wadnr_stage_vect = stratum %>% str_to_lower()) %>%
  mutate(stratum = factor(stratum, levels = seral_classes)) %>% select(-stratum)
site_data = st_join(site_data, vect_patches)

# Origin year is calculated from multiple data sources using the following logic: Origin year is reported at the pixel-level (1/10th ac scale). The default value is the predicted origin year, based on RS-FRIS models. These models rely on remotely sensed data (LiDAR and DAP) and are constructed primarily from height-to-age relationships. The default RS-FRIS predicted ages are overwritten using the following data, if available:
# - Tree core data from DNR's historic inventory (FRIS) is used for stands whose origin year is 1900 or earlier.
# - FRIS tree core data from younger stands (post-1900) or FRIS data based on stratified samples is used only if RS-FRIS data is not available.
# - Information from the Land Resource Manager system (LRM) is used for completed harvests ('TEMP_RET_REM', 'TEMP_RET_1ST', 'VRH', 'SEEDTREE_INT', 'CLEAR_CUT', 'PATCH_REGEN', 'LANDUSE_CONV', 'SEEDTREE_REM', 'SHELTER_INT', 'SHELTER_REM'). STAND_ORIGIN_DT was used if populated, otherwise FMA_DT (<forest management activity date?>).
# Data sources are mapped in the "Combined Origin Year Data Source" raster.
rast_origin = load_raster(paste0(dir_rsfris_version, '/RS_FRIS_ORIGIN_YEAR.img'))
summary((terra::values(rast_origin)))

# Flag patches of missing origin year data
# TODO: Impute these patches from canopy cover, size class, canopy layers, and surrounding classes (see Powell vegetation stages white paper). Also consult ESRI World Imagery Wayback for visual inspection.
rast_origin_missing = rast_origin
values(rast_origin_missing)[values(rast_origin_missing) <= 2020] = NA
values(rast_origin)[values(rast_origin) > 2020] = NA

rast_age = round(2020 - rast_origin)

# Classify stand developmental stage classes according to O'Hara et al. 1996, Oliver and Larson 1996, and Spies 1997
stage_classes = c("stand initiation", "stem exclusion", "understory reinitiation", "old forest")
rast_stage_class = classify(rast_age, matrix(c(
  0,   25,  1,   # "stand initiation"
  25,  80,  2,   # "stem exclusion"
  80, 200,  3,   # "understory reinitiation"
  200, Inf, 4    # "old forest"
), ncol = 3, byrow = TRUE), include.lowest = TRUE, right = FALSE)
rast_stage_class = as.factor(rast_stage_class)
unique(terra::values(rast_stage_class))
site_data$stage = extract(rast_stage_class, vect(site_data))[, 2]
site_data$stage = factor(site_data$stage, labels = stage_classes)

# Classify stand seral stage classes
wadnr_classes = c("stand initiation", "canopy closure", "stem exclusion", "mature/old forest") # classify origin raster based on number of years since 2020
rast_wadnr_class = classify(rast_age, matrix(c(
  0,   15,  1,   # "stand initiation"
  15,  25,  2,   # "canopy closure"
  25,  80,  3,   # "stem exclusion"
  80, Inf, 4     # "mature/old forest"
), ncol = 3, byrow = TRUE), include.lowest = TRUE, right = FALSE)
rast_wadnr_class =as.factor(rast_wadnr_class)
unique(terra::values(rast_wadnr_class))
site_data$wadnr_stage_rast = extract(rast_wadnr_class, vect(site_data))[, 2]
site_data$wadnr_stage_rast = factor(site_data$wadnr_stage_rast, labels = wadnr_classes)

# TODO: Resolve the small number of sites that are ambiguously classified due to origin age near class boundaries:
# Dz303i > age 24, stand init / canopy closure / WADNR stemex
# Ca263i > age 86, understory reinit, WADNR stemex
# Dp166i > age 24, stand init / canopy closure / WANDR stemex

# From here forward, we use "rast_stage_class" as the classification scheme to delineate patches
summary(site_data$stage)
mapview(rast_stage_class) + mapview(aru_sites, label = aru_sites$site) + mapview(st_buffer(aru_sites, 100), col.regions = 'transparent', lwd = 2)

# Further delineate patch boundaries by roads, waterbodies, and watercourses

# Paved roads. Unpaved roads do not constitute a boundary of a patch because they are narrow, rarely driven and likely experienced by the birds as gaps in the forest.
roads_wadnr = st_read('data/environment/GIS Data/T3roads.shp') %>%
  st_zm(drop = TRUE) %>% janitor::clean_names() %>% st_transform(crs = crs_m) %>% select(road_usgs1, road_surfa, geometry) %>%
  mutate(geometry = st_cast(geometry, "MULTILINESTRING"))
roads_wsdot = st_read('data/environment/WSDOT_-_Local_Agency_Public_Road_Routes/WSDOT_-_Local_Agency_Public_Road_Routes.shp') %>% filter(RouteIdent %in% c("400000220i", "031265969i")) %>% st_transform(crs_m) %>% select(geometry) # get Hoh Mainline Road / Clearwater Road from WSDOT
roads_wsdot$road_usgs1 = 'Primary Highway'
roads_wsdot$road_surfa = NA
roads = rbind(roads_wsdot, roads_wadnr)
mapview(roads, zcol = 'road_usgs1') + mapview(aru_sites, label = aru_sites$site)

# "The Hoh-Clearwater Mainline is our only double-lane mainline, and it is about 26 feet wide for the asphalt surface. If we’re looking at right of way widths, we could easily assume a minimum width of about 50 feet for a 12-foot road, probably a good 60-80 feet for a 14-20 foot wide road, and about 100 feet for the Hoh Mainline. 100 feet might be good for US 101 as well for right of way width, and maybe about 30 feet for actual road surface width."
paved_primary_roads   = st_make_valid(roads %>% filter(road_usgs1 %in% c("Primary Highway")))
road_half_width_primary = max(100 / 2 * conv_ft_m, res(rast_stage_class)[1] / 2)
# mapview(paved_primary_roads) + mapview(st_union(st_buffer(paved_primary_roads, road_half_width_primary)))
# "Our minimum road surface width is going to be 12 feet. That would cover most of our roads. Main arterials/single lane mainlines will have a minimum surface width of 14 feet, with a few up to 20 feet."
paved_secondary_roads = st_make_valid(roads %>% filter(road_usgs1 %in% c("Light-Duty Road")))
road_half_width_secondary = max(30 / 2 * conv_ft_m, res(rast_stage_class)[1] / 2)
# mapview(paved_secondary_roads) + mapview(st_union(st_buffer(paved_secondary_roads, road_half_width_secondary)))

template = rast(ext(rast_stage_class), resolution = res(rast_stage_class), crs = crs(rast_stage_class))
paved_roads_primary_buffered = (st_buffer(paved_primary_roads, dist = road_half_width_primary))
paved_roads_secondary_buffered = (st_buffer(paved_secondary_roads, dist = road_half_width_secondary))
rast_paved_roads_primary = rasterize(vect(paved_roads_primary_buffered), template, field = 5, background = NA, touches=TRUE)
rast_paved_roads_secondary = rasterize(vect(paved_roads_secondary_buffered), template, field = 5, background = NA, touches=TRUE)
# mapview(rast_paved_roads_primary) + mapview(paved_roads_primary_buffered)
# mapview(rast_paved_roads_secondary) + mapview(paved_roads_secondary_buffered)
# mapview(rast_paved_roads_primary) + mapview(rast_paved_roads_secondary)
rast_updated = rast_stage_class
rast_updated = cover(rast_paved_roads_primary, rast_updated)
rast_updated = cover(rast_paved_roads_secondary, rast_updated)

# Watercourses/waterbodies (rivers and streams) from Type 1-3. These support distinct riparian vegetation which constitutes different habitat. Type 4 (non-fish-bearing streams) and type 5 are not boundaries because they are small features that do not change the vegetation, support less aquatic biota, and often are not permanent.
watercourses = st_read('data/environment/DNR_Hydrography/DNR_Hydrography_-_Watercourses_-_Forest_Practices_Regulation/DNR_Hydrography_-_Watercourses_-_Forest_Practices_Regulation.shp')
watercourses = watercourses %>% st_crop(st_transform(study_area, st_crs(watercourses))) %>% st_transform(crs_m) %>% janitor::clean_names() %>% select(geometry, sl_wtrty_c)

# mapview(watercourses, zcol = 'sl_wtrty_c')

boundary_watercourses = watercourses %>% filter(sl_wtrty_c %in% c(1, 2, 3))
boundary_watercourses$sl_wtrty_c = 1
watercourse_half_width = res(rast_stage_class)[1] / 2
watercourses_buffered = (st_buffer(boundary_watercourses, dist = watercourse_half_width))
rast_watercourses = rasterize(vect(watercourses_buffered), template, field = 6, background = NA, touches=TRUE)
# mapview(rast_watercourses) + mapview(watercourses_buffered)
rast_updated = cover(rast_watercourses, rast_updated)

waterbodies = st_read('data/environment/DNR_Hydrography/DNR_Hydrography_-_Water_Bodies_-_Forest_Practices_Regulation/DNR_Hydrography_-_Water_Bodies_-_Forest_Practices_Regulation.shp')
waterbodies = waterbodies %>% st_crop(st_transform(study_area, st_crs(waterbodies))) %>% st_transform(crs_m) %>% janitor::clean_names()
# mapview(waterbodies, zcol = 'sl_wtrty_c')

boundary_waterbodies = waterbodies %>% filter(sl_wtrty_c %in% c(1))
rast_waterbodies = rasterize(vect(boundary_waterbodies), template, field = 6, background = NA, touches=TRUE)
# mapview(rast_waterbodies) + mapview(boundary_waterbodies)
rast_updated = cover(rast_waterbodies, rast_updated)
mapview(rast_updated)

# Clean raster by generalizing minimum patch area and width
if (overwrite_rast_cover_cache) {
  
  r = rast_updated
  cell_res <- res(r)[1] # cell resolution (in meters)
  min_area <- 0.785 * 1e4 # hectares to square meters
  min_width <- 50 # in meters
  min_cells_area <- ceiling(min_area / (cell_res^2))  # min number of cells by area
  min_cells_width <- ceiling(min_width / cell_res)    # min number of cells by width
  
  # Identify discrete contiguous patches with unique ids
  patch_list <- list()
  start_id <- 1  # Starting patch ID to ensure uniqueness
  cls = sort(na.omit(unique(values(r))))
  for (cl in cls) {
    r_class <- mask(r, r == cl, maskvalues = FALSE) # Mask only the current class
    r_patches <- patches(r_class, directions = 8) # Compute patches for this class only
    # Reclassify patch IDs to be globally unique
    max_patch_id <- global(r_patches, "max", na.rm = TRUE)[1,1]
    if (!is.na(max_patch_id)) {
      r_patches <- classify(r_patches, cbind(1:max_patch_id, start_id:(start_id + max_patch_id - 1)))
      start_id <- start_id + max_patch_id
    }
    patch_list[[as.character(cl)]] <- r_patches
  }
  
  # Merge patch ids for each raster type into one raster
  patch_ids <- cover(patch_list[[1]], patch_list[[2]])
  patch_ids <- cover(patch_ids, patch_list[[3]])
  patch_ids <- cover(patch_ids, patch_list[[4]])
  patch_ids <- cover(patch_ids, patch_list[[5]])
  patch_ids <- cover(patch_ids, patch_list[[6]])
  
  # Create a raster stack (cover class and patch id)
  patch_stack <- c(patch_ids, r)
  names(patch_stack) <- c("patch_id", "cover_class")
  
  mapview(patch_stack[["patch_id"]], col.regions = viridis) + 
    mapview(patch_stack[["cover_class"]])
  
  # Identify patches with maximum width less than minimum width ##################################################
  narrow_patches <- c()
  pids = unique(na.omit(as.vector(values(patch_stack[['patch_id']]))))
  for (pid in pids) {
    print(pid)
    m = trim(classify(patch_stack[['patch_id']], cbind(pid, 1), others = NA)) # Mask patch
    m_vals <- as.matrix(m, wide=TRUE)
    
    if (all(is.na(m_vals))) next  # Skip empty masks
    
    max_h_run <- max(apply(m_vals, 1, function(row) { # Horizontal runs (per row)
      rle_row <- rle(row)
      max(c(0, rle_row$lengths[rle_row$values == 1]), na.rm = TRUE)
    }), na.rm = TRUE)
    max_v_run <- max(apply(m_vals, 2, function(col) { # Vertical runs (per column)
      rle_col <- rle(col)
      max(c(0, rle_col$lengths[rle_col$values == 1]), na.rm = TRUE)
    }), na.rm = TRUE)
    if (max_h_run < min_cells_width || max_v_run < min_cells_width) { # Keep if both directions have max run < 3
      narrow_patches <- c(narrow_patches, pid)
    }
  }
  
  # Identify patch ids with negligible patch width
  patch_stack[["patch_narrow"]] = classify(patch_stack[["patch_id"]], rcl = cbind(narrow_patches, narrow_patches), others = NA)
  patch_narrow_ids = unique(patch_stack[["patch_narrow"]])[,1]
  
  # Reclassify those cells to value 0
  values(patch_stack[["cover_class"]])[which(values(!is.na(patch_stack[['patch_narrow']])))] = 0
  mapview(patch_stack[['cover_class']])
  
  # Identify patches with area less than minimum area ##################################################
  
  # Identify patches with negligible area
  freq_table = freq(patch_stack[["patch_id"]])
  small_patches = freq_table[freq_table$count < min_cells_area, "value"]
  patch_stack[["patch_small"]] = classify(patch_stack[["patch_id"]], rcl = cbind(small_patches, small_patches), others = NA)
  patch_small_ids = unique(patch_stack[["patch_small"]])[,1]
  
  # Reclassify those cells to value 0
  values(patch_stack[["cover_class"]])[which(values(!is.na(patch_stack[['patch_small']])))] = 0
  mapview(patch_stack[['cover_class']])
  
  r_zero_alt <- patch_stack[['cover_class']]
  values(r_zero_alt)[values(r_zero_alt) != 0] <- NA
  mapview(r_zero_alt)
  p = patches(r_zero_alt, directions = 8, values=TRUE)
  
  # Reclassify (generalize) negligible patches as the mode of adjacent nearest neighbors
  rast_cover_clean <- patch_stack[["cover_class"]]
  i = 1
  patch_ids_to_generalize = unique(na.omit(values(p)))
  total = length(patch_ids_to_generalize)
  for (small_patch_id in patch_ids_to_generalize) {
    print(round(i / total, 3))
    
    # mapview(trim(classify(patch_stack[["patch_id"]], cbind(small_patch_id, 1), others = NA)))
    
    patch_cells = which(values(p) == small_patch_id)
    
    # Find adjacent cells to the patch
    adj_cells = adjacent(rast_stage_class, cells = patch_cells, directions = 4, pairs = TRUE)
    
    # Remove pairs where neighbor is also part of the same patch
    neighbor_cells = adj_cells[, 2]
    neighbor_cells = unique(neighbor_cells[!neighbor_cells %in% patch_cells])
    
    # Get land cover class values of neighbor cells
    neighbor_classes = values(r)[neighbor_cells]
    
    # Remove NA neighbors (e.g. outside raster or masked)
    neighbor_classes = neighbor_classes[!is.na(neighbor_classes)]
    
    if (length(neighbor_classes) == 0) {
      majority_class = NA
    } else {
      majority_class = as.numeric(names(sort(table(neighbor_classes), decreasing = TRUE)[1]))
    }
    
    # Replace the small patch cells with the majority class
    values(rast_cover_clean)[patch_cells] <- majority_class
    
    i = i + 1
  }
  
  mapview(rast_cover_clean) + mapview(patch_stack[['cover_class']])

  dir.create(dirname(path_rast_cover_clean_out), recursive = TRUE, showWarnings = FALSE)
  writeRaster(rast_cover_clean, path_rast_cover_clean_out, overwrite=TRUE)
  
} else { # overwrite_rast_cover_cache is FALSE
  rast_cover_clean = rast(path_rast_cover_clean_out)
}
rast_cover_clean = as.factor(rast_cover_clean)

mapview(rast_cover_clean,
        alpha.regions = 1.0,
        col.regions = c('#90c6bd', '#3c8273', '#d8c18a', '#9b652b', 'darkgray', '#6495ED'))

##############################################################################
# Local plot scale covariates

# Plot-level spatial scale buffer
plot_buffer = 100 # 100 meters

path_data_plot = 'data/environment/PAM_PreHarvest_Habitat_results_DD_WD_TM.xlsx'

# Load and clean plot-level data
data_plot = list(
  readxl::read_xlsx(path_data_plot, sheet = 2, skip = 1) %>% janitor::clean_names(),
  readxl::read_xlsx(path_data_plot, sheet = 4, skip = 1) %>% janitor::clean_names(),
  readxl::read_xlsx(path_data_plot, sheet = 5, skip = 1) %>% janitor::clean_names(),
  readxl::read_xlsx(path_data_plot, sheet = 6, skip = 1) %>% janitor::clean_names()
) %>%
  reduce(full_join, by = c("station", "strata")) %>% rename(site = station, stratum = strata) %>%
  mutate( tag = str_extract(site, "_.*$") %>% str_remove("^_"), site = str_remove(site, "_.*$"))
data_plot = data_plot %>% group_by(site) %>% summarise(across(everything(), ~ coalesce(.[!is.na(.)][1], NA))) # coalesce duplicate site entries

# Mark sites that were surveyed for habitat data in-person
site_data[site_data$site %in% data_plot$site, 'hs'] = TRUE
mapview(site_data, zcol = 'hs')

# Elevation [m] (derived from WA state LiDAR flights 3x3-ft raster grid)
elevation = read.csv('data/environment/site_elevation_ft.csv')
elevation$plot_elev_rs = elevation$elev_ft * conv_ft_m
site_data = site_data %>% left_join(elevation %>% select(site, plot_elev_rs), by = 'site')

# Age (mean and cv) [#]
site_data$age_point = extract(rast_age, vect(site_data))[, 2]
site_data$age_mean = as.numeric(compute_raster_buffer_value_func(rast_age, site_data, plot_buffer, mean)[,2])
site_data$age_cv = as.numeric(compute_raster_buffer_value_func(rast_age, site_data, plot_buffer, compute_cv)[,2])
hist(site_data$age_point, breaks = seq(0, max(site_data$age) + 10, by = 10))
hist(site_data$age_mean, breaks = seq(0, max(site_data$age) + 10, by = 10))
table(site_data$age_point)
table(round(site_data$age_mean))

ggplot(site_data, aes(x = stage, y = age_mean)) + geom_boxplot() + theme_minimal()

# Cover class [categorical]
site_data$stage

# Thinning treatment [categorical]
poly_thinning_treatment = st_read('data/environment/GIS Data/Forest Development Strata/ThinAfter94NoHarvSinceClipByInitBuf3.shp') %>% 
  st_transform(crs_m) %>% select(TECHNIQUE_, FMA_DT, FMA_STATUS) %>% janitor::clean_names() %>%
  mutate(technique = technique %>% str_to_lower(), fma_status = fma_status %>% str_to_lower()) %>%
  rename(thinning_treatment = technique, thinning_status = fma_status, thinning_date = fma_dt)
site_data = st_join(site_data, poly_thinning_treatment)
table(site_data$thinning_treatment, useNA = 'ifany')

# Basal area [m2/ha] TODO: confirm if this is a mean value at local level
site_data = site_data %>% left_join(data_plot %>% select(site, plot_ba_hs = ba_ha_all), by = 'site')

rast_ba = load_raster(paste0(dir_rsfris_version, '/RS_FRIS_BA_EXPORT.tif'))
site_data$plot_ba_rs = as.numeric(compute_raster_buffer_value_func(rast_ba, site_data, plot_buffer, mean)[,2]) * conv_sqftAcre_m2Ha
site_data$plot_ba_sd_rs = as.numeric(compute_raster_buffer_value_func(rast_ba, site_data, plot_buffer, sd)[,2]) * conv_sqftAcre_m2Ha

ggplot(site_data, aes(x = plot_ba_hs, y = plot_ba_rs, label = site)) +
  geom_point() + geom_text_repel(size = 2) + geom_abline(slope = 1) +
  xlim(0, max(site_data$plot_ba_hs, site_data$plot_ba_rs, na.rm = TRUE)) +
  ylim(0, max(site_data$plot_ba_hs, site_data$plot_ba_rs, na.rm = TRUE)) +
  ggtitle('Basal area [m2/ha]')
mapview(rast_ba) +
  mapview(site_data, zcol = 'plot_ba_rs', legend = T)

# Tree density (large, dbh > 10 cm) [# trees/ha]
site_data = site_data %>% left_join(data_plot %>% select(site, plot_treeden_gt10cmDbh_hs = large_per_hectare_all), by = 'site')

rast_treeden_gt4inDbh = load_raster(paste0(dir_rsfris_version, '/RS_FRIS_TREE_ACRE_4.img'))
site_data$plot_treeden_gt4inDbh_rs = as.numeric(compute_raster_buffer_value_func(rast_treeden_gt4inDbh, site_data, plot_buffer, mean)[,2]) * conv_acre_hectare

ggplot(site_data, aes(x = plot_treeden_gt10cmDbh_hs, y = plot_treeden_gt4inDbh_rs, label = site)) +
  geom_point() + geom_text_repel(size = 2) + geom_abline(slope = 1) +
  xlim(0, max(site_data$plot_treeden_gt10cmDbh_hs, site_data$plot_treeden_gt4inDbh_rs, na.rm = TRUE)) +
  ylim(0, max(site_data$plot_treeden_gt10cmDbh_hs, site_data$plot_treeden_gt4inDbh_rs, na.rm = TRUE)) +
  ggtitle('Tree density (large, dbh > 10 cm) [# trees/ha]')

mapview(rast_treeden_gt4inDbh) +
  mapview(site_data, zcol = 'plot_treeden_gt4inDbh', legend = T)

# Tree density (small, dbh < 10 cm) [# trees/ha]
site_data = site_data %>% left_join(data_plot %>% select(site, plot_treeden_lt10cmDbh_hs = small_per_hectare), by = 'site')

rast_treeden_all = load_raster(paste0(dir_rsfris_version, '/RS_FRIS_TREE_ACRE.img'))
site_data$plot_treeden_all_rs = as.numeric(compute_raster_buffer_value_func(rast_treeden_all, site_data, plot_buffer, mean)[,2]) * conv_acre_hectare
site_data$plot_treeden_lt4inDbh_rs = site_data$plot_treeden_all_rs - site_data$plot_treeden_gt4inDbh_rs

ggplot(site_data, aes(x = plot_treeden_lt10cmDbh_hs, y = plot_treeden_lt4inDbh_rs, label = site)) +
  geom_point() + geom_text_repel(size = 2) + geom_abline(slope = 1) +
  xlim(0, max(site_data$plot_treeden_lt10cmDbh_hs, site_data$plot_treeden_lt4inDbh_rs, na.rm = TRUE)) +
  ylim(0, max(site_data$plot_treeden_lt10cmDbh_hs, site_data$plot_treeden_lt4inDbh_rs, na.rm = TRUE)) +
  ggtitle('Tree density (small, dbh < 10 cm) [# trees/ha]')

# Total tree density (all sizes) [# trees/ha]
site_data$plot_treeden_all_hs = site_data$plot_treeden_gt10cmDbh_hs + site_data$plot_treeden_lt10cmDbh_hs

ggplot(site_data, aes(x = plot_treeden_all_hs, y = plot_treeden_all_rs, label = site)) +
  geom_point() + geom_text_repel(size = 2) + geom_abline(slope = 1) +
xlim(0, max(site_data$plot_treeden_all_hs, site_data$plot_treeden_all_rs, na.rm = TRUE)) +
  ylim(0, max(site_data$plot_treeden_all_hs, site_data$plot_treeden_all_rs, na.rm = TRUE)) +
  ggtitle('Tree density (total, all sizes) [# trees/ha]')

# Stand density index (Reineke's) [#]
rast_sdi = load_raster(paste0(dir_rsfris_version, '/RS_FRIS_SDI_SUM.img'))
site_data$plot_sdi_rs = as.numeric(compute_raster_buffer_value_func(rast_sdi, site_data, plot_buffer, mean)[,2])

# Tree diameter (mean and SD) [cm]
site_data = site_data %>% left_join(data_plot %>% select(site, plot_qmd_gt10cmDbh_hs = avg_dbh_cm_all), by = 'site')
site_data = site_data %>% left_join(data_plot %>% select(site, plot_qmd_lt10cmDbh_hs = avg_dbh_cm), by = 'site')
site_data$plot_qmd_all_hs = site_data$plot_qmd_gt10cmDbh_hs + site_data$plot_qmd_lt10cmDbh_hs

rast_qmd = load_raster(paste0(dir_rsfris_version, '/RS_FRIS_QMD.img')) #load_raster('data/environment/rs_fris/rs_fris_QMD/RS_FRIS_QMD.tif')
site_data$plot_qmd_rs = as.numeric(compute_raster_buffer_value_func(rast_qmd, site_data, plot_buffer, mean)[,2]) * conv_in_cm

ggplot(site_data, aes(x = plot_qmd_all_hs, y = plot_qmd_rs, label = site)) +
  geom_point() + geom_text_repel(size = 2) + geom_abline(slope = 1) +
  xlim(0, max(site_data$plot_qmd_all_hs, site_data$plot_qmd_rs, na.rm = TRUE)) +
  ylim(0, max(site_data$plot_qmd_all_hs, site_data$plot_qmd_rs, na.rm = TRUE)) +
  ggtitle('Tree diameter (quadratic mean) [cm]')
mapview(rast_qmd) +
  mapview(site_data, zcol = 'plot_qmd_rs', legend = T)

# Tree height (mean and cv) [m] and canopy layers [#]
# TODO: Get high-resolution LiDAR height data from DNR approval procedure
site_data = site_data %>% left_join(data_plot %>% select(site, plot_ht_hs = avg_height_m_all), by = 'site')
site_data = site_data %>% left_join(data_plot %>% select(site, plot_ht_cv_hs = cv_height_all), by = 'site')

rast_htmax = load_raster('data/environment/rs_fris/rs_fris_HTMAX/RS_FRIS_HTMAX.tif') # TODO: Update!
site_data$plot_htmax_rs = as.numeric(compute_raster_buffer_value_func(rast_htmax, site_data, plot_buffer, mean)[,2]) * conv_ft_m
site_data$plot_htmax_cv_rs = as.numeric(compute_raster_buffer_value_func(rast_htmax, site_data, plot_buffer, compute_cv)[,2]) * conv_ft_m * 100

ggplot(site_data, aes(x = plot_ht_hs, y = plot_htmax_rs, label = site)) +
  geom_point() + geom_text_repel(size = 2) + geom_abline(slope = 1) +
  xlim(0, max(site_data$plot_ht_hs, site_data$plot_htmax_rs, na.rm = TRUE)) +
  ylim(0, max(site_data$plot_ht_hs, site_data$plot_htmax_rs, na.rm = TRUE)) +
  ggtitle('Tree height (mean) [m]')

rast_canopy_layers = load_raster('data/environment/rs_fris/rs_fris_CANOPY_LAYERS/RS_FRIS_CANOPY_LAYERS.tif') # TODO: Update!
site_data$canopy_layers_rs = as.numeric(compute_raster_buffer_value_func(rast_canopy_layers, site_data, plot_buffer, mean)[,2])

# Canopy cover and closure [%]
rast_canopy_cover = load_raster(paste0(dir_rsfris_version, '/RS_FRIS_COVER.img'))
site_data$plot_canopy_cover_rs = as.numeric(compute_raster_buffer_value_func(rast_canopy_cover, site_data, plot_buffer, mean)[,2])

rast_canopy_closure = load_raster(paste0(dir_rsfris_version, '/RS_FRIS_CLOSURE.img'))
site_data$plot_canopy_closure_rs = as.numeric(compute_raster_buffer_value_func(rast_canopy_closure, site_data, plot_buffer, mean)[,2])

mapview(rast_canopy_cover) +
  mapview(site_data, zcol = 'plot_canopy_cover_rs', legend = T)

# Height to live crown [m], length of live crown (live crown depth) [m], and live crown ratio [#]
site_data = site_data %>% left_join(data_plot %>% select(site, plot_hlc_hs = avg_hlc_all), by = 'site')
site_data = site_data %>% left_join(data_plot %>% select(site, plot_llc_hs = avg_llc_all), by = 'site')
site_data = site_data %>% left_join(data_plot %>% select(site, plot_llc_hs = avg_lcr_all), by = 'site')

# Density of all snags [# snags/ha]
site_data = site_data %>% left_join(data_plot %>% select(site, plot_snagden_hs = snags_ha), by = 'site')

rast_snagden = load_raster(paste0(dir_rsfris_version, '/RS_FRIS_SNAG_ACRE_15.img'))
site_data$plot_snagden_rs = as.numeric(compute_raster_buffer_value_func(rast_snagden, site_data, plot_buffer, mean)[,2]) * conv_acre_hectare

ggplot(site_data, aes(x = plot_snagden_hs, y = plot_snagden_rs, label = site)) +
  geom_point() + geom_text_repel(size = 2) + geom_abline(slope = 1) +
  xlim(0, max(site_data$plot_snagden_hs, site_data$plot_snagden_rs, na.rm = TRUE)) +
  ylim(0, max(site_data$plot_snagden_hs, site_data$plot_snagden_rs, na.rm = TRUE)) +
  ggtitle('Snag density [#/ha]')

# Downed wood volume [m3/ha]
site_data = site_data %>% left_join(data_plot %>% select(site, plot_downvol_hs = vol_alldown_m3), by = 'site')

rast_downvol = load_raster(paste0(dir_rsfris_version, '/RS_FRIS_CFVOL_DDWM.img'))
site_data$plot_downvol_rs = as.numeric(compute_raster_buffer_value_func(rast_downvol, site_data, plot_buffer, mean)[,2]) * (conv_ft3_m3 / conv_hectare_acre)

ggplot(site_data, aes(x = plot_downvol_hs, y = plot_downvol_rs, label = site)) +
  geom_point() + geom_text_repel(size = 2) + geom_abline(slope = 1) +
  xlim(0, max(site_data$plot_downvol_hs, site_data$plot_downvol_rs, na.rm = TRUE)) +
  ylim(0, max(site_data$plot_downvol_hs, site_data$plot_downvol_rs, na.rm = TRUE)) +
  ggtitle('Downed wood volume [m3/ha]')

# Understory vegetation cover [%] and volume [m3/ha]
site_data = site_data %>% left_join(data_plot %>% select(site, plot_understory_cover = per_cover_total_understory), by = 'site')
site_data = site_data %>% left_join(data_plot %>% select(site, plot_understory_vol = shrub_layer_vol), by = 'site')

# Tree species richness, evenness, and diversity [#]
# Calculated from density of each species
tree_gte10cm_obs =  data_plot %>% select(site,
  large_per_hectare_psme, # Douglas fir
  large_per_hectare_thpl, # Western redcedar
  large_per_hectare_abam, # Pacific silver fir
  large_per_hectare_tshe, # Western hemlock
  large_per_hectare_alru, # Red alder
  large_per_hectare_pisi  # Sitka spruce
)
tree_gte10cm_richness  = specnumber(tree_gte10cm_obs %>% select(-site))
tree_gte10cm_diversity = diversity(tree_gte10cm_obs %>% select(-site), index = 'shannon')
tree_gte10cm_evenness  = tree_gte10cm_diversity / log(tree_gte10cm_richness)

tree_density_lt10cm = data_plot %>% pull(small_per_hectare)
tree_lt10cm_obs = data.frame(
  site = data_plot %>% pull(site),
  small_per_hectare_psme = tree_density_lt10cm * data_plot %>% pull(percent_psme) / 100,
  small_per_hectare_thpl = tree_density_lt10cm * data_plot %>% pull(percent_thpl) / 100,
  small_per_hectare_abam = tree_density_lt10cm * data_plot %>% pull(percent_abam) / 100,
  small_per_hectare_tshe = tree_density_lt10cm * data_plot %>% pull(percent_tshe) / 100,
  small_per_hectare_alru = tree_density_lt10cm * data_plot %>% pull(percent_alru) / 100,
  small_per_hectare_pisi = tree_density_lt10cm * data_plot %>% pull(percent_pisi) / 100
)
tree_lt10cm_richness  = specnumber(tree_lt10cm_obs %>% select(-site))
tree_lt10cm_diversity = diversity(tree_lt10cm_obs %>% select(-site), index = 'shannon')
tree_lt10cm_evenness  = tree_lt10cm_diversity / log(tree_lt10cm_richness)

tree_all_obs = data.frame(
  site = data_plot %>% pull(site),
  all_per_hectare_psme = tree_gte10cm_obs$large_per_hectare_psme + tree_lt10cm_obs$small_per_hectare_psme,
  all_per_hectare_thpl = tree_gte10cm_obs$large_per_hectare_thpl + tree_lt10cm_obs$small_per_hectare_thpl,
  all_per_hectare_abam = tree_gte10cm_obs$large_per_hectare_abam + tree_lt10cm_obs$small_per_hectare_abam,
  all_per_hectare_tshe = tree_gte10cm_obs$large_per_hectare_tshe + tree_lt10cm_obs$small_per_hectare_tshe,
  all_per_hectare_alru = tree_gte10cm_obs$large_per_hectare_alru + tree_lt10cm_obs$small_per_hectare_alru,
  all_per_hectare_pisi = tree_gte10cm_obs$large_per_hectare_pisi + tree_lt10cm_obs$small_per_hectare_pisi
)
tree_all_richness  = specnumber(tree_all_obs %>% select(-site))
tree_all_diversity = diversity(tree_all_obs %>% select(-site), index = 'shannon')
tree_all_evenness  = tree_all_diversity / log(tree_all_richness)

tree_div_metrics = data.frame(
  site = data_plot %>% pull(site),
  tree_all_richness,  tree_gte10cm_richness,  tree_lt10cm_richness,
  tree_all_diversity, tree_gte10cm_diversity, tree_lt10cm_diversity,
  tree_all_evenness,  tree_gte10cm_evenness,  tree_lt10cm_evenness,
  tree_all_density_psme = tree_all_obs$all_per_hectare_psme,
  tree_gte10cm_density_psme = tree_gte10cm_obs$large_per_hectare_psme,
  tree_gte10cm_density_thpl = tree_gte10cm_obs$large_per_hectare_thpl,
  tree_gte10cm_density_abam = tree_gte10cm_obs$large_per_hectare_abam,
  tree_gte10cm_density_tshe = tree_gte10cm_obs$large_per_hectare_tshe,
  tree_gte10cm_density_alru = tree_gte10cm_obs$large_per_hectare_alru,
  tree_gte10cm_density_pisi = tree_gte10cm_obs$large_per_hectare_pisi,
  tree_all_density_thpl = tree_all_obs$all_per_hectare_thpl,
  tree_all_density_abam = tree_all_obs$all_per_hectare_abam,
  tree_all_density_tshe = tree_all_obs$all_per_hectare_tshe,
  tree_all_density_alru = tree_all_obs$all_per_hectare_alru,
  tree_all_density_pisi = tree_all_obs$all_per_hectare_pisi
)

site_data = site_data %>% left_join(tree_div_metrics, by = 'site')


ggplot(site_data, aes(x = `stage`, y = tree_all_richness)) +
  geom_violin() +
  stat_summary(fun = mean, geom = "point") +
  ggtitle("Tree species richness across stages") + theme_minimal()

# psme - Douglas fir
# thpl - Western redcedar
# abam - Pacific silver fir
# tshe - Western hemlock
# alru - Red alder
# pisi - Sitka spruce
#
# Stand initiation: small alru pioneers newly-disturbed habitat as an early-seral species, otherwise plantation psme and tshe are establishing
# Stem exclusion: alru is outcompeted by psme and tshe and mostly disappears, tree size is more uniform
# Understory reinit: tshe dominates, abam establishing, tree size less uniform
# Old forest: tshe is climax species, abam established
metric = "tree_all_density_" # "tree_all_density_", "tree_gte10cm_density_"
df_long <- site_data %>% select(site, `stage`, starts_with(metric)) %>%
  pivot_longer(cols = starts_with(metric), names_to = "species", values_to = "density") %>%
  mutate(species = gsub(metric, "", species))
ggplot(df_long, aes(x = species, y = density, fill = species)) +
  geom_boxplot() +
  facet_wrap(~ stage, scales = "free_y") +
  coord_cartesian(ylim = c(0, 1250)) + 
  labs(title = "Tree species density by stage") +
  theme_minimal()

# Inspect specific variables by stage
ggplot(site_data, aes(x = stage, y = canopy_layers_rs)) +
  geom_boxplot() +
  theme_minimal()

# Distance to water (type 1-3 and all types) [m]
dist_watercourses_major = st_distance(site_data, boundary_watercourses)
dist_watercourses_major = apply(dist_watercourses_major, 1, min)
site_data$dist_watercourses_major = dist_watercourses_major
mapview(site_data, zcol = "dist_watercourses_major") + mapview(boundary_watercourses)

dist_watercourses_all = st_distance(site_data, watercourses)
dist_watercourses_all = apply(dist_watercourses_all, 1, min)
site_data$dist_watercourses_all = dist_watercourses_all
mapview(site_data, zcol = "dist_watercourses_all") + mapview(watercourses)

# Distance to nearest edge [m]
# For each site, find its patch, then find the distance to the nearest non-patch cell
patch_ids <- patches(rast_cover_clean, directions=8, values=TRUE)
site_data$dist_nearest_edge = NA
for (i in 1:nrow(site_data)) {
  d = site_data[i,]
  cover_class = extract(rast_cover_clean, vect(d))[,2]
  pid = extract(patch_ids, vect(d))[,2]
  m = trim(classify(patch_ids, cbind(pid, 1), others = NA))
  na_cells <- which(is.na(values(m)))
  na_coords <- xyFromCell(m, na_cells)
  d_coords <- st_coordinates(d)
  distances <- sqrt((na_coords[,1] - d_coords[1])^2 + (na_coords[,2] - d_coords[2])^2)
  min_dist <- min(distances)
  site_data[i, 'dist_nearest_edge'] = min_dist
  # mapview(d) + mapview(m) + mapview(st_buffer(d, 100), col.regions = 'transparent', lwd = 2)
}
mapview(rast_cover_clean) + mapview(site_data, label = site_data$site, zcol = 'dist_nearest_edge') + mapview(st_buffer(site_data, 100), col.regions = 'transparent', lwd = 2)

##############################################################################
# Focal patch scale and home range neighborhood scale covariates (limited to the extent of the species home range)

# For now, just calculate at a fixed radius plot scale (e.g. median predicted radius of 200 m)
# TODO: Calculate per species according to predicted home range size
home_range_buffer = 300 # meters

# "Edge influences on distributions of organisms or factors affecting organisms (such as predation and nest parasitism) are concentrated within 50m of the edge." (Kremaster L. & Bunnell F. L. (1999) Edge effects: Theory, evidence and implications to management of western North American forests. In: Forest Fragmentation: Wildlife and Management Implications (eds J. Wisniewski, J. A. Rochelle & L. Lehmann) pp. 117–53. Leiden, Boston, MA.)
core_area_buffer = 50 # meters
core_area_buffer_ncells = round(core_area_buffer / res(rast_cover_clean)[1], 0)

for (i in 1:nrow(site_data)) {
  site = site_data[i,]
  cover_class = extract(rast_cover_clean, vect(site))[,2]
  
  homerange_and_edge_buffer = st_buffer(site, home_range_buffer + core_area_buffer) # additional buffer to ensure that edges are retained in mask
  homerange_buffer = st_buffer(site, home_range_buffer)
  homerange_and_edge_crop = crop(rast_cover_clean, vect(homerange_and_edge_buffer))
  homerange_and_edge = mask(homerange_and_edge_crop, vect(homerange_and_edge_buffer))
  
  patch_ids = patches(homerange_and_edge, directions=8, values=TRUE)
  pid = extract(patch_ids, vect(site))[,2]
  focal_patch = trim(classify(patch_ids, cbind(pid, 1), others = NA))
  focal_patch = crop(focal_patch, vect(homerange_buffer))
  focal_patch = mask(focal_patch, vect(homerange_buffer))
  
  # Area [ha]
  lsm_p_area(focal_patch, directions = 8)[1, 'value'] # in hectares (ha)
  global(!is.na(focal_patch), sum, na.rm = TRUE) * prod(res(focal_patch)) # in sq meters
  mapview(site) +
    mapview(homerange_and_edge) + mapview(homerange_buffer, col.regions = 'transparent', lwd = 2) + mapview(focal_patch)
  
  ncells_focal_patch = sum(!is.na(values(focal_patch)))
  (focal_patch_area = ncells_focal_patch * prod(res(focal_patch)) * conv_m2_ha)
  
  # Core area [ha]
  focal_patch_and_edge_buffer = trim(classify(patch_ids, cbind(pid, 1), others = NA))
  focal_patch_inv = ifel(is.na(focal_patch_and_edge_buffer), 1, NA)
  focal_patch_edge_dist = distance(focal_patch_inv)
  focal_patch_edge_dist = ifel(focal_patch_edge_dist == 0, NA, focal_patch_edge_dist)
  focal_patch_core = ifel(focal_patch_edge_dist >= core_area_buffer, 1, NA)
  focal_patch_core = crop(focal_patch_core, vect(homerange_buffer))
  focal_patch_core = mask(focal_patch_core, vect(homerange_buffer))
  mapview(homerange_buffer, col.regions = 'transparent', lwd = 2) + mapview(focal_patch_edge_dist) + mapview(focal_patch) + mapview(focal_patch_and_edge_buffer) + mapview(focal_patch_core)
  
  ncells_focal_patch_core = sum(!is.na(values(focal_patch_core)))
  (focal_patch_core_area = ncells_focal_patch_core * prod(res(focal_patch)) * conv_m2_ha)
  
  homerange_crop = crop(rast_cover_clean, vect(homerange_buffer))
  homerange_cover = mask(homerange_crop, vect(homerange_buffer))
  homerange_cover_forest = homerange_cover
  homerange_cover_forest[!(homerange_cover_forest[] %in% c(1, 2, 3, 4))] <- NA
  mapview(homerange_cover)
  ncells_homerange = sum(!is.na(values(homerange_cover)))
  
  # Focal patch class percentage of home range [%]
  (focal_patch_pcnt = sum(!is.na(values(focal_patch))) / ncells_homerange)
  
  # Focal patch class core area percentage of home range [%]
  (focal_patch_core_pcnt = sum(!is.na(values(focal_patch_core))) / ncells_homerange)
  
  # Proportional abundance of each forest cover class [%]
  (as.numeric(global(homerange_cover == 1, fun = "sum", na.rm = TRUE)) / ncells_homerange)
  (as.numeric(global(homerange_cover == 2, fun = "sum", na.rm = TRUE)) / ncells_homerange)
  (as.numeric(global(homerange_cover == 3, fun = "sum", na.rm = TRUE)) / ncells_homerange)
  (as.numeric(global(homerange_cover == 4, fun = "sum", na.rm = TRUE)) / ncells_homerange)
  
  # Forest cover class richness and evenness (i.e. dominance) [#]
  cover_freq = freq(homerange_cover) %>% select(value, count) %>% filter(value %in% c(1,2,3,4))
  (cover_richness = nrow(cover_freq))
  lsm_l_shei(homerange_cover_forest) # 0 when only one forest cover present, 1 when equally distributed
  
  # Contagion index [#]
  # NOTE: The number of classes to calculate contagion must be >= 2
  lsm_l_contag(homerange_cover, verbose = FALSE)
  
  # Interspersion and juxtaposition index [#]
  # NOTE: The number of classes to calculate iji index must be >= 3
  lsm_l_iji(homerange_cover, verbose = FALSE)
  
  # Aggregation index [#]
  lsm_l_ai(homerange_cover, directions = 8)
  
  # Shape and fractal dimension indices [#]
  # NOTE:The number of classes to calculate fractal dimension index must be >= 10
  lsm_l_shape_mn(homerange_cover, directions = 8)
  lsm_l_pafrac(homerange_cover, directions = 8, verbose = FALSE)
  
  # Nearest neighbor distance, proximity and similarity indices [TODO: UNITS]
  
  # Contrast-weighted edge density [TODO: UNITS]
  
  # Density of roads (paved and all) [TODO: UNITS]
  homerange_roads_paved = st_intersection(st_make_valid(roads %>% filter(road_usgs1 %in% c("Primary Highway", "Light-Duty Road"))), homerange_buffer)
  homerange_roads = st_intersection(st_make_valid(roads), homerange_buffer)
  homerange_roads_paved_length = sum(st_length(homerange_roads_paved))
  homerange_roads_length = sum(st_length(homerange_roads))
  mapview(homerange_roads_paved) + mapview(homerange_roads) + mapview(homerange_buffer) + homerange_cover
  
  # Density of streams (type 1-3 and all types) [TODO: UNITS]
  homerange_streams_major = st_intersection(st_make_valid(watercourses %>% filter(sl_wtrty_c %in% c(1, 2, 3))), homerange_buffer)
  homerange_streams = st_intersection(st_make_valid(watercourses), homerange_buffer)
  homerange_streams_major_length = sum(st_length(homerange_streams_major))
  homerange_streams_length = sum(st_length(homerange_streams))
  mapview(homerange_streams_major) + mapview(homerange_streams) + mapview(homerange_buffer) + homerange_cover
}

