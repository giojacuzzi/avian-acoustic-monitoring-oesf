##############################################################################
# Quantify covariates on occurrence
#
# Note "hs" refers to data collected from in-person habitat surveys, while
# "rs" refers to data derived via remote-sensing imagery.
##############################################################################

# TODO: Ensure that aru site locations are correctly poisitioned within patches (i.e. not near edges so as to get incorrect covariate estimates)

overwrite_rast_cover_cache = FALSE
overwrite_data_plot_scale_cache = TRUE
overwrite_data_homerange_scale_cache = TRUE
path_rast_cover_clean_out = "data/cache/occurrence_covariates/rast_cover_clean.tif"
path_data_plot_scale_out = "data/cache/occurrence_covariates/data_plot_scale.rds"
path_data_homerange_scale_out = "data/cache/occurrence_covariates/data_homerange_scale.rds"

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
study_area = st_buffer(st_as_sfc(st_bbox(aru_sites)), dist = 100) # TODO: Re-run with larger buffer

# Load a raster, cropped, projected, and masked to the study area
load_raster = function(path_rast) {
  r = rast(path_rast)
  sa_r = vect(st_transform(study_area, crs(r)))
  r = crop(r, sa_r)
  r = mask(r, sa_r)
  r = project(r, crs_m_rast)
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

# Base RS-FRIS data (0.1 acre resolution, i.e. ~404m2 or 20.10836 * 20.10836 m grain, roughly 1% of the area of a 100m radius circle)
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
data_plot_scale = aru_sites
data_plot_scale$hs = FALSE # flag for habitat survey data availability

seral_classes = c("stand initiation", "canopy closure", "stem exclusion", "mature/old forest")

# WADNR delineated patch vectors. Manually determined from a combination of aerial imagery and remote sensing.
vect_patches = st_read('data/environment/GIS Data/Forest Development Strata/AgeStrataFixedDIS_RSFRIS20200130.shp') %>%
  st_transform(crs_m) %>% janitor::clean_names() %>% mutate(wadnr_stage_vect = stratum %>% str_to_lower()) %>%
  mutate(stratum = factor(stratum, levels = seral_classes)) %>% select(-stratum)
data_plot_scale = st_join(data_plot_scale, vect_patches)

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
data_plot_scale$stage = terra::extract(rast_stage_class, vect(data_plot_scale))[, 2]
data_plot_scale$stage = factor(data_plot_scale$stage, labels = stage_classes)

# Classify stand seral stage classes
wadnr_classes = c("stand initiation", "canopy closure", "stem exclusion", "mature/old forest") # classify origin raster based on number of years since 2020
rast_wadnr_class = classify(rast_age, matrix(c(
  0,   15,  1,   # "stand initiation"
  15,  25,  2,   # "canopy closure"
  25,  80,  3,   # "stem exclusion"
  80, Inf, 4     # "mature/old forest"
), ncol = 3, byrow = TRUE), include.lowest = TRUE, right = FALSE)
rast_wadnr_class = as.factor(rast_wadnr_class)
unique(terra::values(rast_wadnr_class))
data_plot_scale$wadnr_stage_rast = terra::extract(rast_wadnr_class, vect(data_plot_scale))[, 2]
data_plot_scale$wadnr_stage_rast = factor(data_plot_scale$wadnr_stage_rast, labels = wadnr_classes)

# Determine thinning treatment
poly_thinning_treatment = st_read('data/environment/GIS Data/Forest Development Strata/ThinAfter94NoHarvSinceClipByInitBuf3.shp') %>% 
  st_transform(crs_m) %>% select(TECHNIQUE_, FMA_DT, FMA_STATUS) %>% janitor::clean_names() %>%
  mutate(technique = technique %>% str_to_lower(), fma_status = fma_status %>% str_to_lower()) %>%
  rename(thinning_treatment = technique, thinning_status = fma_status, thinning_date = fma_dt)
data_plot_scale = st_join(data_plot_scale, poly_thinning_treatment)
data_plot_scale = data_plot_scale %>% mutate(stratum = if_else(!is.na(thinning_treatment), "thinned", stage))
strata = c("stand initiation", "stem exclusion", "thinned", "understory reinitiation", "old forest") 
data_plot_scale$stratum = factor(data_plot_scale$stratum, labels = strata)
table(data_plot_scale$thinning_treatment, useNA = 'ifany')
table(data_plot_scale$stratum, useNA = 'ifany')

# TODO: Resolve the small number of sites that are ambiguously classified due to origin age near class boundaries:
# Dz303i > age 24, stand init / canopy closure / WADNR stemex
# Ca263i > age 86, understory reinit, WADNR stemex
# Dp166i > age 24, stand init / canopy closure / WANDR stemex

# From here forward, we use "rast_stage_class" as the classification scheme to delineate patches
# TODO: Also consider instead using a raster with stratum (i.e. classifying thinned stands as their own cover types)
mapview(vect_patches) + mapview(poly_thinning_treatment)
summary(data_plot_scale$stage)
mapview(rast_stage_class) + mapview(aru_sites, label = aru_sites$site) + mapview(st_buffer(aru_sites, 100), col.regions = 'transparent', lwd = 2)

# Further delineate patch boundaries by roads, waterbodies, and watercourses

# Paved roads. Unpaved roads do not constitute a boundary of a patch because they are narrow, rarely driven and likely experienced by the birds as gaps in the forest.
roads_wadnr = st_read('data/environment/GIS Data/roads/T3roads.shp') %>%
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

  message('Saving raster cover data cache ', path_rast_cover_clean_out)
  dir.create(dirname(path_rast_cover_clean_out), recursive = TRUE, showWarnings = FALSE)
  writeRaster(rast_cover_clean, path_rast_cover_clean_out, overwrite=TRUE)
  
} else { # overwrite_rast_cover_cache is FALSE
  message('Loading raster cover data from cache ', path_rast_cover_clean_out)
  rast_cover_clean = rast(path_rast_cover_clean_out)
}
rast_cover_clean = as.factor(rast_cover_clean)

mapview(rast_cover_clean,
        alpha.regions = 1.0,
        col.regions = c('#90c6bd', '#3c8273', '#d8c18a', '#9b652b', 'darkgray', '#6495ED'))

##############################################################################
# Local plot scale covariates

if (overwrite_data_plot_scale_cache) {

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
    mutate(tag = str_extract(site, "_.*$") %>% str_remove("^_"), site = str_remove(site, "_.*$"))
  
  # Coalesce duplicate site entries
  # TODO: Double check this makes sense -- are there any "moved" sites?
  data_plot = data_plot %>% group_by(site) %>% summarise(across(everything(), ~ coalesce(.[!is.na(.)][1], NA)))
  
  # Remove sites with missing data
  # TODO: Re-add these sites when data are available
  data_plot = data_plot %>% filter(!site %in% c("Bz217i", "Bz219i"))
  
  # Mark sites that were surveyed for habitat data in-person
  intersect(data_plot_scale$site, data_plot$site)
  length(intersect(data_plot_scale$site, data_plot$site))
  setdiff(data_plot$site, data_plot_scale$site) # This should be 0
  data_plot_scale[data_plot_scale$site %in% data_plot$site, 'hs'] = TRUE
  mapview(data_plot_scale, zcol = 'hs')
  
  # Elevation [m] (derived from WA state LiDAR flights 3x3-ft raster grid)
  elevation = read.csv('data/environment/site_elevation_ft.csv')
  elevation$plot_elev_rs = elevation$elev_ft * conv_ft_m
  data_plot_scale = data_plot_scale %>% left_join(elevation %>% select(site, plot_elev_rs), by = 'site')
  
  # Age (mean and cv) [#]
  data_plot_scale$age_point = terra::extract(rast_age, vect(data_plot_scale))[, 2]
  data_plot_scale$age_mean = as.numeric(compute_raster_buffer_value_func(rast_age, data_plot_scale, plot_buffer, mean)[,2])
  data_plot_scale$age_cv = as.numeric(compute_raster_buffer_value_func(rast_age, data_plot_scale, plot_buffer, compute_cv)[,2])
  hist(data_plot_scale$age_point, breaks = seq(0, max(data_plot_scale$age_point) + 10, by = 10))
  hist(data_plot_scale$age_mean, breaks = seq(0, max(data_plot_scale$age_mean) + 10, by = 10))
  table(data_plot_scale$age_point)
  table(round(data_plot_scale$age_mean))
  
  ggplot(data_plot_scale, aes(x = stage, y = age_mean)) + geom_boxplot() + theme_minimal()
  
  # Cover class [categorical]
  data_plot_scale$stage
  
  # Thinning treatment [categorical]
  data_plot_scale$thinning_treatment
  
  # Basal area [m2/ha] TODO: confirm if this is a mean value at local level
  data_plot_scale = data_plot_scale %>% left_join(data_plot %>% select(site, plot_ba_hs = ba_ha_all), by = 'site')
  
  rast_ba = load_raster(paste0(dir_rsfris_version, '/RS_FRIS_BA_EXPORT.tif'))
  data_plot_scale$plot_ba_rs = as.numeric(compute_raster_buffer_value_func(rast_ba, data_plot_scale, plot_buffer, mean)[,2]) * conv_sqftAcre_m2Ha
  data_plot_scale$plot_ba_sd_rs = as.numeric(compute_raster_buffer_value_func(rast_ba, data_plot_scale, plot_buffer, sd)[,2]) * conv_sqftAcre_m2Ha
  
  ggplot(data_plot_scale, aes(x = plot_ba_hs, y = plot_ba_rs, label = site)) +
    geom_point() + geom_text_repel(size = 2) + geom_abline(slope = 1) +
    xlim(0, max(data_plot_scale$plot_ba_hs, data_plot_scale$plot_ba_rs, na.rm = TRUE)) +
    ylim(0, max(data_plot_scale$plot_ba_hs, data_plot_scale$plot_ba_rs, na.rm = TRUE)) +
    ggtitle('Basal area [m2/ha]')
  mapview(rast_ba) +
    mapview(data_plot_scale, zcol = 'plot_ba_rs', legend = T)
  
  # Tree density (large, dbh > 10 cm) [# trees/ha]
  data_plot_scale = data_plot_scale %>% left_join(data_plot %>% select(site, plot_treeden_gt10cmDbh_hs = large_per_hectare_all), by = 'site')
  
  rast_treeden_gt4inDbh = load_raster(paste0(dir_rsfris_version, '/RS_FRIS_TREE_ACRE_4.img'))
  data_plot_scale$plot_treeden_gt4inDbh_rs = as.numeric(compute_raster_buffer_value_func(rast_treeden_gt4inDbh, data_plot_scale, plot_buffer, mean)[,2]) * conv_acre_hectare
  
  ggplot(data_plot_scale, aes(x = plot_treeden_gt10cmDbh_hs, y = plot_treeden_gt4inDbh_rs, label = site)) +
    geom_point() + geom_text_repel(size = 2) + geom_abline(slope = 1) +
    xlim(0, max(data_plot_scale$plot_treeden_gt10cmDbh_hs, data_plot_scale$plot_treeden_gt4inDbh_rs, na.rm = TRUE)) +
    ylim(0, max(data_plot_scale$plot_treeden_gt10cmDbh_hs, data_plot_scale$plot_treeden_gt4inDbh_rs, na.rm = TRUE)) +
    ggtitle('Tree density (large, dbh > 10 cm) [# trees/ha]')
  
  # mapview(rast_treeden_gt4inDbh) +
    # mapview(data_plot_scale, zcol = 'plot_treeden_gt4inDbh', legend = T)
  
  # Tree density (small, dbh < 10 cm) [# trees/ha]
  data_plot_scale = data_plot_scale %>% left_join(data_plot %>% select(site, plot_treeden_lt10cmDbh_hs = small_per_hectare), by = 'site')
  
  rast_treeden_all = load_raster(paste0(dir_rsfris_version, '/RS_FRIS_TREE_ACRE.img'))
  data_plot_scale$plot_treeden_all_rs = as.numeric(compute_raster_buffer_value_func(rast_treeden_all, data_plot_scale, plot_buffer, mean)[,2]) * conv_acre_hectare
  data_plot_scale$plot_treeden_lt4inDbh_rs = data_plot_scale$plot_treeden_all_rs - data_plot_scale$plot_treeden_gt4inDbh_rs
  
  ggplot(data_plot_scale, aes(x = plot_treeden_lt10cmDbh_hs, y = plot_treeden_lt4inDbh_rs, label = site)) +
    geom_point() + geom_text_repel(size = 2) + geom_abline(slope = 1) +
    xlim(0, max(data_plot_scale$plot_treeden_lt10cmDbh_hs, data_plot_scale$plot_treeden_lt4inDbh_rs, na.rm = TRUE)) +
    ylim(0, max(data_plot_scale$plot_treeden_lt10cmDbh_hs, data_plot_scale$plot_treeden_lt4inDbh_rs, na.rm = TRUE)) +
    ggtitle('Tree density (small, dbh < 10 cm) [# trees/ha]')
  
  # Total tree density (all sizes) [# trees/ha]
  data_plot_scale$plot_treeden_all_hs = data_plot_scale$plot_treeden_gt10cmDbh_hs + data_plot_scale$plot_treeden_lt10cmDbh_hs
  
  ggplot(data_plot_scale, aes(x = plot_treeden_all_hs, y = plot_treeden_all_rs, label = site)) +
    geom_point() + geom_text_repel(size = 2) + geom_abline(slope = 1) +
  xlim(0, max(data_plot_scale$plot_treeden_all_hs, data_plot_scale$plot_treeden_all_rs, na.rm = TRUE)) +
    ylim(0, max(data_plot_scale$plot_treeden_all_hs, data_plot_scale$plot_treeden_all_rs, na.rm = TRUE)) +
    ggtitle('Tree density (total, all sizes) [# trees/ha]')
  
  # Stand density index (Reineke's) [#]
  rast_sdi = load_raster(paste0(dir_rsfris_version, '/RS_FRIS_SDI_SUM.img'))
  data_plot_scale$plot_sdi_rs = as.numeric(compute_raster_buffer_value_func(rast_sdi, data_plot_scale, plot_buffer, mean)[,2])
  
  # Tree diameter (mean and SD) [cm]
  data_plot_scale = data_plot_scale %>% left_join(
    data_plot %>% select(site, plot_qmd_gt10cmDbh_hs = avg_dbh_cm_all) %>%
      mutate(plot_qmd_gt10cmDbh_hs = replace_na(plot_qmd_gt10cmDbh_hs, 0.0)),
    by = 'site')
  data_plot_scale = data_plot_scale %>% left_join(
    data_plot %>% select(site, plot_qmd_lt10cmDbh_hs = avg_dbh_cm) %>%
      mutate(plot_qmd_lt10cmDbh_hs = replace_na(plot_qmd_lt10cmDbh_hs, 0.0)),
    by = 'site')
  data_plot_scale$plot_qmd_all_hs = data_plot_scale$plot_qmd_gt10cmDbh_hs + data_plot_scale$plot_qmd_lt10cmDbh_hs
  
  rast_qmd = load_raster(paste0(dir_rsfris_version, '/RS_FRIS_QMD.img')) #load_raster('data/environment/rs_fris/rs_fris_QMD/RS_FRIS_QMD.tif')
  data_plot_scale$plot_qmd_rs = as.numeric(compute_raster_buffer_value_func(rast_qmd, data_plot_scale, plot_buffer, mean)[,2]) * conv_in_cm
  
  ggplot(data_plot_scale, aes(x = plot_qmd_all_hs, y = plot_qmd_rs, label = site)) +
    geom_point() + geom_text_repel(size = 2) + geom_abline(slope = 1) +
    xlim(0, max(data_plot_scale$plot_qmd_all_hs, data_plot_scale$plot_qmd_rs, na.rm = TRUE)) +
    ylim(0, max(data_plot_scale$plot_qmd_all_hs, data_plot_scale$plot_qmd_rs, na.rm = TRUE)) +
    ggtitle('Tree diameter (quadratic mean) [cm]')
  mapview(rast_qmd) +
    mapview(data_plot_scale, zcol = 'plot_qmd_rs', legend = T)
  
  # Tree height (mean and cv) [m] and canopy layers [#]
  # TODO: Get high-resolution LiDAR height data from DNR approval procedure
  data_plot_scale = data_plot_scale %>% left_join(
    data_plot %>% select(site, plot_ht_hs = avg_height_m_all) %>%
      mutate(plot_ht_hs = replace_na(plot_ht_hs, 0.0)),
    by = 'site')
  data_plot_scale = data_plot_scale %>% left_join(
    data_plot %>% select(site, plot_ht_cv_hs = cv_height_all) %>%
      mutate(plot_ht_cv_hs = replace_na(plot_ht_cv_hs, 0.0)),
    by = 'site')
  
  rast_htmax = load_raster(paste0(dir_rsfris_version, '/RS_FRIS_HTMAX.img')) # NOTE: feet
  max_ht_m = max(values(rast_htmax), na.rm = TRUE) * conv_ft_m
  data_plot_scale$plot_htmax_rs = as.numeric(compute_raster_buffer_value_func(rast_htmax, data_plot_scale, plot_buffer, mean)[,2]) * conv_ft_m
  data_plot_scale$plot_htmax_cv_rs = as.numeric(compute_raster_buffer_value_func(rast_htmax, data_plot_scale, plot_buffer, compute_cv)[,2]) * conv_ft_m * 100
  
  ggplot(data_plot_scale, aes(x = plot_ht_hs, y = plot_htmax_rs, label = site)) +
    geom_point() + geom_text_repel(size = 2) + geom_abline(slope = 1) +
    xlim(0, max(data_plot_scale$plot_ht_hs, data_plot_scale$plot_htmax_rs, na.rm = TRUE)) +
    ylim(0, max(data_plot_scale$plot_ht_hs, data_plot_scale$plot_htmax_rs, na.rm = TRUE)) +
    ggtitle('Tree height (mean) [m]')
  
  ggplot(data_plot_scale, aes(x = plot_ht_cv_hs, y = plot_htmax_cv_rs, label = site)) +
    geom_point() + geom_text_repel(size = 2) + geom_abline(slope = 1) +
    xlim(0, max(data_plot_scale$plot_ht_cv_hs, data_plot_scale$plot_htmax_cv_rs, na.rm = TRUE)) +
    ylim(0, max(data_plot_scale$plot_ht_cv_hs, data_plot_scale$plot_htmax_cv_rs, na.rm = TRUE)) +
    ggtitle('Tree height (cv) [m]')
  
  rast_canopy_layers = load_raster(paste0(dir_rsfris_version, '/RS_FRIS_CANOPY_LAYERS.img'))
  data_plot_scale$canopy_layers_rs = as.numeric(compute_raster_buffer_value_func(rast_canopy_layers, data_plot_scale, plot_buffer, mean)[,2])
  
  # Canopy cover and closure [%]
  rast_canopy_cover = load_raster(paste0(dir_rsfris_version, '/RS_FRIS_COVER.img'))
  data_plot_scale$plot_canopy_cover_rs = as.numeric(compute_raster_buffer_value_func(rast_canopy_cover, data_plot_scale, plot_buffer, mean)[,2])
  
  rast_canopy_closure = load_raster(paste0(dir_rsfris_version, '/RS_FRIS_CLOSURE.img'))
  data_plot_scale$plot_canopy_closure_rs = as.numeric(compute_raster_buffer_value_func(rast_canopy_closure, data_plot_scale, plot_buffer, mean)[,2])
  
  mapview(rast_canopy_cover) +
    mapview(data_plot_scale, zcol = 'plot_canopy_cover_rs', legend = T)
  
  # Height to live crown [m], length of live crown (live crown depth) [m], and live crown ratio [#]
  data_plot_scale = data_plot_scale %>% left_join(
    data_plot %>% select(site, plot_hlc_hs = avg_hlc_all) %>%
      mutate(plot_hlc_hs = replace_na(plot_hlc_hs, 0.0)), by = 'site')
  data_plot_scale = data_plot_scale %>% left_join(
    data_plot %>% select(site, plot_llc_hs = avg_llc_all) %>%
      mutate(plot_llc_hs = replace_na(plot_llc_hs, 0.0)), by = 'site')
  data_plot_scale = data_plot_scale %>% left_join(
    data_plot %>% select(site, plot_lcr_hs = avg_lcr_all) %>%
      mutate(plot_lcr_hs = replace_na(plot_lcr_hs, 0.0)), by = 'site')
  
  # Density of all snags [# snags/ha]
  data_plot_scale = data_plot_scale %>% left_join(data_plot %>% select(site, plot_snagden_hs = snags_ha) %>% mutate(plot_snagden_hs = replace_na(plot_snagden_hs, 0.0)), by = 'site')
  
  rast_snagden = load_raster(paste0(dir_rsfris_version, '/RS_FRIS_SNAG_ACRE_15.img'))
  data_plot_scale$plot_snagden_rs = as.numeric(compute_raster_buffer_value_func(rast_snagden, data_plot_scale, plot_buffer, mean)[,2]) * conv_acre_hectare
  
  ggplot(data_plot_scale, aes(x = plot_snagden_hs, y = plot_snagden_rs, label = site)) +
    geom_point() + geom_text_repel(size = 2) + geom_abline(slope = 1) +
    xlim(0, max(data_plot_scale$plot_snagden_hs, data_plot_scale$plot_snagden_rs, na.rm = TRUE)) +
    ylim(0, max(data_plot_scale$plot_snagden_hs, data_plot_scale$plot_snagden_rs, na.rm = TRUE)) +
    ggtitle('Snag density [#/ha]')
  
  # Downed wood volume [m3/ha]
  data_plot_scale = data_plot_scale %>% left_join(data_plot %>% select(site, plot_downvol_hs = vol_alldown_m3) %>% mutate(plot_downvol_hs = replace_na(plot_downvol_hs, 0.0)), by = 'site')
  
  rast_downvol = load_raster(paste0(dir_rsfris_version, '/RS_FRIS_CFVOL_DDWM.img'))
  data_plot_scale$plot_downvol_rs = as.numeric(compute_raster_buffer_value_func(rast_downvol, data_plot_scale, plot_buffer, mean)[,2]) * (conv_ft3_m3 / conv_hectare_acre)
  
  ggplot(data_plot_scale, aes(x = plot_downvol_hs, y = plot_downvol_rs, label = site)) +
    geom_point() + geom_text_repel(size = 2) + geom_abline(slope = 1) +
    xlim(0, max(data_plot_scale$plot_downvol_hs, data_plot_scale$plot_downvol_rs, na.rm = TRUE)) +
    ylim(0, max(data_plot_scale$plot_downvol_hs, data_plot_scale$plot_downvol_rs, na.rm = TRUE)) +
    ggtitle('Downed wood volume [m3/ha]')
  
  # Understory vegetation cover [%] and volume [m3/ha]
  data_plot_scale = data_plot_scale %>% left_join(data_plot %>% select(site, plot_understory_cover = per_cover_total_understory) %>% mutate(plot_understory_cover = replace_na(plot_understory_cover, 0.0)), by = 'site')
  data_plot_scale = data_plot_scale %>% left_join(data_plot %>% select(site, plot_understory_vol = shrub_layer_vol) %>% mutate(plot_understory_vol = replace_na(plot_understory_vol, 0.0)), by = 'site')
  
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
    tree_gte10cm_density_psme = tree_gte10cm_obs$large_per_hectare_psme,
    tree_gte10cm_density_thpl = tree_gte10cm_obs$large_per_hectare_thpl,
    tree_gte10cm_density_abam = tree_gte10cm_obs$large_per_hectare_abam,
    tree_gte10cm_density_tshe = tree_gte10cm_obs$large_per_hectare_tshe,
    tree_gte10cm_density_alru = tree_gte10cm_obs$large_per_hectare_alru,
    tree_gte10cm_density_pisi = tree_gte10cm_obs$large_per_hectare_pisi,
    tree_all_density_psme = tree_all_obs$all_per_hectare_psme,
    tree_all_density_thpl = tree_all_obs$all_per_hectare_thpl,
    tree_all_density_abam = tree_all_obs$all_per_hectare_abam,
    tree_all_density_tshe = tree_all_obs$all_per_hectare_tshe,
    tree_all_density_alru = tree_all_obs$all_per_hectare_alru,
    tree_all_density_pisi = tree_all_obs$all_per_hectare_pisi
  )
  
  data_plot_scale = data_plot_scale %>% left_join(tree_div_metrics, by = 'site')
  
  ggplot(data_plot_scale, aes(x = `stratum`, y = tree_all_diversity)) +
    geom_violin() +
    stat_summary(fun = mean, geom = "point") +
    ggtitle("Tree shannon diversity across strata") + theme_minimal()
  
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
  df_long <- data_plot_scale %>% select(site, `stage`, starts_with(metric)) %>%
    pivot_longer(cols = starts_with(metric), names_to = "species", values_to = "density") %>%
    mutate(species = gsub(metric, "", species))
  ggplot(df_long, aes(x = species, y = density, fill = species)) +
    geom_boxplot() +
    facet_wrap(~ stage, scales = "free_y") +
    coord_cartesian(ylim = c(0, 1250)) + 
    labs(title = "Tree species density by stage") +
    theme_minimal()
  
  # Inspect specific variables by stage
  ggplot(data_plot_scale, aes(x = stage, y = canopy_layers_rs)) +
    geom_boxplot() +
    theme_minimal()
  
  # Distance to water (type 1-3 and all types) [m]
  dist_watercourses_major = st_distance(data_plot_scale, boundary_watercourses)
  dist_watercourses_major = apply(dist_watercourses_major, 1, min)
  data_plot_scale$dist_watercourses_major = dist_watercourses_major
  mapview(data_plot_scale, zcol = "dist_watercourses_major") + mapview(boundary_watercourses)
  
  dist_watercourses_all = st_distance(data_plot_scale, watercourses)
  dist_watercourses_all = apply(dist_watercourses_all, 1, min)
  data_plot_scale$dist_watercourses_all = dist_watercourses_all
  mapview(data_plot_scale, zcol = "dist_watercourses_all") + mapview(watercourses)
  
  # Distance to nearest edge [m]
  # For each site, find its patch, then find the distance to the nearest non-patch cell
  patch_ids <- patches(rast_cover_clean, directions=8, values=TRUE)
  data_plot_scale$dist_nearest_edge = NA
  for (i in 1:nrow(data_plot_scale)) {
    d = data_plot_scale[i,]
    cover_class = terra::extract(rast_cover_clean, vect(d))[,2]
    pid = terra::extract(patch_ids, vect(d))[,2]
    m = trim(classify(patch_ids, cbind(pid, 1), others = NA))
    na_cells <- which(is.na(values(m)))
    na_coords <- xyFromCell(m, na_cells)
    d_coords <- st_coordinates(d)
    distances <- sqrt((na_coords[,1] - d_coords[1])^2 + (na_coords[,2] - d_coords[2])^2)
    min_dist <- min(distances)
    data_plot_scale[i, 'dist_nearest_edge'] = min_dist
    # mapview(d) + mapview(m) + mapview(st_buffer(d, 100), col.regions = 'transparent', lwd = 2)
  }
  mapview(rast_cover_clean) + mapview(data_plot_scale, label = data_plot_scale$site, zcol = 'dist_nearest_edge') + mapview(st_buffer(data_plot_scale, 100), col.regions = 'transparent', lwd = 2)
  
  message('Saving plot scale data cache ', path_data_plot_scale_out)
  dir.create(dirname(path_data_plot_scale_out), recursive = TRUE, showWarnings = FALSE)
  saveRDS(data_plot_scale, path_data_plot_scale_out)

} else { # overwrite_data_plot_scale_cache is FALSE
  message('Loading plot scale data from cache ', path_data_plot_scale_out)
  data_plot_scale = readRDS(path_data_plot_scale_out)
}

##############################################################################
# Focal patch scale and home range neighborhood scale covariates (limited to the extent of the species home range)

if (overwrite_data_homerange_scale_cache) {

  data_homerange_scale = list()
  
  # For now, just calculate at a fixed radius plot scale (e.g. median predicted radius of 200 m)
  # TODO: Calculate per species according to predicted home range size
  homeranges = data.frame(
    scale = c('Plot', 'Median', 'Mean', 'Nearmax'),
    size    = c(100, # plot scale
                204, # median
                468, # mean
                1000) # meters
  )
  for (i in 1:nrow(homeranges)) {
    scale = homeranges[i,'scale']
    homerange_buffer_size = homeranges[i,'size']
    message("Scale ", scale, ", home range buffer size ", homerange_buffer_size)
    
    data_homerange_scale_species = data.frame()
  
    # "Edge influences on distributions of organisms or factors affecting organisms (such as predation and nest parasitism) are concentrated within 50m of the edge." (Kremaster L. & Bunnell F. L. (1999) Edge effects: Theory, evidence and implications to management of western North American forests. In: Forest Fragmentation: Wildlife and Management Implications (eds J. Wisniewski, J. A. Rochelle & L. Lehmann) pp. 117–53. Leiden, Boston, MA.)
    core_area_buffer = 50 # meters
    core_area_buffer_ncells = round(core_area_buffer / res(rast_cover_clean)[1], 0)
    
    sites = 1:nrow(data_plot_scale)
    for (j in sites) {
      site = data_plot_scale[j,]
      message(site$site, ' (', j / length(sites), ')')
      cover_class = terra::extract(rast_cover_clean, vect(site))[,2]
      
      homerange_and_edge_buffer = st_buffer(site, homerange_buffer_size + core_area_buffer) # additional buffer to ensure that edges are retained in mask
      homerange_buffer = st_buffer(site, homerange_buffer_size)
      homerange_and_edge_crop = crop(rast_cover_clean, vect(homerange_and_edge_buffer))
      homerange_and_edge = mask(homerange_and_edge_crop, vect(homerange_and_edge_buffer))
      
      patch_ids = patches(homerange_and_edge, directions=8, values=TRUE)
      pid = terra::extract(patch_ids, vect(site))[,2]
      focal_patch = trim(classify(patch_ids, cbind(pid, 1), others = NA))
      focal_patch = crop(focal_patch, vect(homerange_buffer))
      focal_patch = mask(focal_patch, vect(homerange_buffer))
      
      # Focal patch area [ha]
      # lsm_p_area(focal_patch, directions = 8)[1, 'value'] # in hectares (ha)
      # global(!is.na(focal_patch), sum, na.rm = TRUE) * prod(res(focal_patch)) # in sq meters
      # mapview(site) + mapview(homerange_and_edge) + mapview(homerange_buffer, col.regions = 'transparent', lwd = 2) + mapview(focal_patch)
      
      ncells_focal_patch = sum(!is.na(values(focal_patch)))
      (focal_patch_area = ncells_focal_patch * prod(res(focal_patch)) * conv_m2_ha)
      
      # Focal patch core area [ha]
      focal_patch_and_edge_buffer = trim(classify(patch_ids, cbind(pid, 1), others = NA))
      focal_patch_inv = ifel(is.na(focal_patch_and_edge_buffer), 1, NA)
      focal_patch_edge_dist = distance(focal_patch_inv)
      focal_patch_edge_dist = ifel(focal_patch_edge_dist == 0, NA, focal_patch_edge_dist)
      focal_patch_core = ifel(focal_patch_edge_dist >= core_area_buffer, 1, NA)
      focal_patch_core = crop(focal_patch_core, vect(homerange_buffer))
      focal_patch_core = mask(focal_patch_core, vect(homerange_buffer))
      # mapview(homerange_buffer, col.regions = 'transparent', lwd = 2) + mapview(focal_patch_edge_dist) + mapview(focal_patch) + mapview(focal_patch_and_edge_buffer) + mapview(focal_patch_core)
      
      ncells_focal_patch_core = sum(!is.na(values(focal_patch_core)))
      (focal_patch_core_area = ncells_focal_patch_core * prod(res(focal_patch)) * conv_m2_ha)
      
      homerange_crop = crop(rast_cover_clean, vect(homerange_buffer))
      homerange_cover = mask(homerange_crop, vect(homerange_buffer))
      homerange_cover_forest = homerange_cover
      homerange_cover_forest[!(homerange_cover_forest[] %in% c(1, 2, 3, 4))] <- NA
      # mapview(homerange_cover)
      ncells_homerange = sum(!is.na(values(homerange_cover)))
      
      # Focal patch percentage of home range [%]
      (focal_patch_pcnt = sum(!is.na(values(focal_patch))) / ncells_homerange)
      
      # Focal patch core percentage of home range [%]
      (focal_patch_core_pcnt = sum(!is.na(values(focal_patch_core))) / ncells_homerange)
  
      # Focal patch euclidean nearest neighbor distance [m]
      # Quantifies habitat isolation
      # Approaches 0 as the distance to the nearest neighbor decreases. Minimum is constrained by the cell size.
      focal_cover_class = as.numeric(site$stage)
      matching_cover = rast_cover_clean
      matching_cover[matching_cover != focal_cover_class] = NA
      focal_patch_extended = extend(focal_patch, matching_cover)
      matching_neighbors = mask(matching_cover, focal_patch_extended, maskvalue=NA, inverse=TRUE)
      focal_patch_distance = distance(focal_patch_extended)
      matching_neighbor_distance = mask(focal_patch_distance, matching_neighbors, maskvalue=NA, inverse=FALSE)
      # mapview(matching_neighbors) + mapview(focal_patch) + mapview(focal_patch_distance) + mapview(matching_neighbor_distance)
      (focal_patch_isolation = min(values(matching_neighbor_distance), na.rm = TRUE))
      
      # Cover class richness [#]
      # Quantifies the number of cover types present in home range
      cover_freq = freq(homerange_cover) %>% select(value, count) %>% mutate(value = as.integer(value))
      cover_freq = data.frame(value = 1:6) %>%
        left_join(cover_freq, by = "value") %>%
        mutate(count = ifelse(is.na(count), 0, count))
      cover_forest_freq = cover_freq %>% filter(value %in% c(1,2,3,4))
      (cover_richness = cover_freq %>% filter(count != 0) %>% nrow())
      (cover_forest_richness = cover_forest_freq %>% filter(count != 0) %>% nrow())
      
      # Cover class evenness (i.e. dominance) [#]
      # Quantifies the degree of evenness versus dominance in cover type distribution 
      # 0 when only one cover type is present, 1 when types are equally distributed
      cover_evenness = lsm_l_shei(homerange_cover) %>% pull(value)
      cover_forest_evenness = lsm_l_shei(homerange_cover_forest) %>% pull(value)
      
      # Cover class diversity (shannon) [#]
      # Quantifies both the richness and evenness of cover type distributions (inversely related to contagion)
      # 0 when only one cover type is present and increases as the number of classes increases while the proportions are equally distributed
      cover_diversity = lsm_l_shdi(homerange_cover) %>% pull(value)
      cover_forest_diversity = lsm_l_shdi(homerange_cover_forest) %>% pull(value)
      
      # Proportional abundance of each cover class [%]
      (prop_abund_1 = cover_freq %>% filter(value == 1) %>% pull(count) / ncells_homerange)
      (prop_abund_2 = cover_freq %>% filter(value == 2) %>% pull(count) / ncells_homerange)
      (prop_abund_3 = cover_freq %>% filter(value == 3) %>% pull(count) / ncells_homerange)
      (prop_abund_4 = cover_freq %>% filter(value == 4) %>% pull(count) / ncells_homerange)
      (prop_abund_5 = cover_freq %>% filter(value == 5) %>% pull(count) / ncells_homerange)
      (prop_abund_6 = cover_freq %>% filter(value == 6) %>% pull(count) / ncells_homerange)
      
      # Aggregation index [#]
      # Quantifies degree of habitat contiguity versus fragmentation
      # Equals 0 when the patch types are maximally disaggregated (i.e., when there are no like adjacencies); AI increases as the landscape is increasingly aggregated and equals 100 when the landscape consists of a single patch.
      aggregation_idx = lsm_l_ai(homerange_cover, directions = 8) %>% pull(value)
      
      # Shape index [#]
      # Quantifies patch shape complexity
      # Equals 1 if all patches are squares. Increases, without limit, as the shapes of patches become more complex.
      shape_idx = lsm_l_shape_mn(homerange_cover, directions = 8) %>% pull(value)
      
      # Contrast-weighted edge density [m/ha]
      # The density of patch edges weighted by their contrast
      # Equals 0 when there is no edge in the landscape (i.e. landscape consists of a single patch). Increases as the amount of edge in the landscape increases and/or as the contrast in edges increase (i.e. contrast weight approaches 1).
      homerange_height = mask(crop(rast_htmax, homerange_and_edge_buffer), homerange_and_edge_buffer)
      # mapview(patch_ids) + mapview(homerange_height)
      
      height_aligned <- resample(homerange_height, patch_ids, method = "bilinear")
      zonal_means <- zonal(height_aligned, patch_ids, fun = "median", na.rm = TRUE) # TODO: consider min instead
      
      if (nrow(zonal_means) > 1) {
        # Equals the sum of the lengths (m) of each edge segment in the landscape multiplied by the
        # corresponding contrast weight, divided by the total landscape area (m2), converted to hectares.
        colnames(zonal_means) = c('class', 'height')
        # mapview(classify(patch_ids, rcl = zonal_means))
        edge_lengths = get_adjacencies(patch_ids, neighbourhood = 8, what = "unlike", upper = FALSE)[[1]] * res(patch_ids)[1]
        heights <- zonal_means$height
        names(heights) <- zonal_means$class
        height_diff_matrix <- outer(heights, heights, FUN = function(x, y) abs(x - y))
        
        # Weights are derived from proportion of difference in height relative to the maximum potential difference in height (Hou and Walz 2016, Huang et al. 2014)
        # d = 0 --> 0 m difference
        # d = 1 --> max m difference (i.e. maximum tree height of entire landscape)
        weights = height_diff_matrix / max_ht_m
        homerange_rast_area = ncells_homerange * res(homerange_cover)[1] * res(homerange_cover)[2] * conv_m2_ha
        (cw_edge_density = sum(edge_lengths * weights, na.rm = TRUE) / homerange_rast_area)
      } else { 
        # Only one patch type
        (cw_edge_density = 0.0)
      }
      cw_edge_density = set_units(cw_edge_density, 'm/ha')
      
      homerange_buffer_area_ha = set_units(st_area(homerange_buffer), ha)
      
      # Density of roads (paved and all) [m/ha]
      homerange_roads_paved = st_intersection(st_geometry(st_make_valid(roads %>% filter(road_usgs1 %in% c("Primary Highway", "Light-Duty Road")))), homerange_buffer)
      homerange_roads = st_intersection(st_geometry(st_make_valid(roads)), homerange_buffer)
      homerange_roads_paved_length = sum(st_length(homerange_roads_paved))
      homerange_roads_length = sum(st_length(homerange_roads))
      (density_roads_paved = homerange_roads_paved_length / homerange_buffer_area_ha)
      (density_roads = homerange_roads_length / homerange_buffer_area_ha)
      # mapview(homerange_roads_paved) + mapview(homerange_roads) + mapview(homerange_buffer) + homerange_cover
      
      # Density of streams (type 1-3 and all types) [m/ha]
      homerange_streams_major = st_intersection(st_geometry(st_make_valid(watercourses %>% filter(sl_wtrty_c %in% c(1, 2, 3)))), homerange_buffer)
      homerange_streams = st_intersection(st_geometry(st_make_valid(watercourses)), homerange_buffer)
      homerange_streams_major_length = sum(st_length(homerange_streams_major))
      homerange_streams_length = sum(st_length(homerange_streams))
      (density_streams_major = homerange_streams_major_length / homerange_buffer_area_ha)
      (density_streams = homerange_streams_length / homerange_buffer_area_ha)
      # mapview(homerange_streams_major) + mapview(homerange_streams) + mapview(homerange_buffer) + homerange_cover
      
      data_homerange_scale_species = rbind(data_homerange_scale_species, data.frame(
        site = site$site,
        buffer = homerange_buffer_size,
        focal_patch_pcnt,
        focal_patch_core_pcnt,
        focal_patch_isolation,
        cover_richness,
        cover_forest_richness,
        cover_evenness,
        cover_forest_evenness,
        cover_diversity,
        cover_forest_diversity,
        prop_abund_1,
        prop_abund_2,
        prop_abund_3,
        prop_abund_4,
        prop_abund_5,
        prop_abund_6,
        aggregation_idx,
        shape_idx,
        cw_edge_density,
        density_roads_paved,
        density_roads,
        density_streams_major,
        density_streams
      ))
    }
    
    # Sanity check results
    if (FALSE) {
      site_check = data_homerange_scale_species %>% slice_max(focal_patch_isolation, n = 1, with_ties = FALSE) %>% pull(site)
      site_check = data_plot_scale %>% filter(site == site_check)
      site_check_buffer = st_buffer(site_check, homerange_buffer_size)
      mapview(site_check) + mapview(mask(crop(rast_cover_clean, site_check_buffer), site_check_buffer))
      
      site_check = data_homerange_scale_species %>% slice_min(focal_patch_isolation, n = 1, with_ties = FALSE) %>% pull(site)
      site_check = data_plot_scale %>% filter(site == site_check)
      site_check_buffer = st_buffer(site_check, homerange_buffer_size)
      mapview(site_check) + mapview(mask(crop(rast_cover_clean, site_check_buffer), site_check_buffer))
    }
    
    data_homerange_scale[[scale]] = data_homerange_scale_species
  }
  message('Saving homerange data cache ', path_data_homerange_scale_out)
  dir.create(dirname(path_data_homerange_scale_out), recursive = TRUE, showWarnings = FALSE)
  saveRDS(data_homerange_scale, path_data_homerange_scale_out)
  
} else { # overwrite_data_homerange_scale_cache is FALSE
  message('Loading homerange scale data from cache ', path_data_homerange_scale_out)
  data_homerange_scale = readRDS(path_data_homerange_scale_out)
}

#############################################################################################################
# Data inspection

### Plot scale candidate variables

data_plot_scale_candidates = data_plot_scale %>% st_drop_geometry() %>% select(where(is.numeric))

cor_matrix_plot_scale = cor(data_plot_scale_candidates, use = "pairwise.complete.obs", method = "pearson")
cor_matrix_plot_scale[lower.tri(cor_matrix_plot_scale, diag = TRUE)] = NA
collinearity_candidates = subset(as.data.frame(as.table(cor_matrix_plot_scale)), !is.na(Freq) & abs(Freq) >= 0.8)

# Reduce list of candidate variables (highly correlated, less preferable, irrelevant, etc.)
vars_to_drop = c(
  # Highly correlated with age_mean, which is a better representation of stand age
  'age_point',
  # Highly correlated with plot_treeden_all_rs
  'plot_treeden_lt4inDbh_rs',
  # Individual tree species composition is better summarized with taxonomic diversity metrics
  'tree_gte10cm_density_psme', 'tree_gte10cm_density_thpl', 'tree_gte10cm_density_abam', 'tree_gte10cm_density_tshe', 'tree_gte10cm_density_alru', 'tree_gte10cm_density_pisi',
  'tree_all_density_thpl', 'tree_all_density_abam', 'tree_all_density_tshe', 'tree_all_density_alru', 'tree_all_density_pisi', 'tree_all_density_psme',
  # TODO: Basal area (ba), height, qmd, and stand density are highly correlated.
  # Josh and Teddy recommend independent variables for 1) tree density and 2) tree size (QMD).
  # Dan, however, notes that integrative metrics that incorporate both tree density and (mean) tree size do the best job in these analyses.
  # Drop SDI (which is highly correlated with qmd, ba, and ht) in favor of a raw density measurement.
 'plot_sdi_rs',
 # Canopy cover and closure are highly correlated. Drop closure in favor of cover here because it is likely a more accurate measurement and is more widely used in landscape-scale analyses.
 'plot_canopy_closure_rs'
)
data_plot_scale_candidates = data_plot_scale_candidates %>% select(-all_of(vars_to_drop))

corrplot::corrplot(cor_matrix_plot_scale, method = "color", tl.cex = 0.8)






ggplot(data_plot_scale, aes(x = `stage`, y = plot_qmd_all_hs)) +
  geom_violin() +
  stat_summary(fun = mean, geom = "point") +
  ggtitle("Tree species richness across stages") + theme_minimal()

df_long <- data_plot_scale %>% select(site, `stage`, plot_qmd_lt10cmDbh_hs, plot_qmd_gt10cmDbh_hs, plot_qmd_all_hs) %>%
  pivot_longer(cols = starts_with("plot_qmd"), names_to = "metric", values_to = "qmd")
ggplot(df_long, aes(x = stage, y = qmd, fill = metric)) +
  geom_boxplot() + theme_minimal()

ggplot(data_plot_scale %>% select(site, `stage`, plot_treeden_lt10cmDbh_hs, plot_treeden_gt10cmDbh_hs, plot_treeden_all_hs) %>%
         pivot_longer(cols = starts_with("plot_treeden"), names_to = "variable", values_to = "value"),
       aes(x = stage, y = value, fill = variable)) +
  geom_boxplot() + coord_cartesian(ylim = c(0, 2500)) + theme_minimal() +
  ggtitle('Tree density (habitat surveys)')

ggplot(data_plot_scale %>% select(site, `stage`, plot_treeden_lt4inDbh_rs, plot_treeden_gt4inDbh_rs, plot_treeden_all_rs) %>%
         pivot_longer(cols = starts_with("plot_treeden"), names_to = "variable", values_to = "value"),
       aes(x = stage, y = value, fill = variable)) +
  geom_boxplot() + theme_minimal() +
  ggtitle('Tree density (remote sensing)')

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
df_long <- data_plot_scale %>% select(site, age_mean, `stage`, starts_with(metric)) %>%
  pivot_longer(cols = starts_with(metric), names_to = "species", values_to = "density") %>%
  mutate(species = gsub(metric, "", species))
ggplot(df_long, aes(x = species, y = log(density), fill = species)) +
  geom_boxplot() +
  facet_wrap(~ stage, scales = "free_y") +
  # coord_cartesian(ylim = c(0, 1250)) + 
  labs(title = "Tree species density by stage") +
  theme_minimal()

ggplot(data_plot_scale, aes(x = `stage`, y = plot_elev_rs, fill = `stage`)) +
  geom_boxplot()

# Inspect specific variables by stage
ggplot(data_plot_scale, aes(x = stage, y = canopy_layers_rs)) +
  geom_boxplot() +
  theme_minimal()







### Homerange scales
hist(data_homerange_scale[["Nearmax"]]$focal_patch_pcnt)
hist(data_homerange_scale[["Nearmax"]]$cover_forest_richness)
hist(data_homerange_scale[["Plot"]]$cw_edge_density)

