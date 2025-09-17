##############################################################################
# Quantify covariates on occurrence
#
# Note "hs" refers to data collected from in-person habitat surveys, while
# "rs" refers to data derived via remote-sensing imagery.
##############################################################################

# TODO: Ensure that aru site locations are correctly positioned within patches (i.e. not near edges so as to get incorrect covariate estimates)

overwrite_rast_cover_cache = FALSE
overwrite_data_plot_scale_cache = FALSE
overwrite_data_homerange_scale_cache = FALSE
path_rast_cover_clean_out = "data/cache/occurrence_covariates/rast_cover_clean.tif"
path_data_plot_scale_out = "data/cache/occurrence_covariates/data_plot_scale.rds"
path_data_homerange_scale_out = "data/cache/occurrence_covariates/data_homerange_scale.rds"

library(progress)
library(tidyverse)
library(sf)
library(terra)
library(mapview)
options(mapview.maxpixels = 2117676)
library(ggrepel)
theme_set(theme_minimal())
library(landscapemetrics)
library(units)
library(viridis)
library(vegan)

##############################################################################
# Study area, sites, boundaries, and helper functions

crs_m = 32610 # coordinate reference system EPSG:32610, UTM Zone 10N, (meters)
crs_m_rast = "EPSG:32610"

species_trait_data = read.csv("data/cache/species_traits/species_traits.csv")

# ARU locations (sites)
aru_sites = st_read('data/environment/GIS Data/AcousticStations.shp') %>%
  st_drop_geometry() %>% janitor::clean_names() %>% as.data.frame() %>%
  st_as_sf(coords = c("utm_e", "utm_n"), crs = crs_m) %>%
  select(name, ces, treatment, geometry) %>% rename(site = name)
mapview(aru_sites, label = aru_sites$name)

# Study area boundary buffer
study_area = st_buffer(st_as_sfc(st_bbox(aru_sites)), dist = 10000)

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
  region = st_buffer(points, dist = buffer_distance)
  buffer_values = terra::extract(r, vect(region), fun = func, na.rm = TRUE)
}

# Compute summary statistics for the given numeric vector
# Column 1: mean
# Column 2: standard deviation
# Column 3: coefficient of variation
summary_stats = function(x, na.rm = TRUE, conversion_factor = 1) {
  if (na.rm) x = x[!is.na(x)]
  if (length(x) == 0) return(c(mean = NA, sd = NA, cv = NA))
  x = x * conversion_factor
  mu = mean(x)
  sigma = sd(x)
  if (mu == 0) {
    cv = NA
  } else {
    cv = sigma / mu
  }
  return(data.frame(mean = mu, sd = sigma, cv = cv))
}

# Conversion factors
conv_in_to_cm =                  2.54 # inches to centimeters
conv_ft_to_m =                 0.3048 # feet to meters
conv_ft2peracre_to_m2perha = 0.229568 # square feet/acre to square meters/ha
conv_ft3_to_m3 = 0.0283168            # cubic feet to cubic meters
conv_ft3peracre_to_m3perha = 0.069968 # cubic feet/acre to cubic meters/hectare
conv_peracre_to_perha = 2.47105       # units per acreto units per hectare
conv_perha_to_peracre = 0.404686      # units per hectare to units per acre
conv_m2_to_ha         = 0.0001        # square meters to ha # TODO: CHECK THIS

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
# TODO: Impute patches identified on Desktop pdf from canopy cover, size class, canopy layers, and surrounding classes (see Powell vegetation stages white paper). Also consult ESRI World Imagery Wayback for visual inspection.
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
data_plot_scale$wadnr_stage_rast = factor(data_plot_scale$wadnr_stage_rast)

# Fill in missing patches with WADNR predictions
rast_impute = mask(rasterize(vect(vect_patches), rast_origin_missing, field = "wadnr_stage_vect"), rast_origin_missing)
rast_impute = classify(rast_impute, matrix(c(
  1, 2, 1,    # "initiation" > "stand initiation"
  3, 4, 2,    # "stemex" > "stem exclusion"
  0, 1, 2,    # "canclose" > "stem exclusion"
  2, 3, 3     # "understory reinitiation"
), ncol = 3, byrow = TRUE), include.lowest = TRUE, right = FALSE)
mapview(rast_stage_class) + mapview(rast_impute)
rast_stage_class = cover(rast_impute, rast_stage_class)

# Determine thinning treatment
poly_thinning_treatment = st_read('data/environment/GIS Data/Forest Development Strata/ThinAfter94NoHarvSinceClipByInitBuf3.shp') %>% 
  st_transform(crs_m) %>% select(TECHNIQUE_, FMA_DT, FMA_STATUS) %>% janitor::clean_names() %>%
  mutate(technique = technique %>% str_to_lower(), fma_status = fma_status %>% str_to_lower()) %>%
  rename(thinning_treatment = technique, thinning_status = fma_status, thinning_date = fma_dt)
sum(poly_thinning_treatment$thinning_status == "completed") == nrow(poly_thinning_treatment) # check that all treatments are completed
commrcl_thin = poly_thinning_treatment %>% filter(thinning_treatment == 'commrcl_thin') %>% st_union()
variabl_thin = poly_thinning_treatment %>% filter(thinning_treatment == 'variabl_thin') %>% st_union()
commrcl_sf = st_sf(
  thinning_status = TRUE,
  thinning_treatment = "commerical_thin",
  geometry = poly_thinning_treatment %>% filter(thinning_treatment == 'commrcl_thin') %>% st_union()
)
variabl_sf = st_sf(
  thinning_status = TRUE,
  thinning_treatment = "variable_thin",
  geometry = poly_thinning_treatment %>% filter(thinning_treatment == 'variabl_thin') %>% st_union()
)
poly_thinning_treatment = bind_rows(commrcl_sf, variabl_sf)

data_plot_scale = st_join(data_plot_scale, poly_thinning_treatment)
table(data_plot_scale$thinning_treatment, useNA = 'ifany')

# Create a "stratum" rast from rast_stage_class,
rast_stratum = rast_stage_class

# Overlap rast_stratum with thinning treatment as an additional class
rast_thinning = rasterize(vect(poly_thinning_treatment), rast(ext(rast_stratum), resolution = res(rast_stratum), crs = crs(rast_stratum)), field = "thinning_status", background = NA)
values(rast_thinning)[!is.na(values(rast_thinning))] = 5
rast_stratum = cover(rast_thinning, rast_stratum)
data_plot_scale$stratum = terra::extract(rast_stratum, vect(data_plot_scale))[, 2]
table(data_plot_scale$stratum, useNA = 'ifany')

# TODO: Resolve the small number of sites that are ambiguously classified due to origin age near class boundaries:
# Dz303i > age 24, stand init / canopy closure / WADNR stemex
# Ca263i > age 86, understory reinit, WADNR stemex
# Dp166i > age 24, stand init / canopy closure / WANDR stemex

# From here forward, we use "rast_stratum" as the classification scheme to delineate patches
mapview(vect_patches) + mapview(poly_thinning_treatment)
summary(data_plot_scale$stratum)
mapview(rast_stratum) +
  mapview(aru_sites, label = aru_sites$site)

# Further delineate patch boundaries by roads, waterbodies, and watercourses

# Paved roads. Unpaved roads do not constitute a boundary of a patch because they are narrow, rarely driven and likely experienced by the birds as gaps in the forest.
roads_wadnr = st_read('data/environment/GIS Data/roads/T3roads.shp') %>%
  st_zm(drop = TRUE) %>% janitor::clean_names() %>% st_transform(crs = crs_m) %>% select(road_usgs1, road_surfa, geometry) %>%
  mutate(geometry = st_cast(geometry, "MULTILINESTRING"))
roads_wsdot = st_read('data/environment/WSDOT_-_Local_Agency_Public_Road_Routes/WSDOT_-_Local_Agency_Public_Road_Routes.shp') %>% filter(RouteIdent %in% c("400000220i", "031265969i")) %>% st_transform(crs_m) %>% select(geometry) # get Hoh Mainline Road / Clearwater Road from WSDOT
roads_wsdot$road_usgs1 = 'Primary Highway'
roads_wsdot$road_surfa = NA
roads = rbind(roads_wsdot, roads_wadnr)
# mapview(roads, zcol = 'road_usgs1') + mapview(aru_sites, label = aru_sites$site)

# "The Hoh-Clearwater Mainline is our only double-lane mainline, and it is about 26 feet wide for the asphalt surface. If we’re looking at right of way widths, we could easily assume a minimum width of about 50 feet for a 12-foot road, probably a good 60-80 feet for a 14-20 foot wide road, and about 100 feet for the Hoh Mainline. 100 feet might be good for US 101 as well for right of way width, and maybe about 30 feet for actual road surface width."
paved_primary_roads   = st_make_valid(roads %>% filter(road_usgs1 %in% c("Primary Highway")))
road_half_width_primary = max(100 / 2 * conv_ft_to_m, res(rast_stratum)[1] / 2)
# mapview(paved_primary_roads) + mapview(st_union(st_buffer(paved_primary_roads, road_half_width_primary)))
# "Our minimum road surface width is going to be 12 feet. That would cover most of our roads. Main arterials/single lane mainlines will have a minimum surface width of 14 feet, with a few up to 20 feet."
paved_secondary_roads = st_make_valid(roads %>% filter(road_usgs1 %in% c("Light-Duty Road")))
road_half_width_secondary = max(30 / 2 * conv_ft_to_m, res(rast_stratum)[1] / 2)
# mapview(paved_secondary_roads) + mapview(st_union(st_buffer(paved_secondary_roads, road_half_width_secondary)))

template = rast(ext(rast_stratum), resolution = res(rast_stratum), crs = crs(rast_stratum))
paved_roads_primary_buffered = (st_buffer(paved_primary_roads, dist = road_half_width_primary))
paved_roads_secondary_buffered = (st_buffer(paved_secondary_roads, dist = road_half_width_secondary))
rast_paved_roads_primary = rasterize(vect(paved_roads_primary_buffered), template, field = 6, background = NA, touches=TRUE)
rast_paved_roads_secondary = rasterize(vect(paved_roads_secondary_buffered), template, field = 6, background = NA, touches=TRUE)
# mapview(rast_paved_roads_primary) + mapview(paved_roads_primary_buffered)
# mapview(rast_paved_roads_secondary) + mapview(paved_roads_secondary_buffered)
# mapview(rast_paved_roads_primary) + mapview(rast_paved_roads_secondary)
rast_updated = rast_stratum
rast_updated = cover(rast_paved_roads_primary, rast_updated)
rast_updated = cover(rast_paved_roads_secondary, rast_updated)

# Watercourses/waterbodies (rivers and streams) from Type 1-3. These support distinct riparian vegetation which constitutes different habitat. Type 4 (non-fish-bearing streams) and type 5 are not boundaries because they are small features that do not change the vegetation, support less aquatic biota, and often are not permanent.
watercourses = st_read('data/environment/DNR_Hydrography/DNR_Hydrography_-_Watercourses_-_Forest_Practices_Regulation/DNR_Hydrography_-_Watercourses_-_Forest_Practices_Regulation.shp')
watercourses = watercourses %>% st_crop(st_transform(study_area, st_crs(watercourses))) %>% st_transform(crs_m) %>% janitor::clean_names() %>% select(geometry, sl_wtrty_c)
# mapview(watercourses, zcol = 'sl_wtrty_c')

boundary_watercourses = watercourses %>% filter(sl_wtrty_c %in% c(1, 2, 3))
boundary_watercourses$sl_wtrty_c = 1
watercourse_half_width = res(rast_stratum)[1] / 2
watercourses_buffered = (st_buffer(boundary_watercourses, dist = watercourse_half_width))
rast_watercourses = rasterize(vect(watercourses_buffered), template, field = 7, background = NA, touches=TRUE)
# mapview(rast_watercourses) + mapview(watercourses_buffered)
rast_updated = cover(rast_watercourses, rast_updated)

waterbodies = st_read('data/environment/DNR_Hydrography/DNR_Hydrography_-_Water_Bodies_-_Forest_Practices_Regulation/DNR_Hydrography_-_Water_Bodies_-_Forest_Practices_Regulation.shp')
waterbodies = waterbodies %>% st_crop(st_transform(study_area, st_crs(waterbodies))) %>% st_transform(crs_m) %>% janitor::clean_names()
# mapview(waterbodies, zcol = 'sl_wtrty_c')

boundary_waterbodies = waterbodies %>% filter(sl_wtrty_c %in% c(1))
rast_waterbodies = rasterize(vect(boundary_waterbodies), template, field = 7, background = NA, touches=TRUE)
# mapview(rast_waterbodies) + mapview(boundary_waterbodies)

# Fill missing water bodies off coast
ocean_polyfill = st_polygon(list(rbind(
  c(-180, -90), c(-124.42, -90), c(-124.42, 47.63), c(-180, 47.63), c(-180, -90)
))) |> st_sfc(crs = 4326) |> st_transform(crs(rast_waterbodies))
p = terra::project(vect(ocean_polyfill), crs(rast_waterbodies))
mask_r = rasterize(p, rast_waterbodies, field=1, background=NA, touches=TRUE)
rast_waterbodies[!is.na(mask_r[])] <- 7

rast_updated = cover(rast_waterbodies, rast_updated)
mapview(rast_updated)

# TODO: Cover impervious surfaces

# Impute any remaining missing patches with cover class 2 (stem exclusion)
rast_updated[is.na(rast_updated[])] <- 2
mapview(rast_updated)

# Crop raster to the study area (bounded by maximum home range size)
max_homerange_buffers = st_buffer(aru_sites, max(species_trait_data$home_range_radius_m))
study_region_bbox = st_bbox(st_union(max_homerange_buffers))

r_cropped = crop(rast_updated, ext(study_region_bbox$xmin, study_region_bbox$xmax, study_region_bbox$ymin, study_region_bbox$ymax))

mapview(r_cropped) +
  mapview(rast_updated) +
  mapview(aru_sites, label = aru_sites$site) +
  mapview(st_buffer(aru_sites, 100), col.regions = 'transparent', lwd = 2) +
  mapview(max_homerange_buffers, col.regions = 'transparent', lwd = 2) +
  mapview(study_region_bbox)

rast_updated = r_cropped

# Clean raster by generalizing minimum patch area and width
if (overwrite_rast_cover_cache) {
  
  message("Generating cover cache (current time ", time_start <- Sys.time(), ")")
  
  min_species_homerange_area_m2 = round(pi * min(species_trait_data$home_range_radius_m)^2,0)

  r = rast_updated
  cell_res = res(r)[1] # cell resolution (in meters)
  min_area = min_species_homerange_area_m2 # OR 0.785 * 1e4 --> e.g. 0.785 hectares to square meters
  min_width = 50 # in meters
  min_cells_area = ceiling(min_area / (cell_res^2))  # min number of cells by area
  min_cells_width = ceiling(min_width / cell_res)    # min number of cells by width
  
  # Identify discrete contiguous patches with unique ids
  message("Identifying discrete contiguous patch candidates")
  patch_list = list()
  start_id = 1  # Starting patch ID to ensure uniqueness
  cls = sort(na.omit(unique(values(r))))
  for (cl in cls) {
    message("Class: ", cl)
    r_class = mask(r, r == cl, maskvalues = FALSE) # Mask only the current class
    r_patches = patches(r_class, directions = 8) # Compute patches for this class only
    # Reclassify patch IDs to be globally unique
    max_patch_id = global(r_patches, "max", na.rm = TRUE)[1,1]
    if (!is.na(max_patch_id)) {
      r_patches = classify(r_patches, cbind(1:max_patch_id, start_id:(start_id + max_patch_id - 1)))
      start_id = start_id + max_patch_id
    }
    patch_list[[as.character(cl)]] <- r_patches
  }
  
  # Merge patch ids for each raster type into one raster
  patch_ids = cover(patch_list[[1]], patch_list[[2]])
  patch_ids = cover(patch_ids, patch_list[[3]])
  patch_ids = cover(patch_ids, patch_list[[4]])
  patch_ids = cover(patch_ids, patch_list[[5]])
  patch_ids = cover(patch_ids, patch_list[[6]])
  patch_ids = cover(patch_ids, patch_list[[7]])
  
  # Create a raster stack (cover class and patch id)
  patch_stack = c(patch_ids, r)
  names(patch_stack) = c("patch_id", "cover_class")
  
  mapview(patch_stack[["patch_id"]], col.regions = viridis) + 
    mapview(patch_stack[["cover_class"]])
  
  # Identify patches with maximum width less than minimum width ##################################################
  message("Identifying patch candidates with negligible width")
  narrow_patches = c()
  pids = unique(na.omit(as.vector(values(patch_stack[['patch_id']]))))
  for (pid in pids) {
    print(pid)
    m = trim(classify(patch_stack[['patch_id']], cbind(pid, 1), others = NA)) # Mask patch
    m_vals = as.matrix(m, wide=TRUE)
    
    if (all(is.na(m_vals))) next  # Skip empty masks
    
    max_h_run = max(apply(m_vals, 1, function(row) { # Horizontal runs (per row)
      rle_row = rle(row)
      max(c(0, rle_row$lengths[rle_row$values == 1]), na.rm = TRUE)
    }), na.rm = TRUE)
    max_v_run = max(apply(m_vals, 2, function(col) { # Vertical runs (per column)
      rle_col = rle(col)
      max(c(0, rle_col$lengths[rle_col$values == 1]), na.rm = TRUE)
    }), na.rm = TRUE)
    if (max_h_run < min_cells_width || max_v_run < min_cells_width) { # Keep if both directions have max run < 3
      narrow_patches = c(narrow_patches, pid)
    }
  }
  
  # Identify patch ids with negligible patch width
  patch_stack[["patch_narrow"]] = classify(patch_stack[["patch_id"]], rcl = cbind(narrow_patches, narrow_patches), others = NA)
  patch_narrow_ids = unique(patch_stack[["patch_narrow"]])[,1]
  
  # Reclassify those cells to value 0
  values(patch_stack[["cover_class"]])[which(values(!is.na(patch_stack[['patch_narrow']])))] = 0
  mapview(patch_stack[['cover_class']])
  
  # Identify patches with area less than minimum area ##################################################
  message("Identifying patch candidates with negligible area")
  
  # Identify patches with negligible area
  freq_table = freq(patch_stack[["patch_id"]])
  small_patches = freq_table[freq_table$count < min_cells_area, "value"]
  patch_stack[["patch_small"]] = classify(patch_stack[["patch_id"]], rcl = cbind(small_patches, small_patches), others = NA)
  patch_small_ids = unique(patch_stack[["patch_small"]])[,1]
  
  # Reclassify those cells to value 0
  values(patch_stack[["cover_class"]])[which(values(!is.na(patch_stack[['patch_small']])))] = 0
  mapview(patch_stack[['cover_class']])
  
  r_zero_alt = patch_stack[['cover_class']]
  values(r_zero_alt)[values(r_zero_alt) != 0] <- NA
  mapview(r_zero_alt)
  p = patches(r_zero_alt, directions = 8, values=TRUE)
  
  # Reclassify (generalize) negligible patches as the mode of adjacent nearest neighbors
  message("Reclassifying negligible patches as the mode of adjacent nearest neighbors")
  rast_cover_clean = patch_stack[["cover_class"]]
  i = 1
  patch_ids_to_generalize = unique(na.omit(values(p)))
  total = length(patch_ids_to_generalize)
  for (small_patch_id in patch_ids_to_generalize) {
    print(round(i / total, 3))
    
    # mapview(trim(classify(patch_stack[["patch_id"]], cbind(small_patch_id, 1), others = NA)))
    
    patch_cells = which(values(p) == small_patch_id)
    
    # Find adjacent cells to the patch
    adj_cells = adjacent(rast_stratum, cells = patch_cells, directions = 4, pairs = TRUE)
    
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
  rast_cover_clean[is.nan(rast_cover_clean)] = NA
  mapview(rast_cover_clean) + mapview(patch_stack[['cover_class']])

  message('Saving raster cover data cache ', path_rast_cover_clean_out)
  dir.create(dirname(path_rast_cover_clean_out), recursive = TRUE, showWarnings = FALSE)
  writeRaster(rast_cover_clean, path_rast_cover_clean_out, overwrite=TRUE)
  message("Finished generating cover cache (", round(as.numeric(difftime(Sys.time(), time_start, units = 'mins')), 2), " minutes)")
  
} else { # overwrite_rast_cover_cache is FALSE
  message('Loading raster cover data from cache ', path_rast_cover_clean_out)
  rast_cover_clean = rast(path_rast_cover_clean_out)
}
rast_cover_clean = as.factor(rast_cover_clean)

mapview(rast_cover_clean,
        alpha.regions = 1.0,
        col.regions = c('#90c6bd', '#3c8273', '#d8c18a', '#9b652b', '#b2675e', 'darkgray', '#6495ed')) +
  mapview(aru_sites, label = aru_sites$site) +
  mapview(st_buffer(aru_sites, 100), col.regions = 'transparent', lwd = 2) +
  mapview(max_homerange_buffers, col.regions = 'transparent', lwd = 2)

##############################################################################
# Local plot scale covariates

if (overwrite_data_plot_scale_cache) {
  
  message("Generating plot scale data cache (current time ", time_start <- Sys.time(), ")")

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
  
  # Mark sites that were surveyed for habitat data in-person
  intersect(data_plot_scale$site, data_plot$site)
  length(intersect(data_plot_scale$site, data_plot$site))
  setdiff(data_plot$site, data_plot_scale$site) # This should be 0
  data_plot_scale[data_plot_scale$site %in% data_plot$site, 'hs'] = TRUE
  mapview(data_plot_scale, zcol = 'hs')
  
  # Elevation [m] (derived from WA state LiDAR flights 3x3-ft raster grid)
  elevation = read.csv('data/environment/site_elevation_ft.csv')
  elevation$elev = elevation$elev_ft * conv_ft_to_m
  data_plot_scale = data_plot_scale %>% left_join(elevation %>% select(site, elev), by = 'site')
  
  # Age (mean and cv) [#]
  data_plot_scale$age_point = terra::extract(rast_age, vect(data_plot_scale))[, 2]
  stats_age = compute_raster_buffer_value_func(rast_age, data_plot_scale, plot_buffer, summary_stats)
  data_plot_scale$age_mean = as.numeric(stats_age[,2])
  data_plot_scale$age_cv   = as.numeric(stats_age[,4])
  hist(data_plot_scale$age_point, breaks = seq(0, max(data_plot_scale$age_point) + 10, by = 10))
  hist(data_plot_scale$age_mean, breaks = seq(0, max(data_plot_scale$age_mean) + 10, by = 10))
  table(data_plot_scale$age_point)
  table(round(data_plot_scale$age_mean))
  
  ggplot(data_plot_scale, aes(x = stage, y = age_mean)) + geom_boxplot() + theme_minimal()
  
  # Cover class [categorical]
  data_plot_scale$stage
  
  # Thinning treatment [categorical]
  data_plot_scale$thinning_treatment
  
  # Basal area (all live trees) [m2/ha] TODO: confirm if this is a mean value at local plot level
  data_plot_scale = data_plot_scale %>% left_join(data_plot %>% select(site, plot_ba_hs = ba_ha_all), by = 'site')
  
  # Tree density (all live trees large, dbh > 10 cm) [# trees/ha]
  data_plot_scale = data_plot_scale %>% left_join(data_plot %>% select(site, plot_treeden_gt10cmDbh_hs = large_per_hectare_all), by = 'site')
  
  # Tree density (small, dbh < 10 cm) [# trees/ha]
  data_plot_scale = data_plot_scale %>% left_join(data_plot %>% select(site, plot_treeden_lt10cmDbh_hs = small_per_hectare), by = 'site')
  
  # Total tree density (all sizes) [# trees/ha]
  data_plot_scale$plot_treeden_all_hs = data_plot_scale$plot_treeden_gt10cmDbh_hs + data_plot_scale$plot_treeden_lt10cmDbh_hs
  
  # Tree quadratic mean diameter (large, dbh > 10cm) [cm]
  data_plot_scale = data_plot_scale %>% left_join(
    data_plot %>% select(site, plot_qmd_gt10cmDbh_hs = avg_dbh_cm_all) %>%
      mutate(plot_qmd_gt10cmDbh_hs = replace_na(plot_qmd_gt10cmDbh_hs, 0.0)),
    by = 'site')
  # Tree quadratic mean diameter (small, dbh < 10 cm) [cm]
  data_plot_scale = data_plot_scale %>% left_join(
    data_plot %>% select(site, plot_qmd_lt10cmDbh_hs = avg_dbh_cm) %>%
      mutate(plot_qmd_lt10cmDbh_hs = replace_na(plot_qmd_lt10cmDbh_hs, 0.0)),
    by = 'site')
  # Tree quadratic mean diameter (all) [cm]
  data_plot_scale$plot_qmd_all_hs = data_plot_scale$plot_qmd_gt10cmDbh_hs + data_plot_scale$plot_qmd_lt10cmDbh_hs
  
  # Tree height (mean and cv) [m]
  # TODO: Get high-resolution LiDAR height data from DNR approval procedure
  data_plot_scale = data_plot_scale %>% left_join(
    data_plot %>% select(site, plot_ht_hs = avg_height_m_all) %>%
      mutate(plot_ht_hs = replace_na(plot_ht_hs, 0.0)),
    by = 'site')
  data_plot_scale = data_plot_scale %>% left_join(
    data_plot %>% select(site, plot_ht_cv_hs = cv_height_all) %>%
      mutate(plot_ht_cv_hs = replace_na(plot_ht_cv_hs, 0.0)),
    by = 'site')
  
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
  
  # Downed wood volume [m3/ha]
  data_plot_scale = data_plot_scale %>% left_join(data_plot %>% select(site, plot_downvol_hs = vol_alldown_m3) %>% mutate(plot_downvol_hs = replace_na(plot_downvol_hs, 0.0)), by = 'site')
  
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
  plot_tree_gte10cm_richness  = specnumber(tree_gte10cm_obs %>% select(-site))
  plot_tree_gte10cm_diversity = diversity(tree_gte10cm_obs %>% select(-site), index = 'shannon')
  plot_tree_gte10cm_evenness  = plot_tree_gte10cm_diversity / log(plot_tree_gte10cm_richness)
  
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
  plot_tree_lt10cm_richness  = specnumber(tree_lt10cm_obs %>% select(-site))
  plot_tree_lt10cm_diversity = diversity(tree_lt10cm_obs %>% select(-site), index = 'shannon')
  plot_tree_lt10cm_evenness  = plot_tree_lt10cm_diversity / log(plot_tree_lt10cm_richness)
  
  tree_all_obs = data.frame(
    site = data_plot %>% pull(site),
    all_per_hectare_psme = tree_gte10cm_obs$large_per_hectare_psme + tree_lt10cm_obs$small_per_hectare_psme,
    all_per_hectare_thpl = tree_gte10cm_obs$large_per_hectare_thpl + tree_lt10cm_obs$small_per_hectare_thpl,
    all_per_hectare_abam = tree_gte10cm_obs$large_per_hectare_abam + tree_lt10cm_obs$small_per_hectare_abam,
    all_per_hectare_tshe = tree_gte10cm_obs$large_per_hectare_tshe + tree_lt10cm_obs$small_per_hectare_tshe,
    all_per_hectare_alru = tree_gte10cm_obs$large_per_hectare_alru + tree_lt10cm_obs$small_per_hectare_alru,
    all_per_hectare_pisi = tree_gte10cm_obs$large_per_hectare_pisi + tree_lt10cm_obs$small_per_hectare_pisi
  )
  plot_tree_all_richness  = specnumber(tree_all_obs %>% select(-site))
  plot_tree_all_diversity = diversity(tree_all_obs %>% select(-site), index = 'shannon')
  plot_tree_all_evenness  = plot_tree_all_diversity / log(plot_tree_all_richness)
  
  tree_div_metrics = data.frame(
    site = data_plot %>% pull(site),
    plot_tree_all_richness,  plot_tree_gte10cm_richness,  plot_tree_lt10cm_richness,
    plot_tree_all_diversity, plot_tree_gte10cm_diversity, plot_tree_lt10cm_diversity,
    plot_tree_all_evenness,  plot_tree_gte10cm_evenness,  plot_tree_lt10cm_evenness,
    plot_tree_gte10cm_density_psme = tree_gte10cm_obs$large_per_hectare_psme,
    plot_tree_gte10cm_density_thpl = tree_gte10cm_obs$large_per_hectare_thpl,
    plot_tree_gte10cm_density_abam = tree_gte10cm_obs$large_per_hectare_abam,
    plot_tree_gte10cm_density_tshe = tree_gte10cm_obs$large_per_hectare_tshe,
    plot_tree_gte10cm_density_alru = tree_gte10cm_obs$large_per_hectare_alru,
    plot_tree_gte10cm_density_pisi = tree_gte10cm_obs$large_per_hectare_pisi,
    plot_tree_all_density_psme = tree_all_obs$all_per_hectare_psme,
    plot_tree_all_density_thpl = tree_all_obs$all_per_hectare_thpl,
    plot_tree_all_density_abam = tree_all_obs$all_per_hectare_abam,
    plot_tree_all_density_tshe = tree_all_obs$all_per_hectare_tshe,
    plot_tree_all_density_alru = tree_all_obs$all_per_hectare_alru,
    plot_tree_all_density_pisi = tree_all_obs$all_per_hectare_pisi
  )
  
  data_plot_scale = data_plot_scale %>% left_join(tree_div_metrics, by = 'site')
  
  ggplot(data_plot_scale, aes(x = as.factor(`stratum`), y = plot_tree_all_diversity)) +
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
  metric = "plot_tree_all_density_" # "tree_all_density_", "tree_gte10cm_density_"
  dps_long = data_plot_scale %>% select(site, `stage`, starts_with(metric)) %>%
    pivot_longer(cols = starts_with(metric), names_to = "species", values_to = "density") %>%
    mutate(species = gsub(metric, "", species))
  ggplot(dps_long, aes(x = species, y = density, fill = species)) +
    geom_boxplot() +
    facet_wrap(~ stage, scales = "free_y") +
    coord_cartesian(ylim = c(0, 1250)) + 
    labs(title = "Tree species density by stage") +
    theme_minimal()
  
  # Distance to roads
  dist_road_paved = st_distance(data_plot_scale, paved_primary_roads)
  dist_road_paved = apply(dist_road_paved, 1, min)
  data_plot_scale$dist_road_paved = dist_road_paved
  mapview(data_plot_scale, zcol = "dist_road_paved") + mapview(paved_primary_roads)
  
  dist_roads_all = st_distance(data_plot_scale, roads)
  dist_roads_all = apply(dist_roads_all, 1, min)
  data_plot_scale$dist_roads_all = dist_roads_all
  mapview(data_plot_scale, zcol = "dist_roads_all") + mapview(paved_primary_roads)
  
  # Distance to water (type 1-3 and all types) [m]
  dist_watercourse_major = st_distance(data_plot_scale, boundary_watercourses)
  dist_watercourse_major = apply(dist_watercourse_major, 1, min)
  data_plot_scale$dist_watercourse_major = dist_watercourse_major
  mapview(data_plot_scale, zcol = "dist_watercourse_major") + mapview(boundary_watercourses)
  
  dist_watercourse_all = st_distance(data_plot_scale, watercourses)
  dist_watercourse_all = apply(dist_watercourse_all, 1, min)
  data_plot_scale$dist_watercourse_all = dist_watercourse_all
  mapview(data_plot_scale, zcol = "dist_watercourse_all") + mapview(watercourses)
  
  # # Distance to nearest edge [m]
  # # For each site, find its patch, then find the distance to the nearest non-patch cell
  # patch_ids = patches(rast_cover_clean, directions=8, values=TRUE)
  # data_plot_scale$dist_nearest_edge = NA
  # for (i in 1:nrow(data_plot_scale)) {
  #   d = data_plot_scale[i,]
  #   cover_class = terra::extract(rast_cover_clean, vect(d))[,2]
  #   pid = terra::extract(patch_ids, vect(d))[,2]
  #   m = trim(classify(patch_ids, cbind(pid, 1), others = NA))
  #   na_cells = which(is.na(values(m)))
  #   na_coords = xyFromCell(m, na_cells)
  #   d_coords = st_coordinates(d)
  #   distances = sqrt((na_coords[,1] - d_coords[1])^2 + (na_coords[,2] - d_coords[2])^2)
  #   min_dist = min(distances)
  #   data_plot_scale[i, 'dist_nearest_edge'] = min_dist
  #   # mapview(d) + mapview(m) + mapview(st_buffer(d, 100), col.regions = 'transparent', lwd = 2)
  # }
  # mapview(rast_cover_clean) + mapview(data_plot_scale, label = data_plot_scale$site, zcol = 'dist_nearest_edge') + mapview(st_buffer(data_plot_scale, 100), col.regions = 'transparent', lwd = 2)
  
  message('Saving plot scale data cache ', path_data_plot_scale_out)
  dir.create(dirname(path_data_plot_scale_out), recursive = TRUE, showWarnings = FALSE)
  saveRDS(data_plot_scale, path_data_plot_scale_out)
  message("Finished generating plot scale data cache (", round(as.numeric(difftime(Sys.time(), time_start, units = 'mins')), 2), " minutes)")

} else { # overwrite_data_plot_scale_cache is FALSE
  message('Loading plot scale data from cache ', path_data_plot_scale_out)
  data_plot_scale = readRDS(path_data_plot_scale_out)
}

##############################################################################
# Focal patch scale and home range neighborhood scale covariates (limited to the extent of the species home range)

if (overwrite_data_homerange_scale_cache) {
  
  message("Generating homerange scale data cache (current time ", time_start <- Sys.time(), ")")

  data_homerange_scale = list()
  
  homeranges = species_trait_data %>% select(common_name, home_range_radius_m)
  message("Homerange buffer min: ",    min(homeranges$home_range_radius_m))
  message("Homerange buffer median: ", buff_median <- round(median(homeranges$home_range_radius_m),0))
  message("Homerange buffer mean: ",   buff_mean <- round(mean(homeranges$home_range_radius_m),0))
  message("Homerange buffer max: ",    buff_max <- round(max(homeranges$home_range_radius_m),0))
  homeranges = rbind(data.frame(
    common_name = c("plot", "median", "mean", "max"),
    home_range_radius_m = c(100, buff_median, buff_mean, buff_max)
  ), homeranges %>% mutate(home_range_radius_m = round(home_range_radius_m, 0)))
  
  # DEBUG
  # Overwrite homeranges with the median and mean range only for now
  # homeranges = data.frame(
  #   common_name = c("min", "median", "mean", "max"),
  #   home_range_radius_m = c(100, buff_median, buff_mean, buff_max)
  # )
  # DEBUG
  
  ## Load forest structure rasters
  message("Loading forest structure rasters")
  
  # Basal area [ft2/acre] --> [m2/ha]
  rast_ba = load_raster(paste0(dir_rsfris_version, '/RS_FRIS_BA_EXPORT.tif')) * 0.229568
  # Tree density (total) [#/acre] --> [#/ha]
  rast_treeden_all = load_raster(paste0(dir_rsfris_version, '/RS_FRIS_TREE_ACRE.img'))  * conv_peracre_to_perha
  # Tree density (DBH > 4 in) [#/acre]-->[#/ha]
  rast_treeden_gt4inDbh = load_raster(paste0(dir_rsfris_version, '/RS_FRIS_TREE_ACRE_4.img')) * conv_peracre_to_perha
  # Tree quadratic mean diameter [in]-->[cm]
  rast_qmd = load_raster(paste0(dir_rsfris_version, '/RS_FRIS_QMD.img')) * conv_in_to_cm
  # Tree height (max) [ft]-->[m]
  rast_htmax = load_raster(paste0(dir_rsfris_version, '/RS_FRIS_HTMAX.img')) * conv_ft_to_m
  # Canopy layers [#]
  rast_canopy_layers = load_raster(paste0(dir_rsfris_version, '/RS_FRIS_CANOPY_LAYERS.img'))
  # Canopy cover [%]
  rast_canopy_cover = load_raster(paste0(dir_rsfris_version, '/RS_FRIS_COVER.img'))
  # Canopy closure [%]
  rast_canopy_closure = load_raster(paste0(dir_rsfris_version, '/RS_FRIS_CLOSURE.img'))
  # Density of snags > 15" DBH [#/acre] --> [#/ha]
  rast_snagden_gt15inDbh = load_raster(paste0(dir_rsfris_version, '/RS_FRIS_SNAG_ACRE_15.img')) * conv_peracre_to_perha
  # Downed wood volume [ft3/acre] --> [m3/ha]
  rast_downvol = load_raster(paste0(dir_rsfris_version, '/RS_FRIS_CFVOL_DDWM.img')) * 0.069968
  
  # Calculate per scale according to predicted home range size
  message("Calculating homerange variables across all scales")
  for (i in 1:nrow(homeranges)) {
    scale = homeranges[i,'common_name']
    homerange_buffer_size = homeranges[i,'home_range_radius_m']
    message("Scale '", scale, "', home range buffer size: ", homerange_buffer_size)
    
    data_homerange_scale_species = data.frame()
  
    # "Edge influences on distributions of organisms or factors affecting organisms (such as predation and nest parasitism) are concentrated within 50m of the edge." (Kremaster L. & Bunnell F. L. (1999) Edge effects: Theory, evidence and implications to management of western North American forests. In: Forest Fragmentation: Wildlife and Management Implications (eds J. Wisniewski, J. A. Rochelle & L. Lehmann) pp. 117–53. Leiden, Boston, MA.)
    core_area_buffer = 50 # meters
    core_area_buffer_ncells = round(core_area_buffer / res(rast_cover_clean)[1], 0)
    
    sites = 1:nrow(data_plot_scale)
    pb = progress_bar$new(format = "[:bar] :percent :elapsedfull (ETA :eta)", total = nrow(data_plot_scale), clear = FALSE)
    for (j in sites) {
      
      ## Get site cover data within and surrounding homerange buffer
      site = data_plot_scale[j,]
      # message(site$site, ' [j=', j,']', ' (', round(j / length(sites),3), ')')
      
      homerange_and_edge_buffer = st_buffer(site, homerange_buffer_size + core_area_buffer) # additional buffer to ensure that edges are retained in mask
      homerange_buffer = st_buffer(site, homerange_buffer_size)
      homerange_and_edge_crop = crop(rast_cover_clean, vect(homerange_and_edge_buffer))
      homerange_and_edge = mask(homerange_and_edge_crop, vect(homerange_and_edge_buffer))
      homerange_crop = crop(rast_cover_clean, vect(homerange_buffer))
      homerange_cover = mask(homerange_crop, vect(homerange_buffer))
      rast_homerange = ifel(!is.na(homerange_cover), 1, NA)
      # mapview(homerange_cover) + mapview(homerange_and_edge) + mapview(homerange_buffer)

      r = homerange_and_edge
      patch_ids = r 
      values(patch_ids) <- NA 
      start_id = 1  # Starting patch ID to ensure uniqueness
      cls = sort(na.omit(unique(values(r))))
      for (cl in cls) {
        # message("Class: ", cl)
        r_class = mask(r, r == cl, maskvalues = FALSE) # Mask only the current class
        r_patches = patches(r_class, directions = 8) # Compute patches for this class only
        # Reclassify patch IDs to be globally unique
        max_patch_id = global(r_patches, "max", na.rm = TRUE)[1,1]
        if (!is.na(max_patch_id)) {
          r_patches = classify(r_patches, cbind(1:max_patch_id, start_id:(start_id + max_patch_id - 1)))
          start_id = start_id + max_patch_id
        }
        patch_ids = cover(patch_ids, r_patches)
      }
      pid = terra::extract(patch_ids, vect(site))[,2]
      rast_focal_patch = trim(classify(patch_ids, cbind(pid, 1), others = NA))
      rast_focal_patch = crop(rast_focal_patch, vect(homerange_buffer))
      rast_focal_patch = mask(rast_focal_patch, vect(homerange_buffer))
      # mapview(rast_focal_patch) + mapview(homerange_cover)
      
      ## Focal patch scale
      
      # Focal patch cover class
      (focalpatch_cover = terra::extract(rast_cover_clean, vect(site))[,2])
      
      # Focal patch area [ha]
      ncells_focal_patch = sum(!is.na(values(rast_focal_patch)))
      (focal_patch_area = ncells_focal_patch * prod(res(rast_focal_patch)) * conv_m2_to_ha)
      
      # Focal patch core area [ha]
      focal_patch_and_edge_buffer = trim(classify(patch_ids, cbind(pid, 1), others = NA))
      focal_patch_inv = ifel(is.na(focal_patch_and_edge_buffer), 1, NA)
      focal_patch_edge_dist = distance(focal_patch_inv)
      focal_patch_edge_dist = ifel(focal_patch_edge_dist == 0, NA, focal_patch_edge_dist)
      focal_patch_core = ifel(focal_patch_edge_dist >= core_area_buffer, 1, NA)
      focal_patch_core = crop(focal_patch_core, vect(homerange_buffer))
      focal_patch_core = mask(focal_patch_core, vect(homerange_buffer))
      # mapview(homerange_buffer, col.regions = 'transparent', lwd = 2) + mapview(focal_patch_edge_dist) + mapview(rast_focal_patch) + mapview(focal_patch_and_edge_buffer) + mapview(focal_patch_core)
      ncells_focal_patch_core = sum(!is.na(values(focal_patch_core)))
      (focal_patch_core_area = ncells_focal_patch_core * prod(res(rast_focal_patch)) * conv_m2_to_ha)
      
      # Structural variables (across focal patch, then across homerange)
      results_list = list()
      regions = list(
        focalpatch = rast_focal_patch,
        homerange  = rast_homerange
      )
      for (k in seq_along(regions)) {
        region_name = names(regions)[k]
        # message("Calculating forest structural variables for: ", region_name)
        region = resample(regions[[k]], rast_age, method = "near") # align region to raster data
        # print(sum(!is.na(values(region))))
        
        vars = c(
          # Age [#]
          age = summary_stats(values(mask(rast_age, region)), na.rm = TRUE),
          # Basal area [m2/ha]
          ba = summary_stats(values(mask(rast_ba, region)), na.rm = TRUE),
          # Tree density (total) [# trees/ha]
          treeden_all = summary_stats(values(mask(rast_treeden_all, region)), na.rm = TRUE),
          # Tree density (large, dbh > 4 in) [# trees/ha]
          treeden_gt4inDbh = summary_stats(values(mask(rast_treeden_gt4inDbh, region)), na.rm = TRUE),
          # Tree diameter [cm]
          qmd = summary_stats(values(mask(rast_qmd, region)), na.rm = TRUE),
          # Tree height [m]
          htmax = summary_stats(values(mask(rast_htmax, region)), na.rm = TRUE),
          # Canopy layers [#]
          canopy_layers = summary_stats(values(mask(rast_canopy_layers, region)), na.rm = TRUE),
          # Canopy cover [%]
          canopy_cover = summary_stats(values(mask(rast_canopy_cover, region)), na.rm = TRUE),
          # Canopy closure [%]
          canopy_closure = summary_stats(values(mask(rast_canopy_closure, region)), na.rm = TRUE),
          # Density of snags > 15" DBH [#/ha]
          snagden_gt15dbh = summary_stats(values(mask(rast_snagden_gt15inDbh, region)), na.rm = TRUE),
          # Downed wood volume [m3/ha]
          downvol = summary_stats(values(mask(rast_downvol, region)), na.rm = TRUE)
        )
        # format and discard irrelevant variables
        vars = as.data.frame(vars) %>% select(-ends_with(".sd"))
        results_list[[region_name]] = vars
      }
      final_df = do.call(cbind, results_list) %>% janitor::clean_names()
    
      ## Homerange scale metrics
      
      # message("Calculating homerange configuration and composition variables")
      
      homerange_cover_forest = homerange_cover
      homerange_cover_forest[!(homerange_cover_forest[] %in% c(1, 2, 3, 4, 5))] = NA
      # mapview(homerange_cover) + mapview(homerange_cover_forest)
      ncells_homerange = sum(!is.na(values(homerange_cover)))
      
      # Focal patch percentage of home range [%]
      (focalpatch_area_homeange_pcnt = sum(!is.na(values(rast_focal_patch))) / ncells_homerange)
      
      # Focal patch core percentage of home range [%]
      (focalpatch_core_area_homeange_pcnt = sum(!is.na(values(focal_patch_core))) / ncells_homerange)
  
      # Focal patch euclidean nearest neighbor distance [m]
      # Quantifies habitat isolation
      # Approaches 0 as the distance to the nearest neighbor decreases. Minimum is constrained by the cell size.
      (focal_cover_class = site$stratum)
      matching_cover = rast_cover_clean
      matching_cover[matching_cover != focal_cover_class] <- NA
      focal_patch_extended = extend(rast_focal_patch, matching_cover)
      matching_neighbors = mask(matching_cover, focal_patch_extended, maskvalue=NA, inverse=TRUE)
      focal_patch_distance = distance(focal_patch_extended)
      matching_neighbor_distance = mask(focal_patch_distance, matching_neighbors, maskvalue=NA, inverse=FALSE)
      # mapview(matching_neighbors) + mapview(rast_focal_patch) + mapview(focal_patch_distance) + mapview(matching_neighbor_distance)
      (focalpatch_isolation = min(values(matching_neighbor_distance), na.rm = TRUE) - res(matching_neighbor_distance)[1])
      
      # Cover class richness [#]
      # Quantifies the number of cover types present in home range
      cover_freq = freq(homerange_cover) %>% select(value, count) %>% mutate(value = as.integer(value))
      cover_freq = data.frame(value = 1:7) %>%
        left_join(cover_freq, by = "value") %>%
        mutate(count = ifelse(is.na(count), 0, count))
      cover_forest_freq = cover_freq %>% filter(value %in% c(1,2,3,4,5))
      (cover_richness = cover_freq %>% filter(count != 0) %>% nrow())
      (cover_forest_richness = cover_forest_freq %>% filter(count != 0) %>% nrow())
      
      # Cover class evenness (i.e. dominance) [#]
      # Quantifies the degree of evenness versus dominance in cover type distribution 
      # 0 when only one cover type is present, 1 when types are equally distributed
      (cover_evenness = lsm_l_shei(homerange_cover) %>% pull(value))
      (cover_forest_evenness = lsm_l_shei(homerange_cover_forest) %>% pull(value))
      
      # Cover class diversity (shannon) [#]
      # Quantifies both the richness and evenness of cover type distributions (inversely related to contagion)
      # 0 when only one cover type is present and increases as the number of classes increases while the proportions are equally distributed
      (cover_diversity = lsm_l_shdi(homerange_cover) %>% pull(value))
      (cover_forest_diversity = lsm_l_shdi(homerange_cover_forest) %>% pull(value))
      
      # Proportional abundance of each cover class [%]
      (prop_abund_standinit = cover_freq %>% filter(value == 1)     %>% pull(count) / ncells_homerange) # "stand initiation" 0-25 yr
      (prop_abund_stemexcl = cover_freq %>% filter(value == 2)      %>% pull(count) / ncells_homerange) # "stem exclusion" 25-80 yr
      (prop_abund_undstryreinit = cover_freq %>% filter(value == 3) %>% pull(count) / ncells_homerange) # "understory reinitiation" 80-200 yr
      (prop_abund_oldgrowth = cover_freq %>% filter(value == 4)     %>% pull(count) / ncells_homerange) # "old-growth forest" 200+ yr
      (prop_abund_lsog = prop_abund_undstryreinit + prop_abund_oldgrowth) # "late-successional and old growth forest" (80+ yr)
      
      (prop_abund_comthin = cover_freq %>% filter(value == 5) %>% pull(count) / ncells_homerange) # "commercial thinning"
      (prop_abund_roads = cover_freq %>% filter(value == 6)   %>% pull(count) / ncells_homerange) # "roads"
      (prop_abund_water = cover_freq %>% filter(value == 7)   %>% pull(count) / ncells_homerange) # "water"
      
      # Aggregation index [#]
      # Quantifies degree of habitat contiguity versus fragmentation
      # Equals 0 when the patch types are maximally disaggregated (i.e., when there are no like adjacencies); AI increases as the landscape is increasingly aggregated and equals 100 when the landscape consists of a single patch.
      (aggregation_idx = lsm_l_ai(homerange_cover, directions = 8) %>% pull(value))
      
      # Shape index [#]
      # Quantifies patch shape complexity
      # Equals 1 if all patches are squares. Increases, without limit, as the shapes of patches become more complex.
      (shape_idx = lsm_l_shape_mn(homerange_cover, directions = 8) %>% pull(value))
      
      # Contrast-weighted edge density [m/ha]
      # The density of patch edges weighted by their contrast
      # Equals 0 when there is no edge in the landscape (i.e. landscape consists of a single patch). Increases as the amount of edge in the landscape increases and/or as the contrast in edges increase (i.e. contrast weight approaches 1).
      max_ht_m = max(values(rast_htmax), na.rm = TRUE)
      homerange_height = mask(crop(rast_htmax, homerange_and_edge_buffer), homerange_and_edge_buffer)
      # mapview(patch_ids) + mapview(homerange_height)
      
      height_aligned = resample(homerange_height, patch_ids, method = "bilinear")
      zonal_means = zonal(height_aligned, patch_ids, fun = "median", na.rm = TRUE) # TODO: consider min instead
      
      if (nrow(zonal_means) > 1) {
        # Equals the sum of the lengths (m) of each edge segment in the landscape multiplied by the
        # corresponding contrast weight, divided by the total landscape area (m2), converted to hectares.
        colnames(zonal_means) = c('class', 'height')
        # mapview(classify(patch_ids, rcl = zonal_means))
        edge_lengths = get_adjacencies(patch_ids, neighbourhood = 8, what = "unlike", upper = FALSE)[[1]] * res(patch_ids)[1]
        heights = zonal_means$height
        names(heights) = zonal_means$class
        height_diff_matrix = outer(heights, heights, FUN = function(x, y) abs(x - y))
        
        # Weights are derived from proportion of difference in height relative to the maximum potential difference in height (Hou and Walz 2016, Huang et al. 2014)
        # d = 0 --> 0 m difference
        # d = 1 --> max m difference (i.e. maximum tree height of entire landscape)
        weights = height_diff_matrix / max_ht_m
        homerange_rast_area = ncells_homerange * res(homerange_cover)[1] * res(homerange_cover)[2] * conv_m2_to_ha
        (density_edge_cw = sum(edge_lengths * weights, na.rm = TRUE) / homerange_rast_area)
      } else { 
        # Only one patch type
        (density_edge_cw = 0.0)
      }
      density_edge_cw = set_units(density_edge_cw, 'm/ha')
      
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
      
      final_df = cbind(data.frame(
        scale = scale,
        buffer_radius_m = homerange_buffer_size,
        site = site$site,
        focalpatch_cover,
        focalpatch_area_homeange_pcnt,
        focalpatch_core_area_homeange_pcnt,
        focalpatch_isolation,
        cover_richness,
        cover_forest_richness,
        cover_evenness,
        cover_forest_evenness,
        cover_diversity,
        cover_forest_diversity,
        prop_abund_standinit,
        prop_abund_stemexcl,
        prop_abund_undstryreinit,
        prop_abund_oldgrowth,
        prop_abund_lsog,
        prop_abund_comthin,
        prop_abund_roads,
        prop_abund_water,
        aggregation_idx,
        shape_idx,
        density_edge_cw,
        density_roads_paved,
        density_roads,
        density_streams_major,
        density_streams
      ), final_df)
      
      data_homerange_scale_species = rbind(data_homerange_scale_species, final_df)
      pb$tick()
    }
    
    # Sanity check results
    if (FALSE) {
      site_check = data_homerange_scale_species %>% slice_max(focalpatch_isolation, n = 1, with_ties = FALSE) %>% pull(site)
      site_check = data_plot_scale %>% filter(site == site_check)
      site_check_buffer = st_buffer(site_check, homerange_buffer_size)
      mapview(site_check) + mapview(mask(crop(rast_cover_clean, site_check_buffer), site_check_buffer))
      
      site_check = data_homerange_scale_species %>% slice_min(focalpatch_isolation, n = 1, with_ties = FALSE) %>% pull(site)
      site_check = data_plot_scale %>% filter(site == site_check)
      site_check_buffer = st_buffer(site_check, homerange_buffer_size)
      mapview(site_check) + mapview(mask(crop(rast_cover_clean, site_check_buffer), site_check_buffer))
    }
    
    data_homerange_scale[[scale]] = data_homerange_scale_species
  }
  message('Saving homerange data cache ', path_data_homerange_scale_out)
  dir.create(dirname(path_data_homerange_scale_out), recursive = TRUE, showWarnings = FALSE)
  saveRDS(data_homerange_scale, path_data_homerange_scale_out)
  message("Finished generating homerange data cache (", round(as.numeric(difftime(Sys.time(), time_start, units = 'mins')), 2), " minutes)")
  
} else { # overwrite_data_homerange_scale_cache is FALSE
  message('Loading homerange scale data from cache ', path_data_homerange_scale_out)
  data_homerange_scale = readRDS(path_data_homerange_scale_out)
}

#############################################################################################################
# Data inspection at plot scale

data_total = full_join(st_drop_geometry(data_plot_scale), data_homerange_scale[['min']], by = 'site')

comparison_plot = function(x, y, s, title, x_lab = "Habitat survey", y_lab = "Remote sensing") {
  r = cor(x, y, use = "complete.obs", method = "pearson")
  ggplot() +
    geom_abline(slope = 1, color = "gray") +
    geom_point(aes(x = x, y = y)) +
    geom_text_repel(aes(x = x, y = y, label = s), size = 2) +
    xlim(0, max(x, y, na.rm = TRUE)) +
    ylim(0, max(x, y, na.rm = TRUE)) +
    labs(
      title = title, subtitle = paste0("(Pearson r = ", round(r, 2), ")"),
      x = x_lab, y = y_lab
    )
}

comparison_plot(data_total$plot_ba_hs, data_total$homerange_ba_mean, data_total$site, "Basal area [m2/ha]")

comparison_plot(data_total$plot_treeden_gt10cmDbh_hs, data_total$homerange_treeden_gt4in_dbh_mean, data_total$site, "Density (large trees, DBH > 10 cm) [# trees/ha]")

comparison_plot(data_total$plot_treeden_all_hs, data_total$homerange_treeden_all_mean, data_total$site, "Density (all tree sizes) [# trees/ha]")

comparison_plot(data_total$plot_qmd_all_hs, data_total$homerange_qmd_mean, data_total$site, "Quadratic mean diameter [cm]")

comparison_plot(data_total$plot_ht_hs, data_total$homerange_htmax_mean, data_total$site, "Height mean [cm]")

# TODO: plot_ht_cv_hs appears erroneous?
comparison_plot(data_total$plot_ht_cv_hs, data_total$homerange_htmax_cv, data_total$site, "Height CV [cm]")

comparison_plot(data_total$plot_snagden_hs, data_total$homerange_snagden_gt15dbh_mean, data_total$site, "Snag density [cm]")

comparison_plot(data_total$plot_downvol_hs, data_total$homerange_downvol_mean, data_total$site, "Downed wood volume [m3/ha]")

# Inspect specific variables by stage or stratum
ggplot(data_total, aes(x = stage, y = homerange_canopy_layers_mean)) +
  geom_boxplot()

