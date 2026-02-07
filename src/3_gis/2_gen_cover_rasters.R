## 2_gen_cover_rasters.R #########################################################################################
# Generate and inspect raster(s) for landscape cover
#
# CONFIG:
overwrite_rast_cover_cache = FALSE
#
## OUTPUT:
path_rast_cover_clean_out = "data/cache/occurrence_covariates/rast_cover_clean.tif"
#
## INPUT:
# Base RS-FRIS data (0.1 acre resolution, i.e. ~404m2 or 20.10836 * 20.10836 m grain, roughly 1% of the area of a 100m radius circle)
# RS-FRIS 4.0 uses a combination of 2019 and 2020 photogrammetry.
# RS-FRIS 5.0 uses a combination of 2021 and 2022 photogrammetry. 
# TODO: Calculate on a yearly basis!
dir_rsfris_version = 'data/environment/rsfris_v4.0' # Only use 2020 for now
###########################################################################################################

source("src/global.R")
source("src/3_gis/1_preprocess_gis_data.R")

options(mapview.maxpixels = 2117676)

##############################################################################
# Study area, sites, boundaries, and helper functions

# ARU locations (sites)
# Used to define study area
mapview(aru_sites, label = aru_sites$name)

# Load a raster, cropped, projected, and masked to the study area
load_raster = function(path_rast) {
  r = rast(path_rast)
  sa_r = vect(st_transform(study_area, crs(r)))
  r = crop(r, sa_r)
  r = mask(r, sa_r)
  r = project(r, crs_m_rast)
  return(r)
}

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

template = rast(ext(rast_stratum), resolution = res(rast_stratum), crs = crs(rast_stratum))
watercourse_half_width = res(rast_stratum)[1] / 2
watercourses_buffered = (st_buffer(boundary_watercourses, dist = watercourse_half_width))
rast_watercourses = rasterize(vect(watercourses_buffered), template, field = 7, background = NA, touches=TRUE)
# mapview(rast_watercourses) + mapview(watercourses_buffered)

rast_updated = rast_stratum
rast_updated = cover(rast_watercourses, rast_updated)

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

# "The Hoh-Clearwater Mainline is our only double-lane mainline, and it is about 26 feet wide for the asphalt surface. If we’re looking at right of way widths, we could easily assume a minimum width of about 50 feet for a 12-foot road, probably a good 60-80 feet for a 14-20 foot wide road, and about 100 feet for the Hoh Mainline. 100 feet might be good for US 101 as well for right of way width, and maybe about 30 feet for actual road surface width."
road_half_width_primary = max(100 / 2 * conv_ft_to_m, res(rast_stratum)[1] / 2)
# mapview(paved_primary_roads) + mapview(st_union(st_buffer(paved_primary_roads, road_half_width_primary)))
# "Our minimum road surface width is going to be 12 feet. That would cover most of our roads. Main arterials/single lane mainlines will have a minimum surface width of 14 feet, with a few up to 20 feet."
paved_secondary_roads = st_make_valid(roads %>% filter(road_usgs1 %in% c("Light-Duty Road")))
road_half_width_secondary = max(30 / 2 * conv_ft_to_m, res(rast_stratum)[1] / 2)
# mapview(paved_secondary_roads) + mapview(st_union(st_buffer(paved_secondary_roads, road_half_width_secondary)))

paved_roads_primary_buffered = (st_buffer(paved_primary_roads, dist = road_half_width_primary))
paved_roads_secondary_buffered = (st_buffer(paved_secondary_roads, dist = road_half_width_secondary))
rast_paved_roads_primary = rasterize(vect(paved_roads_primary_buffered), template, field = 6, background = NA, touches=TRUE)
rast_paved_roads_secondary = rasterize(vect(paved_roads_secondary_buffered), template, field = 6, background = NA, touches=TRUE)
# mapview(rast_paved_roads_primary) + mapview(paved_roads_primary_buffered)
# mapview(rast_paved_roads_secondary) + mapview(paved_roads_secondary_buffered)
# mapview(rast_paved_roads_primary) + mapview(rast_paved_roads_secondary)
rast_updated = cover(rast_paved_roads_primary, rast_updated)
rast_updated = cover(rast_paved_roads_secondary, rast_updated)

# TODO: Cover impervious surfaces

# Impute any remaining missing patches with cover class 2 (stem exclusion)
# rast_updated[is.na(rast_updated[])] <- 2
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
  
  # TODO: Save an alternative that combines understory reinit and old growth into an LSOG category
  
  # TODO: Save alternative patch ids too
  
  message('Saving raster cover data cache ', path_rast_cover_clean_out)
  dir.create(dirname(path_rast_cover_clean_out), recursive = TRUE, showWarnings = FALSE)
  writeRaster(rast_cover_clean, path_rast_cover_clean_out, overwrite=TRUE)
  message(crayon::green("Finished generating cover cache (", round(as.numeric(difftime(Sys.time(), time_start, units = 'mins')), 2), " minutes)"))
  
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

