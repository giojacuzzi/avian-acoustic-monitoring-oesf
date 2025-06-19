library(tidyverse)
library(sf)
library(terra)
library(mapview)
options(mapview.maxpixels = 9000000)
library(ggrepel)
theme_set(theme_minimal())
library(landscapemetrics)
library(lwgeom)
library(units)

# TODO: GET RASTER DATED FOR 2020!

### Study area

crs_m = 32610 # EPSG:32610, UTM Zone 10N, (meters)

# ARU locations
aru_sites = st_read('data/environment/GIS Data/AcousticStations.shp') %>%
  st_drop_geometry() %>% janitor::clean_names() %>% as.data.frame() %>%
  st_as_sf(coords = c("utm_e", "utm_n"), crs = crs_m) %>%
  select(name, ces, treatment, geometry) %>% rename(site = name)
mapview(aru_sites, label = aru_sites$name)

# site_data = gpx::read_gpx('data/sites/PAMlocations20230510.gpx')$waypoints %>%
#   janitor::clean_names() %>% rename(site = name) %>% select(site, latitude, longitude)
# site_data = st_transform(st_as_sf(site_data, coords = c('longitude', 'latitude'), crs = crs_data), crs = crs_projected)

# Study area bounding buffer
study_area = st_buffer(st_as_sfc(st_bbox(aru_sites)), dist = 100)

# Strata and forest inventory units

  # Seral stage class
  seral_classes = c("initiation", "canclose", "stemex", "mature")
  class_colors = viridis::viridis(4, option = 'D', direction = -1)
  
  poly_seral_class = st_read('data/environment/GIS Data/Forest Development Strata/AgeStrataFixedDIS_RSFRIS20200130.shp') %>%
    st_transform(crs_m) %>% janitor::clean_names() %>% mutate(class_seral = stratum %>% str_to_lower()) %>%
    mutate(stratum = factor(stratum, levels = seral_classes)) %>% select(-stratum)
  
  # rsfris_units = st_set_crs(st_read('data/environment/GIS Data/Forest Development Strata/AgeStrataFixedRSFRIS20200130.shp'), 'NAD83')
  # poly_rsfris = st_read('data/environment/rs_fris/Polygon_RS-FRIS/Polygon_RS-FRIS.shp')
  # poly_rsfris = st_crop(poly_rsfris, st_transform(study_area, st_crs(poly_rsfris))) %>% st_transform(crs_m) %>% janitor::clean_names() %>% select(riu_id, age, land_cov_c, land_cov_n, geometry)
  
  # Thinning treatment
  poly_thinning_treatment = st_read('data/environment/GIS Data/Forest Development Strata/ThinAfter94NoHarvSinceClipByInitBuf3.shp') %>% 
    st_transform(crs_m) %>% select(TECHNIQUE_, FMA_DT, FMA_STATUS) %>% janitor::clean_names() %>%
    mutate(technique = technique %>% str_to_lower(), fma_status = fma_status %>% str_to_lower()) %>%
    rename(thinning_treatment = technique, thinning_status = fma_status, thinning_date = fma_dt)

mapview(study_area, col.regions = "transparent", color = "black", lwd = 2) +
  mapview(poly_seral_class, col.regions = class_colors) +
  mapview(aru_sites, label = aru_sites$site) + mapview(st_buffer(aru_sites, 100), col.regions = 'transparent', lwd = 2)

# Define patches by seral class, waterbodies, and road intersection (all road types)
patches = poly_seral_class

waterbodies = st_read('data/environment/DNR_Hydrography/DNR_Hydrography_-_Water_Bodies_-_Forest_Practices_Regulation/DNR_Hydrography_-_Water_Bodies_-_Forest_Practices_Regulation.shp')
waterbodies = waterbodies %>% st_crop(st_transform(study_area, st_crs(waterbodies))) %>% st_transform(crs_m) %>% janitor::clean_names()

patches = st_difference(patches, st_union(waterbodies))

# Roads
roads_wadnr = st_read('data/environment/GIS Data/T3roads.shp') %>%
  st_zm(drop = TRUE) %>% janitor::clean_names() %>% st_transform(crs = crs_m) %>% select(road_usgs1, road_surfa, geometry) %>%
  mutate(geometry = st_cast(geometry, "MULTILINESTRING"))

# Get Hoh Mainline Road / Clearwater Road from WSDOT
roads_wsdot = st_read('data/environment/WSDOT_-_Local_Agency_Public_Road_Routes/WSDOT_-_Local_Agency_Public_Road_Routes.shp') %>% filter(RouteIdent %in% c("400000220i", "031265969i")) %>% st_transform(crs_m) %>% select(geometry)
roads_wsdot$road_usgs1 = 'Primary Highway'
roads_wsdot$road_surfa = NA

roads = rbind(roads_wsdot, roads_wadnr)

mapview(roads, zcol = 'road_usgs1') +
  mapview(aru_sites, label = aru_sites$site) + mapview(st_buffer(aru_sites, 100), col.regions = 'transparent', lwd = 2)

roads_highway    = roads %>% filter(road_usgs1 %in% c('Primary Highway'))
roads_lightduty  = roads %>% filter(road_usgs1 %in% c('Light-Duty Road'))
roads_unimproved = roads %>% filter(!(road_usgs1 %in% c('Primary Highway', 'Light-Duty Road')))
roads_highway_buff    = st_buffer(st_union(roads_highway), dist = 5) # half of road width (m)
roads_lightduty_buff  = st_buffer(st_union(roads_lightduty), dist = 3)
roads_unimproved_buff = st_buffer(st_union(roads_unimproved), dist = 2)
roads_buff = st_union(st_union(roads_highway_buff, roads_lightduty_buff), roads_unimproved_buff)

patches = st_difference(patches, st_union(roads_buff))

# Re-id patches as individual polygons (not multipolygons)
patches = patches %>% st_cast("POLYGON")

# "We defined minimum patch size as 0.785 ha and > 50 m wide in the narrowest dimension. This minimum area corresponds roughly to the smallest estimated home range size of any bird species found in the study area (Brown 1985)."
# TODO: Rather than erasing the fragements, consider recursively merging these tiny fragments into the nearest polygon with the smallest area until there are none left
patches = patches %>% mutate(area_m2 = as.numeric(st_area(.)))
bbox_dims = patches %>%
  st_geometry() %>%
  map(~ st_bbox(.x)) %>% map_dfr(~ {
    data.frame(
      width  = .x["xmax"] - .x["xmin"],
      height = .x["ymax"] - .x["ymin"]
    )
  }) %>%
  mutate(min_dim = pmin(width, height))

patches = patches %>%
  bind_cols(bbox_dims) %>%
  filter(area_m2 > 7850, min_dim > 50)

# Add patch id to sites
patches = patches %>% mutate(patch_id = row_number())

site_data = aru_sites
site_data = st_join(site_data, patches[, c("patch_id")])

n_na = sum(is.na(site_data$patch_id))
if (n_na > 0) {
  warning(sprintf("Removing %d sites(s) not located in any patch", n_na))
  print(site_data %>% filter(is.na(patch_id)) %>% pull(site))
  site_data = site_data %>% filter(!is.na(patch_id))
}

# Classify sites based on the seral class and thinning prescription polygons in which they reside
site_data = st_join(site_data, poly_seral_class)
table(site_data$class_seral, useNA = 'ifany')
site_data = st_join(site_data, poly_thinning_treatment)
table(site_data$thinning_treatment, useNA = 'ifany')

# Count number of sites in each patch
point_counts = site_data %>%
  st_drop_geometry() %>%
  count(patch_id, name = "nsites")
patches = patches %>%
  left_join(point_counts, by = "patch_id") %>%
  mutate(nsites = replace_na(nsites, 0))

mapview(study_area, col.regions = 'transparent', color = "black", lwd = 2) +
  mapview(patches, zcol = 'class_seral', col.regions = class_colors, label = patches$patch_id) +
  mapview(aru_sites, label = aru_sites$site) + mapview(st_buffer(aru_sites, 100), col.regions = 'transparent', lwd = 2)

## Helper functions =================================================

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

## Plot data ============================================================

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
plot_data = plot_data %>% group_by(site) %>% summarise(across(everything(), ~ coalesce(.[!is.na(.)][1], NA))) # coalesce duplicate site entries

### Plot scale variables ==============================================================================================

# Elevation (m)

  # Local: TODO

  # Regional: TODO

# Origin (age)

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

# Total basal area (m2/ha)

  # Local: TODO
  plot_data$ba_ha_all

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
  
  # Local: TODO
  plot_data$avg_height_m_all
  plot_data$cv_height_all
  plot_data$avg_hlc_all
  
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
  plot_data$snags_ha
  
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
  plot_data$avg_understory_heights
  plot_data$cv_avg_understory_heights
  plot_data$per_cover_hgf
  plot_data$per_cover_all_shrubs
  plot_data$per_cover_total_understory
  plot_data$shrub_layer_vol
  
  # Regional: TODO
  
# Dead and downed woody material
  
  # Local: TODO
  plot_data$vol_alldown_m3
  
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
  plot_data$ba_ha_psme
  
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

### Patch scale variables ==========================================================================
  
  # Area (ha)
  patches$area_m2 = st_area(patches)
  
  # Core area index (%)
  
    core_edge_buffer = 100 # meters inward from patch edge
    core_areas = st_buffer(patches, dist = -core_edge_buffer)
    patches$core_area = st_area(core_areas)
    patches$core_area_index = as.numeric(patches$core_area) / as.numeric(patches$patch_area) * 100
    
    mapview(core_areas) + mapview(patches, zcol = 'patch_area')
  
  # Edge contrast index (%)
  id = 1741 # mature surrounded by stemex
  id = 891 # init surrounded by mostly stemex
  patches$edge_contrast_index = NA
  patch_ids_to_evaluate = patches[patches$nsites > 0, 'patch_id'] %>% pull(patch_id)
  for (id in patch_ids_to_evaluate) {
    print(id)
    patch = st_make_valid(patches %>% filter(patch_id == id))
    # mapview(patch) + mapview(rast_canopy_cover)
    
    edge_buffer <- 50  # meters
    inner_edge = st_difference(patch, st_buffer(patch, dist = -edge_buffer))
    outer_edge = st_difference(st_buffer(patch, dist = edge_buffer), patch)

    v_inner <- vect(inner_edge)
    v_outer <- vect(outer_edge)
    inner_canopy_cover <- mask(rast_canopy_cover, v_inner)
    outer_canopy_cover <- mask(rast_canopy_cover, v_outer)
    vals_inner <- values(inner_canopy_cover, na.rm = TRUE)
    vals_outer <- values(outer_canopy_cover, na.rm = TRUE)
    mean_inner <- mean(vals_inner)
    mean_outer <- mean(vals_outer)
    
    # Mean edge contrast index
    patches[patches$patch_id == id, 'edge_contrast_index'] = abs(mean_inner - mean_outer)
  }
  mapview(patches, zcol = 'edge_contrast_index') +
    mapview(rast_canopy_cover)

### Landscape scale variables (in a buffered radius) ===============================================

buffer_size = 500
buffer_area = pi * buffer_size^2
  
# Distance to edge (m)

  site_data = site_data %>% mutate(min_dist_m = NA)
  for (i in seq_len(nrow(site_data))) {
    site = site_data[i, ]
    containing_patch = patches %>% filter(patch_id == site$patch_id)
    if (nrow(containing_patch) == 1) {
      site_data[i, 'min_dist_m'] <- st_distance(st_geometry(site), st_boundary(containing_patch))
    } else {
      site_data[i, 'min_dist_m'] <- NA
    }
  }
  
# Proportional abundance of each patch class (%)
  site_data = site_data %>% mutate(
      propab_initiation_500m = NA,
      propab_canclose_500m = NA,
      propab_stemex_500m = NA,
      propab_mature_500m = NA
    )
  for (i in seq_len(nrow(site_data))) {
    site = site_data[i, ]
    site_buff = st_intersection(patches, st_buffer(site, dist = buffer_size))
    site_buff = site_buff %>% mutate(area_m2 = st_area(geometry))
    
    area_by_class = site_buff %>%
      group_by(class_seral) %>%
      summarise(total_area_m2 = sum(area_m2)) %>%
      ungroup()
    area_by_class = area_by_class %>%
      mutate(proportional_abundance = as.numeric(total_area_m2) / buffer_area)
    # mapview(area_by_class, zcol = 'proportional_abundance)
    
    site_data[i, 'propab_initiation_500m'] = area_by_class %>%
      filter(class_seral == 'initiation') %>%
      pull(proportional_abundance) %>% { if (length(.) == 0) 0 else . }
    site_data[i, 'propab_canclose_500m'] = area_by_class %>%
      filter(class_seral == 'canclose') %>%
      pull(proportional_abundance) %>% { if (length(.) == 0) 0 else . }
    site_data[i, 'propab_stemex_500m'] = area_by_class %>%
      filter(class_seral == 'stemex') %>%
      pull(proportional_abundance) %>% { if (length(.) == 0) 0 else . }
    site_data[i, 'propab_mature_500m'] = area_by_class %>%
      filter(class_seral == 'mature') %>%
      pull(proportional_abundance) %>% { if (length(.) == 0) 0 else . }
  }

# Patch diversity/evenness/richness (#)
  
  site_data = site_data %>% mutate(
    richness = NA, diversity_simpson = NA, diversity_shannon = NA, evenness_shannon = NA
  )
  for (i in seq_len(nrow(site_data))) {
    site = site_data[i, ]
    propab = c( # proportional abundances
      site$propab_initiation_500m,
      site$propab_canclose_500m,
      site$propab_stemex_500m,
      site$propab_mature_500m
    )
    propab_nonzero = propab[propab > 0] # avoid invalid products with log(0)
    
    S = length(propab_nonzero) # Richness (count of non-zero forest cover classes)
    
    D = 1 - sum(propab_nonzero^2) # Simpson diversity index
    
    H = -sum(propab_nonzero * log(propab_nonzero)) # Shannon diversity index
    
    E = if (S > 1) H / log(S) else 0 # Shannon evenness index
    
    site_data[i, 'richness'] = S
    site_data[i, 'diversity_simpson'] = D
    site_data[i, 'diversity_shannon'] = H
    site_data[i, 'evenness_shannon'] = E
  }

# Patch similarity/isolation index (#)
  
  # TODO

# Edge density: TODO - incorporate waterbodies as well
  
  site_data = site_data %>% mutate(edges_m_per_m2 = NA)
  for (i in seq_len(nrow(site_data))) {
    print(i)
    site = site_data[i, ]
    site_buff = st_buffer(site, dist = buffer_size)
    site_patches = st_intersection(poly_seral_class, site_buff) # NOTE: poly_seral_class, not patches, to avoid double-counting roads as two edges
    
    # Get patch edges, including roads, but remove artificial edges produced on the buffer boundary
    patch_edges = st_boundary(site_patches)
    patch_edges = st_difference(patch_edges, st_boundary(site_buff))
    road_edges = st_union(st_intersection(roads, site_buff))
    n_patch_edges = length(st_geometry(patch_edges))
    n_road_edges = length(st_geometry(road_edges))
    edges = st_sf(geometry = st_sfc(), crs = crs_m)
    if (n_patch_edges == 0 & n_road_edges > 0) {
      edges = road_edges # contiguous patch
    } else if (n_patch_edges > 0 & n_road_edges == 0) {
      edges = patch_edges
    } else if (n_patch_edges > 0 & n_road_edges > 0) {
      edges = st_union(patch_edges, road_edges)
    }
    
    # mapview(site) + mapview(site_patches) + mapview(edges)

    edge_density = as.numeric(sum(st_length(edges))) / as.numeric(buffer_area) # * 10000 for m/ha
    site_data$edges_m_per_m2[i] = edge_density
  }
  
# TODO: Consider contrast-weighted edge density

# Density of roads (m/m2)
  
  buffer_size = 500
  buffer_area = pi * buffer_size^2
  roads_paved_buffers = st_intersection(roads_paved, st_buffer(site_data, dist = buffer_size))
  
  roads_paved_buffers = roads_paved_buffers %>% # calculate total paved road length per site
    mutate(roads_paved_length = st_length(geometry)) %>%
    group_by(site) %>%
    summarize(roads_paved_total_length_m = sum(roads_paved_length)) %>%
    ungroup()
  roads_paved_buffers = roads_paved_buffers %>% # calculate density of paved roads
    mutate(roads_paved_m_per_m2 = as.numeric(roads_paved_total_length_m) / buffer_area)
  
  site_data = site_data %>% left_join(roads_paved_buffers %>% st_drop_geometry(), by = 'site')
  
  # roads_usgs = st_read('/Users/giojacuzzi/repos/avian-acoustic-monitoring-oesf/data/environment/TRAN_Washington_State_Shape/Shape/Trans_RoadSegment_0.shp')
  # roads_usgs = roads_usgs %>% st_crop(st_transform(study_area, st_crs(roads_usgs))) %>% st_transform(crs_m)

# Density of streams (km/km2)
  
  watercourses = st_read('data/environment/DNR_Hydrography/DNR_Hydrography_-_Watercourses_-_Forest_Practices_Regulation/DNR_Hydrography_-_Watercourses_-_Forest_Practices_Regulation.shp')
  watercourses = watercourses %>% st_crop(st_transform(study_area, st_crs(watercourses))) %>% st_transform(crs_m) %>% janitor::clean_names() %>% select(geometry, wc_cart_1)

  buffer_size = 500
  buffer_area = pi * buffer_size^2
  watercourse_buffers = st_intersection(watercourses, st_buffer(aru_sites, dist = buffer_size))

  watercourse_buffers = watercourse_buffers %>% # calculate total watercourse length per site
    mutate(watercourse_length = st_length(geometry)) %>%
    group_by(site) %>%
    summarize(watercourse_total_length_m = sum(watercourse_length)) %>%
    ungroup()
  watercourse_buffers = watercourse_buffers %>% # calculate density of watercourses
    mutate(watercourse_m_per_m2 = as.numeric(watercourse_total_length_m) / buffer_area)
  
  site_data = site_data %>% left_join(watercourse_buffers %>% st_drop_geometry(), by = 'site')

####################################################################################################
# Local habitat data

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







######## DEBUG

grain = c(20.10835, 20.10835) # 0.1 acre spatial grain (~404.35 m²)
grain = 3 # 2 meter spatial grain

# Define rasterized patches according to seral class and road intersection
rast_seral_class = rasterize(vect(patches), rast(vect(patches), resolution = grain), field = 'class_seral')
check_landscape(rast_seral_class)
lsm_l_np(rast_seral_class) # number of patches
lsm_c_np(rast_seral_class) # number of patches per class

rast_patches = get_patches(rast_seral_class)

p_area_m2 = lsm_p_area(rast_seral_class)


# check for patches of class 1
rast_patches_0 <- get_patches(rast_seral_class, class = 0)[[1]][[1]]
rast_patches_1 <- get_patches(rast_seral_class, class = 1)[[1]][[1]]
rast_patches_2 <- get_patches(rast_seral_class, class = 2)[[1]][[1]]
rast_patches_3 <- get_patches(rast_seral_class, class = 3)[[1]][[1]]

rast_patches = get_patches(rast_seral_class)
mapview(rast_patches$layer_1$class_0, col.regions = 'red') +
  mapview(rast_patches$layer_1$class_1, col.regions = 'blue') +
  mapview(rast_patches$layer_1$class_2, col.regions = 'green') +
  mapview(rast_patches$layer_1$class_3, col.regions = 'yellow')

p <- get_patches(rast_seral_class, class = 1)[[1]][[1]]
length(unique(na.omit(values(p))))
max_val <- max(values(p), na.rm = TRUE)
plot(p, col = rainbow(max_val), main = "Individual patches")


p_area_m2 = lsm_p_area(rast_seral_class)

cores <- calculate_lsm(rast_seral_class, what = "lsm_p_core")

patches_list <- get_patches(rast_seral_class, return_raster = TRUE)

# Example: take the first patch layer (patches for class 1)
patch_raster <- patches_list[[1]][1]

# Apply negative buffer (shrink edges = simulate core area)
core_area <- terra::buffer(patch_raster, width = -100)  # 100 meters inward

plot(core_area, main = "Buffered Core Area - Class 1")

r_cat <- classify(rast_canopy_cover, rbind(c(-Inf, 50, 0), c(50, Inf, 1)), include.lowest = TRUE)