#### SHARED CONFIG SCRIPT FOR UPDATED OCC VARS

## 3_calc_occurrence_vars.R ###############################################################################
# Quantify variables for occurrence
# NOTE: "hs" refers to data collected from in-person habitat surveys, while "rs" refers to data derived via remote-sensing imagery.
# ETA: 16 hours
#
# CONFIG:
pnts_name = "landscape" # "sites" for surveyed sites only or "landscape" for a grid of points across the landscape
overwrite_data_plot_scale_cache = TRUE
overwrite_data_homerange_scale_cache = TRUE
cover_classification = "clean_strata_4" # e.g. clean_stage_3, strata_4
t = 2020 # year: 2020, 2021, 2022, 2023
#
## OUTPUT:
path_data_plot_scale_out      = paste0("data/cache/3_gis/3_calc_occurrence_vars/V2_data_plot_scale_", t, "_", cover_classification, "_", pnts_name, ".rds")
path_data_homerange_scale_out = paste0("data/cache/3_gis/3_calc_occurrence_vars/V2_data_homerange_scale_", t, "_", cover_classification, "_", pnts_name, ".rds")
#
## INPUT:
path_rast_cover       = paste0("data/cache/3_gis/2_gen_cover_rasters/rast_cover_", t, "_", cover_classification, ".tif")
path_site_cover_class = paste0("data/cache/3_gis/2_gen_cover_rasters/site_cover_class_", t, "_sf.rds")
path_trait_data = "data/cache/2_traits/1_agg_traits/trait_data.csv"
# Base RS-FRIS data (0.1 acre resolution, i.e. ~404m2 or 20.10836 * 20.10836 m grain, roughly 1% of the area of a 100m radius circle)
# RS-FRIS 4.0 uses a combination of 2019 and 2020 photogrammetry.
# RS-FRIS 5.0 uses a combination of 2021 and 2022 photogrammetry. 
###########################################################################################################

source("src/global.R")

message("Calculating occurrence variables for year ", t)
version = rsfris_version_years %>% filter(year == t) %>% pull(version)
rsfris_version_path = paste0("data/environment/rsfris_study_area/", version)

options(mapview.maxpixels = 2117676)

pnts = switch(pnts_name,
              
              # ARU locations
              sites = {
                path_rast_cover = path_rast_cover
                read_rds(path_site_cover_class)
              },
              
              # Existing landscape
              landscape = {
                path_rast_cover = path_rast_cover
                # debug
                landscape_planning_units_pnts = landscape_planning_units_clean %>% filter(unit == "Upper Clearwater")
                bbox = st_bbox(landscape_planning_units_pnts)
                cellsize = 5000 # distance between points (meters)
                grid = st_sf(st_make_grid(x = st_as_sfc(bbox), cellsize = cellsize, offset = c(bbox["xmin"], bbox["ymin"]), what = "centers"))
                grid = grid[st_within(grid, st_union(landscape_planning_units_pnts), sparse = FALSE), ]
                grid = st_sf(geometry = st_geometry(grid))
                # mapview(bbox) + mapview(landscape_planning_units_pnts) + mapview(grid)
              },
              
              stop("Invalid value: ", pnts_name)
)
message("Retrieved ", nrow(pnts), " points for '", pnts_name, "'")

if (pnts_name == "sites") {
  # Ensure that aru site locations are correctly positioned within patches (i.e. not near edges so as to get incorrect covariate estimates)
  # Overwrite thinning class validity for sites: ca263i (mature), ap022i (standinit), az041i (standinit)
  pnts[pnts$site == "ca263i", ]
  pnts[pnts$site == "ap022i", "stratum_4"] = "standinit"
  pnts[pnts$site == "ap022i", "stratum_5"] = "standinit"
  pnts[pnts$site == "az041i", "stratum_4"] = "standinit"
  pnts[pnts$site == "az041i", "stratum_5"] = "standinit"
  
  # Combine co-located sites
  s = site_key %>% filter(site != site_agg) %>% select(site, site_agg) %>% distinct()
  s_site = unique(s$site)
  s_site_agg = unique(s$site_agg)
  mapview(pnts %>% filter(site %in% c(s_site, s_site_agg)))
  lookup = site_key %>% select(site, site_agg) %>% distinct()
  pnts_alt = pnts %>% left_join(lookup, by = "site")
  pnts = pnts_alt %>% filter(site == site_agg) %>% select(-site_agg)
}

message('Loading raster cover data from cache ', path_rast_cover)
rast_cover = rast(path_rast_cover)
# rast_cover = as.factor(rast_cover)

##############################################################################
# Study area, sites, boundaries, and helper functions

species_trait_data = read.csv(path_trait_data)

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
  if (length(x) == 0) return(c(mean = NA, sd = NA, cv = NA))
  x = x * conversion_factor
  med    = median(x, na.rm = na.rm)
  mu     = mean(x, na.rm = na.rm)
  sigma  = sd(x, na.rm = na.rm)
  if (!is.na(mu) && mu == 0) {
    cv = NA
  } else {
    cv = sigma / mu
  }
  return(data.frame(median = med, mean = mu, sd = sigma, cv = cv))
}

##############################################################################
# Rasters and sf objects

rast_origin = load_raster(paste0(rsfris_version_path, '/ORIGIN_YEAR.tif'))
rast_origin_missing = rast_origin
values(rast_origin_missing)[values(rast_origin_missing) <= t] = NA
values(rast_origin)[values(rast_origin) > t] = NA
rast_age = round(t - rast_origin)
