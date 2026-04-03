# 3_calc_occurrence_vars.R ########################################################
# Quantify variables for occurrence
# NOTE: "hs" refers to data collected from in-person habitat surveys,
#       while "rs" refers to data derived via remote-sensing imagery.
#
# OPTIMIZATION STEP 1: Remove dead patch-delineation code            (verified)
# OPTIMIZATION STEP 2: Parallelise per-point loop with mclapply      (verified)
# OPTIMIZATION STEP 3: Remove redundant homerange_and_edge_buffer    (verified)
# OPTIMIZATION STEP 4: Micro-optimizations inside the worker
#   a) Pre-compute all homerange buffer polygons before mclapply
#   b) vect(homerange_buffer) computed once, used for both crop and mask
#   c) ifel / rast_homerange moved inside the plot-scale block
#   d) freq + dplyr pipeline replaced: ncells_homerange derived from freq(),
#      8 filter+pull calls replaced with a single named-vector lookup
#   e) Pointless single-iteration for(k in seq_along(regions)) loop removed

library(parallel)

## CONFIG:
pnts_name                        = "landscape"
overwrite_data_plot_scale_cache  = TRUE
overwrite_data_homerange_scale_cache = TRUE
cover_classification             = "clean_strata_4"
t                                = 2020

### OUTPUT:
path_data_plot_scale_out      = paste0("data/cache/7_landscape/OPT_landscape_data_plot_scale_",      t, "_", cover_classification, "_", pnts_name, ".rds")
path_data_homerange_scale_out = paste0("data/cache/7_landscape/OPT_landscape_data_homerange_scale_", t, "_", cover_classification, "_", pnts_name, ".rds")

### INPUT:
path_rast_cover       = paste0("data/cache/3_gis/2_gen_cover_rasters/rast_cover_",       t, "_", cover_classification, ".tif")
path_site_cover_class = paste0("data/cache/3_gis/2_gen_cover_rasters/site_cover_class_", t, "_sf.rds")
path_trait_data       = "data/cache/2_traits/1_agg_traits/trait_data.csv"
# Base RS-FRIS data (0.1 acre resolution, i.e. ~404m2 or 20.10836 * 20.10836 m grain)
# RS-FRIS 4.0 uses a combination of 2019 and 2020 photogrammetry.
# RS-FRIS 5.0 uses a combination of 2021 and 2022 photogrammetry.
path_msom = "data/cache/models/V4_msom_V4_nofp_nofp_all.rds"
message("Loading data for multi-species occupancy model ", path_msom)
model_data = readRDS(path_msom)
(species = model_data$species)
rm(model_data)

#################################################################################
source("src/global.R")
message("Calculating occurrence variables for year ", t)

version             = rsfris_version_years %>% filter(year == t) %>% pull(version)
rsfris_version_path = paste0("data/environment/rsfris_study_area/", version)
options(mapview.maxpixels = 2117676)

pnts = switch(
  pnts_name,
  landscape = {
    path_rast_cover = path_rast_cover
    landscape_planning_units_pnts = landscape_planning_units_clean %>% filter(unit == "Upper Clearwater")
    bbox     = st_bbox(landscape_planning_units_pnts)
    cellsize = 2500
    grid     = st_sf(st_make_grid(x = st_as_sfc(bbox), cellsize = cellsize,
                                  offset = c(bbox["xmin"], bbox["ymin"]), what = "centers"))
    grid     = grid[st_within(grid, st_union(landscape_planning_units_pnts), sparse = FALSE), ]
    grid     = st_sf(geometry = st_geometry(grid))
  },
  stop("Invalid value: ", pnts_name)
)

message("Retrieved ", nrow(pnts), " points for '", pnts_name, "'")
mapview(landscape_planning_units_clean) + mapview(pnts)

message('Loading raster cover data from cache ', path_rast_cover)
rast_cover = rast(path_rast_cover)

#################################################################################
species_trait_data = read.csv(path_trait_data)

load_raster = function(path_rast) {
  r    = rast(path_rast)
  sa_r = vect(st_transform(study_area, crs(r)))
  r    = crop(r, sa_r)
  r    = mask(r, sa_r)
  r    = project(r, crs_m_rast)
  return(r)
}

compute_raster_buffer_value_func = function(raster, points, buffer_distance, func) {
  r             = project(raster, crs(points))
  region        = st_buffer(points, dist = buffer_distance)
  buffer_values = terra::extract(r, vect(region), fun = func, na.rm = TRUE)
}

summary_stats = function(x, na.rm = TRUE, conversion_factor = 1) {
  if (length(x) == 0) return(c(mean = NA))
  x  = x * conversion_factor
  mu = mean(x, na.rm = na.rm)
  return(data.frame(mean = mu))
}

#################################################################################
rast_origin         = load_raster(paste0(rsfris_version_path, '/ORIGIN_YEAR.tif'))
rast_origin_missing = rast_origin
values(rast_origin_missing)[values(rast_origin_missing) <= t] = NA
values(rast_origin)[values(rast_origin) > t]                  = NA
rast_age = round(t - rast_origin)

#################################################################################
if (overwrite_data_homerange_scale_cache) {
  
  message("Generating homerange scale data cache (current time ", time_start <- Sys.time(), ")")
  data_homerange_scale = list()
  
  message("Homerange buffer min: ",    min(species_trait_data$home_range_radius_m))
  message("Homerange buffer median: ", buff_median <- round(median(species_trait_data$home_range_radius_m), 0))
  message("Homerange buffer mean: ",   buff_mean   <- round(mean(species_trait_data$home_range_radius_m),   0))
  message("Homerange buffer max: ",    buff_max    <- round(max(species_trait_data$home_range_radius_m),    0))
  
  homeranges = species_trait_data %>%
    select(common_name, home_range_radius_m) %>%
    filter(common_name %in% species) %>%
    mutate(home_range_radius_m = round(home_range_radius_m, 0)) %>%
    rename(scale = common_name)
  
  homeranges = rbind(
    data.frame(scale = c("plot", "median"), home_range_radius_m = c(100, buff_median)),
    homeranges
  )
  
  message("Loading forest structure rasters")
  rast_rsfris = list(
    rast_BAP_HWD     = load_raster(paste0(rsfris_version_path, '/BAP_HWD.tif')),
    rast_QMD_6       = load_raster(paste0(rsfris_version_path, '/QMD_6.tif')),
    rast_TREE_ACRE_6 = load_raster(paste0(rsfris_version_path, '/TREE_ACRE_6.tif'))
  )
  
  n_cores = max(1L, detectCores() - 1L)
  message("Parallelising over ", n_cores, " cores (", nrow(pnts), " points x ",
          nrow(homeranges), " scales)")
  
  all_cover_classes = c("thin", "standinit", "compex", "underdev", "old", "mature", "road_paved", "water")
  
  message("Calculating homerange variables across all scales")
  for (i in 1:nrow(homeranges)) {
    scale                 = homeranges[i, 'scale']
    homerange_buffer_size = homeranges[i, 'home_range_radius_m']
    
    if (homerange_buffer_size < 100) {
      message("Skipping scale '", scale, "', home range buffer size: ", homerange_buffer_size)
      next
    }
    
    message("Scale '", scale, "', home range buffer size: ", homerange_buffer_size)
    
    # 4a: pre-compute all buffer polygons once per scale, outside mclapply.
    # Workers index into this object rather than calling st_buffer themselves.
    homerange_buffers = st_buffer(pnts, homerange_buffer_size)
    
    data_homerange_scale_species_list = mclapply(
      seq_len(nrow(pnts)),
      mc.cores = n_cores,
      FUN = function(j) {
        
        homerange_buffer = homerange_buffers[j, ]
        
        # 4b: convert to SpatVector once; reuse for both crop and mask
        homerange_buffer_vect = vect(homerange_buffer)
        
        homerange_crop = tryCatch({
          crop(rast_cover, homerange_buffer_vect)
        }, error = function(e) {
          if (grepl("extents do not overlap", e$message)) return(NULL)
        })
        
        if (is.null(homerange_crop)) return(NULL)
        
        homerange_cover = mask(homerange_crop, homerange_buffer_vect)
        
        ## Structural variables (plot scale only)
        if (scale == "plot") {
          # 4c: ifel / rast_homerange only created when actually needed
          rast_homerange = ifel(!is.na(homerange_cover), 1, NA)
          # 4e: removed the single-element for(k in seq_along(regions)) loop
          region     = resample(rast_homerange, rast_age, method = "near")
          rsfris_means        = sapply(rast_rsfris, function(r) mean(values(mask(r, region)), na.rm = TRUE))
          names(rsfris_means) = paste0(sub("^rast_", "", names(rsfris_means)), "_mean")
          plot_rsfris_df      = clean_names(as.data.frame(t(rsfris_means)))
        }
        
        ## Homerange scale metrics
        # 4d: freq() called once; ncells derived from it, eliminating the
        #     separate sum(!is.na(values())) call.  Named-vector lookup replaces
        #     8 individual filter+pull operations.
        freq_df          = freq(homerange_cover)
        ncells_homerange = sum(freq_df$count)
        cover_counts     = setNames(freq_df$count, freq_df$value)
        
        get_pcnt = function(cl) {
          if (cl %in% names(cover_counts)) cover_counts[[cl]] / ncells_homerange else 0
        }
        
        final_df = data.frame(
          scale           = scale,
          buffer_radius_m = homerange_buffer_size,
          pcnt_standinit  = get_pcnt("standinit"),
          pcnt_compex     = get_pcnt("compex"),
          # pcnt_underdev   = get_pcnt("underdev"),
          # pcnt_old        = get_pcnt("old"),
          pcnt_mature     = get_pcnt("mature"),
          pcnt_thin       = get_pcnt("thin")
          # pcnt_road_paved = get_pcnt("road_paved"),
          # pcnt_water      = get_pcnt("water")
        )
        
        if (scale == "plot") final_df = cbind(final_df, plot_rsfris_df)
        if (pnts_name == "sites") final_df$site = pnts[j, ] %>% pull(site)
        
        final_df
      }
    )
    
    data_homerange_scale_species = bind_rows(Filter(Negate(is.null), data_homerange_scale_species_list))
    
    if (FALSE) {
      site_check        = data_homerange_scale_species %>% slice_max(focalpatch_isolation, n = 1, with_ties = FALSE) %>% pull(site)
      site_check        = data_plot_scale %>% filter(site == site_check)
      site_check_buffer = st_buffer(site_check, homerange_buffer_size)
      mapview(site_check) + mapview(mask(crop(rast_cover, site_check_buffer), site_check_buffer))
      
      site_check        = data_homerange_scale_species %>% slice_min(focalpatch_isolation, n = 1, with_ties = FALSE) %>% pull(site)
      site_check        = data_plot_scale %>% filter(site == site_check)
      site_check_buffer = st_buffer(site_check, homerange_buffer_size)
      mapview(site_check) + mapview(mask(crop(rast_cover, site_check_buffer), site_check_buffer))
    }
    
    data_homerange_scale[[scale]] = data_homerange_scale_species
  }
  
  dir.create(dirname(path_data_homerange_scale_out), recursive = TRUE, showWarnings = FALSE)
  saveRDS(data_homerange_scale, path_data_homerange_scale_out)
  message(crayon::green("Cached homerange data cache to", path_data_homerange_scale_out,
                        "(", round(as.numeric(difftime(Sys.time(), time_start, units = 'mins')), 2), "min )"))
  
} else {
  message('Loading homerange scale data from cache ', path_data_homerange_scale_out)
  data_homerange_scale = readRDS(path_data_homerange_scale_out)
}