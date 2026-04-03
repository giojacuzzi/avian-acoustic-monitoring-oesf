# DEBUG: trying to speed up homerange variable calc

# 3_calc_occurrence_vars.R ###############################################################################
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
path_data_plot_scale_out      = paste0("data/cache/7_landscape/landscape_data_plot_scale_", t, "_", cover_classification, "_", pnts_name, ".rds")
path_data_homerange_scale_out = paste0("data/cache/7_landscape/landscape_data_homerange_scale_", t, "_", cover_classification, "_", pnts_name, ".rds")
#
## INPUT:
path_rast_cover       = paste0("data/cache/3_gis/2_gen_cover_rasters/rast_cover_", t, "_", cover_classification, ".tif")
path_site_cover_class = paste0("data/cache/3_gis/2_gen_cover_rasters/site_cover_class_", t, "_sf.rds")
path_trait_data = "data/cache/2_traits/1_agg_traits/trait_data.csv"
# Base RS-FRIS data (0.1 acre resolution, i.e. ~404m2 or 20.10836 * 20.10836 m grain, roughly 1% of the area of a 100m radius circle)
# RS-FRIS 4.0 uses a combination of 2019 and 2020 photogrammetry.
# RS-FRIS 5.0 uses a combination of 2021 and 2022 photogrammetry. 

path_msom = "data/cache/models/V4_msom_V4_nofp_nofp_all.rds"
message("Loading data for multi-species occupancy model ", path_msom)
model_data = readRDS(path_msom)
(species = model_data$species)
rm(model_data)

###########################################################################################################

source("src/global.R")

message("Calculating occurrence variables for year ", t)
version = rsfris_version_years %>% filter(year == t) %>% pull(version)
rsfris_version_path = paste0("data/environment/rsfris_study_area/", version)

options(mapview.maxpixels = 2117676)

pnts = switch(pnts_name,
              
              # Existing landscape
              landscape = {
                path_rast_cover = path_rast_cover
                # debug
                landscape_planning_units_pnts = landscape_planning_units_clean %>% filter(unit == "Upper Clearwater")
                bbox = st_bbox(landscape_planning_units_pnts)
                cellsize = 2500 # distance between points (meters)
                grid = st_sf(st_make_grid(x = st_as_sfc(bbox), cellsize = cellsize, offset = c(bbox["xmin"], bbox["ymin"]), what = "centers"))
                grid = grid[st_within(grid, st_union(landscape_planning_units_pnts), sparse = FALSE), ]
                grid = st_sf(geometry = st_geometry(grid))
                # mapview(bbox) + mapview(landscape_planning_units_pnts) + mapview(grid)
              },
              
              stop("Invalid value: ", pnts_name)
)
message("Retrieved ", nrow(pnts), " points for '", pnts_name, "'") 

mapview(landscape_planning_units_clean) + mapview(pnts)

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
summary_stats = function(x, na.rm = TRUE, conversion_factor = 1) {
  if (length(x) == 0) return(c(
    mean = NA
  ))
  x = x * conversion_factor
  mu     = mean(x, na.rm = na.rm)
  return(data.frame(
    mean = mu
  ))
}

##############################################################################
# Rasters and sf objects

rast_origin = load_raster(paste0(rsfris_version_path, '/ORIGIN_YEAR.tif'))
rast_origin_missing = rast_origin
values(rast_origin_missing)[values(rast_origin_missing) <= t] = NA
values(rast_origin)[values(rast_origin) > t] = NA
rast_age = round(t - rast_origin)

##############################################################################
# Focal patch scale and home range neighborhood scale covariates (limited to the extent of the species home range)

if (overwrite_data_homerange_scale_cache) {
  
  message("Generating homerange scale data cache (current time ", time_start <- Sys.time(), ")")
  
  data_homerange_scale = list()
  
  message("Homerange buffer min: ",    min(species_trait_data$home_range_radius_m))
  message("Homerange buffer median: ", buff_median <- round(median(species_trait_data$home_range_radius_m),0))
  message("Homerange buffer mean: ",   buff_mean   <- round(mean(species_trait_data$home_range_radius_m),0))
  message("Homerange buffer max: ",    buff_max    <- round(max(species_trait_data$home_range_radius_m),0))
  # Calculate species-specific ranges
  homeranges = species_trait_data %>% select(common_name, home_range_radius_m) %>%
    filter(common_name %in% species) %>%
    mutate(home_range_radius_m = round(home_range_radius_m, 0)) %>% rename(scale = common_name)
  homeranges = rbind(data.frame(
    scale = c("plot", "median"), # can also do "mean", "max"
    home_range_radius_m = c(100, buff_median) # buff_mean, buff_max
  ), homeranges)

  ## Load forest structure rasters
  message("Loading forest structure rasters")
  
  rast_rsfris = list(
    # BAP_HWD - Percent of trees of a hardwood species [%]
    rast_BAP_HWD = load_raster(paste0(rsfris_version_path, '/BAP_HWD.tif')),
    # QMD_6 - Quadratic mean diameter of trees > 6" DBH [in]-->[cm]
    rast_QMD_6 = load_raster(paste0(rsfris_version_path, '/QMD_6.tif')),
    # TREE_ACRE_6 - Number of trees per acre > 6" DBH [#/acre] --> [#/ha]
    rast_TREE_ACRE_6 = load_raster(paste0(rsfris_version_path, '/TREE_ACRE_6.tif'))
  )
  
  # "Edge influences on distributions of organisms or factors affecting organisms (such as predation and nest parasitism) are concentrated within 50m of the edge." (Kremaster L. & Bunnell F. L. (1999) Edge effects: Theory, evidence and implications to management of western North American forests. In: Forest Fragmentation: Wildlife and Management Implications (eds J. Wisniewski, J. A. Rochelle & L. Lehmann) pp. 117–53. Leiden, Boston, MA.)
  core_area_buffer = 50 # meters
  core_area_buffer_ncells = round(core_area_buffer / res(rast_cover)[1], 0)
  
  # Calculate per scale according to predicted home range size
  message("Calculating homerange variables across all scales")
  for (i in 1:nrow(homeranges)) {
    scale = homeranges[i,'scale']
    homerange_buffer_size = homeranges[i,'home_range_radius_m']
    
    if (homerange_buffer_size < 100) {
      message("Skipping scale '", scale, "', home range buffer size: ", homerange_buffer_size)
      next
    }
    
    message("Scale '", scale, "', home range buffer size: ", homerange_buffer_size)
    
    data_homerange_scale_species = data.frame()
    data_homerange_scale_species_list = list()
    
    # For each point...
    pb = progress_bar$new(format = "[:bar] :percent :elapsedfull (ETA :eta)", total = nrow(pnts), clear = FALSE)
    
    # Rprof("my_latest_profile.out")
    for (j in 1:nrow(pnts)) {
      
      ## Get cover data within and surrounding homerange buffer
      pnt = pnts[j,]
      
      homerange_and_edge_buffer = st_buffer(pnt, homerange_buffer_size + core_area_buffer) # additional buffer to ensure that edges are retained in mask
      homerange_buffer = st_buffer(pnt, homerange_buffer_size)
      homerange_and_edge_crop = tryCatch({
        crop(rast_cover, vect(homerange_and_edge_buffer))
      }, error = function(e) {
        if (grepl("extents do not overlap", e$message)) {
          # warning(paste("Skipping point", j, "- no raster overlap"))
          return(NULL)
        }
      })
      
      # Skip to next iteration if crop failed
      if (is.null(homerange_and_edge_crop)) next
      
      homerange_and_edge = mask(homerange_and_edge_crop, vect(homerange_and_edge_buffer))
      homerange_crop = crop(rast_cover, vect(homerange_buffer))
      homerange_cover = mask(homerange_crop, vect(homerange_buffer))
      rast_homerange = ifel(!is.na(homerange_cover), 1, NA)
      # mapview(homerange_cover) + mapview(homerange_and_edge) + mapview(homerange_buffer)
      
      hr_patch_ids = homerange_and_edge
      values(hr_patch_ids) <- NA
      start_id = 1  # Starting patch ID to ensure uniqueness
      cls = sort(na.omit(unique(values(homerange_and_edge))))
      for (cl in cls) {
        # message("Class: ", cl)
        r_class = mask(homerange_and_edge, homerange_and_edge == cl, maskvalues = FALSE) # Mask only the current class
        r_patches = patches(r_class, directions = 8) # Compute patches for this class only
        # Reclassify patch IDs to be globally unique
        max_patch_id = global(r_patches, "max", na.rm = TRUE)[1,1]
        if (!is.na(max_patch_id)) {
          r_patches = classify(r_patches, cbind(1:max_patch_id, start_id:(start_id + max_patch_id - 1)))
          start_id = start_id + max_patch_id
        }
        hr_patch_ids = cover(hr_patch_ids, r_patches)
      }
      pid = terra::extract(hr_patch_ids, vect(pnt))[,2]
      rast_focal_patch = trim(classify(hr_patch_ids, cbind(pid, 1), others = NA))
      rast_focal_patch = crop(rast_focal_patch, vect(homerange_buffer))
      rast_focal_patch = mask(rast_focal_patch, vect(homerange_buffer))
      # mapview(rast_focal_patch) + mapview(homerange_cover)
      
      ## Focal patch scale
      
      # Focal patch cover class
      # (focalpatch_cover = terra::extract(rast_cover, vect(pnt))[,2])
      
      # Focal patch area [ha]
      ncells_focal_patch = sum(values(rast_focal_patch), na.rm = TRUE)
      (focal_patch_area = ncells_focal_patch * prod(res(rast_focal_patch)) * conv_m2_to_ha)

      # Structural variables (across plot)
      if (scale == "plot") {
        message(yellow("NOTE: Calculating structrual RSFRIS variables for only the plot scale"))
        results_list = list()
        regions = list(
          # focalpatch = rast_focal_patch,
          homerange  = rast_homerange
        )
        for (k in seq_along(regions)) {
          region_name = names(regions)[k]
          # message("Calculating forest structural variables for: ", region_name)
          region = resample(regions[[k]], rast_age, method = "near") # align region to raster data
          # print(sum(!is.na(values(region))))
          
          # Calculate summary stats for all rsfris rasters
          vars = lapply(
            rast_rsfris,
            function(r) summary_stats(values(mask(r, region)), na.rm = TRUE)
          )
          names(vars) = sub("^rast_", "", names(vars))
          
          # format and discard irrelevant variables
          vars = as.data.frame(vars)
          results_list[[region_name]] = vars
        }
        plot_rsfris_df = do.call(cbind, unname(results_list)) %>% clean_names()
      }
      
      ## Homerange scale metrics
      
      # message("Calculating homerange configuration and composition variables")
      
      # homerange_cover_forest = homerange_cover
      # homerange_cover_forest[!(homerange_cover_forest[] %in% c(1, 2, 3, 4, 5))] = NA
      # mapview(homerange_cover) + mapview(homerange_cover_forest)
      ncells_homerange = sum(!is.na(values(homerange_cover)))
      
      cover_freq = tibble(value = c("thin", "standinit", "compex", "underdev", "old", "mature", "road_paved", "water"), count = 0)
      cover_freq = cover_freq %>% rows_update(freq(homerange_cover) %>% select(value, count), by = "value", unmatched = "ignore")
      
      # Proportional abundance of each cover class [%]
      (pcnt_standinit  = cover_freq %>% filter(value == "standinit")  %>% pull(count) / ncells_homerange) # "stand initiation" 0-25 yr
      (pcnt_compex     = cover_freq %>% filter(value == "compex")     %>% pull(count) / ncells_homerange) # "stem exclusion" 25-80 yr
      (pcnt_underdev   = cover_freq %>% filter(value == "underdev")   %>% pull(count) / ncells_homerange) # "understory reinitiation" 80-200 yr
      (pcnt_old        = cover_freq %>% filter(value == "old")        %>% pull(count) / ncells_homerange) # "old-growth forest" 200+ yr
      (pcnt_mature     = cover_freq %>% filter(value == "mature")     %>% pull(count) / ncells_homerange) # "late-successional and old growth forest" (80+ yr)
      (pcnt_thin       = cover_freq %>% filter(value == "thin")       %>% pull(count) / ncells_homerange) # "commercial thinning"
      (pcnt_road_paved = cover_freq %>% filter(value == "road_paved") %>% pull(count) / ncells_homerange) # "roads"
      (pcnt_water      = cover_freq %>% filter(value == "water")      %>% pull(count) / ncells_homerange) # "water"

      final_df = data.frame(
        scale = scale,
        buffer_radius_m = homerange_buffer_size,
        pcnt_standinit,
        pcnt_compex,
        # pcnt_underdev,
        # pcnt_old,
        pcnt_mature,
        pcnt_thin
        # pcnt_road_paved,
        # pcnt_water
      )
      
      if (scale == "plot") {
        final_df = cbind(final_df, plot_rsfris_df)
      }
      
      if (pnts_name == "sites") {
        final_df$site = pnts[j, ] %>% pull(site)
      }
      
      data_homerange_scale_species_list[[j]] = final_df
      pb$tick()
    }
    data_homerange_scale_species = bind_rows(data_homerange_scale_species_list)
    # Rprof(NULL)
    # summaryRprof("my_profile.out")
    
    # Sanity check results
    if (FALSE) {
      site_check = data_homerange_scale_species %>% slice_max(focalpatch_isolation, n = 1, with_ties = FALSE) %>% pull(site)
      site_check = data_plot_scale %>% filter(site == site_check)
      site_check_buffer = st_buffer(site_check, homerange_buffer_size)
      mapview(site_check) + mapview(mask(crop(rast_cover, site_check_buffer), site_check_buffer))
      
      site_check = data_homerange_scale_species %>% slice_min(focalpatch_isolation, n = 1, with_ties = FALSE) %>% pull(site)
      site_check = data_plot_scale %>% filter(site == site_check)
      site_check_buffer = st_buffer(site_check, homerange_buffer_size)
      mapview(site_check) + mapview(mask(crop(rast_cover, site_check_buffer), site_check_buffer))
    }
    
    data_homerange_scale[[scale]] = data_homerange_scale_species
  }
  dir.create(dirname(path_data_homerange_scale_out), recursive = TRUE, showWarnings = FALSE)
  saveRDS(data_homerange_scale, path_data_homerange_scale_out)
  message(crayon::green("Cached homerange data cache to", path_data_homerange_scale_out, "(", round(as.numeric(difftime(Sys.time(), time_start, units = 'mins')), 2), "min )"))
  
} else { # overwrite_data_homerange_scale_cache is FALSE
  message('Loading homerange scale data from cache ', path_data_homerange_scale_out)
  data_homerange_scale = readRDS(path_data_homerange_scale_out)
}

