##############################################################################
# Quantify covariates on occurrence
#
# Note "hs" refers to data collected from in-person habitat surveys, while
# "rs" refers to data derived via remote-sensing imagery.
##############################################################################

source("src/global.R")
source("src/1NEW_preprocess_habitat_data.R")

# TODO: Ensure that aru site locations are correctly positioned within patches (i.e. not near edges so as to get incorrect covariate estimates)

# Inputs
path_rast_cover_clean = "data/cache/occurrence_covariates/rast_cover_clean.tif"
pnts_name = "sites" # e.g. "sites" for ARU locations, ""
pnts = { # sf collection or raster of locations to calculate habitat variables for
  aru_sites
}

# Outputs
overwrite_data_plot_scale_cache = TRUE
overwrite_data_homerange_scale_cache = TRUE
path_data_plot_scale_out = paste0("data/cache/habitat_vars/data_plot_scale_", pnts_name, ".rds")
path_data_homerange_scale_out = paste0("data/cache/habitat_vars/data_homerange_scale_", pnts_name, ".rds")

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

species_trait_data = read.csv("data/cache/species_traits/species_traits.csv")

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

# Base RS-FRIS data (0.1 acre resolution, i.e. ~404m2 or 20.10836 * 20.10836 m grain, roughly 1% of the area of a 100m radius circle)
# RS-FRIS 4.0 uses a combination of 2019 and 2020 photogrammetry.
# RS-FRIS 5.0 uses a combination of 2021 and 2022 photogrammetry. 
# TODO: Calculate on a yearly basis!
dir_rsfris_version = 'data/environment/rsfris_v4.0' # Only use 2020 for now

message('Loading raster cover data from cache ', path_rast_cover_clean)
rast_cover_clean = rast(path_rast_cover_clean)
rast_cover_clean = as.factor(rast_cover_clean)

##############################################################################
# Rasters

rast_origin = load_raster(paste0(dir_rsfris_version, '/RS_FRIS_ORIGIN_YEAR.img'))
rast_origin_missing = rast_origin
values(rast_origin_missing)[values(rast_origin_missing) <= 2020] = NA
values(rast_origin)[values(rast_origin) > 2020] = NA
rast_age = round(2020 - rast_origin)

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
    mutate(tag = str_extract(site, "_.*$") %>% str_remove("^_"), site = str_remove(site, "_.*$")) %>% mutate(site = tolower(site))
  
  # Mark sites that were surveyed for habitat data in-person
  if (pnts_name == "sites") {
    intersect(pnts$site, data_plot$site)
    length(intersect(pnts$site, data_plot$site))
    setdiff(data_plot$site, pnts$site) # This should be 0
    pnts$hs = FALSE
    pnts[pnts$site %in% data_plot$site, 'hs'] = TRUE
    mapview(pnts, zcol = 'hs')
  }
  
  # Elevation [m]
  pnts = elevatr::get_elev_point(pnts) %>% rename(elev = elevation)
  
  # Age (mean and cv) [#]
  pnts$age_point = terra::extract(rast_age, vect(pnts))[, 2]
  stats_age      = compute_raster_buffer_value_func(rast_age, pnts, plot_buffer, summary_stats)
  pnts$age_mean = as.numeric(stats_age[,2])
  pnts$age_cv   = as.numeric(stats_age[,4])
  hist(pnts$age_point, breaks = seq(0, max(pnts$age_point) + 10, by = 10))
  hist(pnts$age_mean, breaks = seq(0, max(pnts$age_mean) + 10, by = 10))
  table(pnts$age_point)
  table(round(pnts$age_mean))
  
  # # Cover class [categorical]
  # pnts$stage
  
  # # Thinning treatment [categorical]
  # pnts$thinning_treatment
  
  # Basal area (all live trees) [m2/ha] TODO: confirm if this is a mean value at local plot level
  pnts = pnts %>% left_join(data_plot %>% select(site, plot_ba_hs = ba_ha_all), by = 'site')
  
  # Tree density (all live trees large, dbh > 10 cm) [# trees/ha]
  pnts = pnts %>% left_join(data_plot %>% select(site, plot_treeden_gt10cmDbh_hs = large_per_hectare_all), by = 'site')
  
  # Tree density (small, dbh < 10 cm) [# trees/ha]
  pnts = pnts %>% left_join(data_plot %>% select(site, plot_treeden_lt10cmDbh_hs = small_per_hectare), by = 'site')
  
  # Total tree density (all sizes) [# trees/ha]
  pnts$plot_treeden_all_hs = pnts$plot_treeden_gt10cmDbh_hs + pnts$plot_treeden_lt10cmDbh_hs
  
  # Tree quadratic mean diameter (large, dbh > 10cm) [cm]
  pnts = pnts %>% left_join(
    data_plot %>% select(site, plot_qmd_gt10cmDbh_hs = avg_dbh_cm_all) %>%
      mutate(plot_qmd_gt10cmDbh_hs = replace_na(plot_qmd_gt10cmDbh_hs, 0.0)),
    by = 'site')
  # Tree quadratic mean diameter (small, dbh < 10 cm) [cm]
  pnts = pnts %>% left_join(
    data_plot %>% select(site, plot_qmd_lt10cmDbh_hs = avg_dbh_cm) %>%
      mutate(plot_qmd_lt10cmDbh_hs = replace_na(plot_qmd_lt10cmDbh_hs, 0.0)),
    by = 'site')
  # Tree quadratic mean diameter (all) [cm]
  pnts$plot_qmd_all_hs = pnts$plot_qmd_gt10cmDbh_hs + pnts$plot_qmd_lt10cmDbh_hs
  
  # Tree height (mean and cv) [m]
  # TODO: Get high-resolution LiDAR height data from DNR approval procedure
  pnts = pnts %>% left_join(
    data_plot %>% select(site, plot_ht_hs = avg_height_m_all) %>%
      mutate(plot_ht_hs = replace_na(plot_ht_hs, 0.0)),
    by = 'site')
  pnts = pnts %>% left_join(
    data_plot %>% select(site, plot_ht_cv_hs = cv_height_all) %>%
      mutate(plot_ht_cv_hs = replace_na(plot_ht_cv_hs, 0.0)),
    by = 'site')
  
  # Height to live crown [m], length of live crown (live crown depth) [m], and live crown ratio [#]
  pnts = pnts %>% left_join(
    data_plot %>% select(site, plot_hlc_hs = avg_hlc_all) %>%
      mutate(plot_hlc_hs = replace_na(plot_hlc_hs, 0.0)), by = 'site')
  pnts = pnts %>% left_join(
    data_plot %>% select(site, plot_llc_hs = avg_llc_all) %>%
      mutate(plot_llc_hs = replace_na(plot_llc_hs, 0.0)), by = 'site')
  pnts = pnts %>% left_join(
    data_plot %>% select(site, plot_lcr_hs = avg_lcr_all) %>%
      mutate(plot_lcr_hs = replace_na(plot_lcr_hs, 0.0)), by = 'site')
  
  # Density of all snags [# snags/ha]
  pnts = pnts %>% left_join(data_plot %>% select(site, plot_snagden_hs = snags_ha) %>% mutate(plot_snagden_hs = replace_na(plot_snagden_hs, 0.0)), by = 'site')
  
  # Downed wood volume [m3/ha]
  pnts = pnts %>% left_join(data_plot %>% select(site, plot_downvol_hs = vol_alldown_m3) %>% mutate(plot_downvol_hs = replace_na(plot_downvol_hs, 0.0)), by = 'site')
  
  # Understory vegetation cover [%] and volume [m3/ha]
  pnts = pnts %>% left_join(data_plot %>% select(site, plot_understory_cover = per_cover_total_understory) %>% mutate(plot_understory_cover = replace_na(plot_understory_cover, 0.0)), by = 'site')
  pnts = pnts %>% left_join(data_plot %>% select(site, plot_understory_vol = shrub_layer_vol) %>% mutate(plot_understory_vol = replace_na(plot_understory_vol, 0.0)), by = 'site')
  
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
  
  pnts = pnts %>% left_join(tree_div_metrics, by = 'site')
  
  # ggplot(pnts, aes(x = as.factor(`stratum`), y = plot_tree_all_diversity)) +
  #   geom_violin() +
  #   stat_summary(fun = mean, geom = "point") +
  #   ggtitle("Tree shannon diversity across strata") + theme_minimal()
  
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
  # metric = "plot_tree_all_density_" # "tree_all_density_", "tree_gte10cm_density_"
  # dps_long = pnts %>% select(site, `stage`, starts_with(metric)) %>%
  #   pivot_longer(cols = starts_with(metric), names_to = "species", values_to = "density") %>%
  #   mutate(species = gsub(metric, "", species))
  # ggplot(dps_long, aes(x = species, y = density, fill = species)) +
  #   geom_boxplot() +
  #   facet_wrap(~ stage, scales = "free_y") +
  #   coord_cartesian(ylim = c(0, 1250)) + 
  #   labs(title = "Tree species density by stage") +
  #   theme_minimal()
  
  # Distance to roads
  dist_road_paved = st_distance(pnts, paved_primary_roads)
  dist_road_paved = apply(dist_road_paved, 1, min)
  pnts$dist_road_paved = dist_road_paved
  mapview(pnts, zcol = "dist_road_paved") + mapview(paved_primary_roads)
  
  dist_roads_all = st_distance(pnts, roads)
  dist_roads_all = apply(dist_roads_all, 1, min)
  pnts$dist_roads_all = dist_roads_all
  mapview(pnts, zcol = "dist_roads_all") + mapview(roads)
  
  # Distance to water (type 1-3 and all types) [m]
  dist_watercourse_major = st_distance(pnts, boundary_watercourses)
  dist_watercourse_major = apply(dist_watercourse_major, 1, min)
  pnts$dist_watercourse_major = dist_watercourse_major
  mapview(pnts, zcol = "dist_watercourse_major") + mapview(boundary_watercourses)
  
  dist_watercourse_all = st_distance(pnts, watercourses)
  dist_watercourse_all = apply(dist_watercourse_all, 1, min)
  pnts$dist_watercourse_all = dist_watercourse_all
  mapview(pnts, zcol = "dist_watercourse_all") + mapview(watercourses)
  
  # # Distance to nearest edge [m]
  # # For each site, find its patch, then find the distance to the nearest non-patch cell
  # patch_ids = patches(rast_cover_clean, directions=8, values=TRUE)
  # pnts$dist_nearest_edge = NA
  # for (i in 1:nrow(pnts)) {
  #   d = pnts[i,]
  #   cover_class = terra::extract(rast_cover_clean, vect(d))[,2]
  #   pid = terra::extract(patch_ids, vect(d))[,2]
  #   m = trim(classify(patch_ids, cbind(pid, 1), others = NA))
  #   na_cells = which(is.na(values(m)))
  #   na_coords = xyFromCell(m, na_cells)
  #   d_coords = st_coordinates(d)
  #   distances = sqrt((na_coords[,1] - d_coords[1])^2 + (na_coords[,2] - d_coords[2])^2)
  #   min_dist = min(distances)
  #   pnts[i, 'dist_nearest_edge'] = min_dist
  #   # mapview(d) + mapview(m) + mapview(st_buffer(d, 100), col.regions = 'transparent', lwd = 2)
  # }
  # mapview(rast_cover_clean) + mapview(pnts, label = pnts$site, zcol = 'dist_nearest_edge') + mapview(st_buffer(pnts, 100), col.regions = 'transparent', lwd = 2)
  
  # TODO:
  message('Saving plot scale data cache ', path_data_plot_scale_out)
  dir.create(dirname(path_data_plot_scale_out), recursive = TRUE, showWarnings = FALSE)
  saveRDS(pnts, path_data_plot_scale_out)
  message("Finished generating plot scale data cache (", round(as.numeric(difftime(Sys.time(), time_start, units = 'mins')), 2), " minutes)")
  
} else { # overwrite_data_plot_scale_cache is FALSE
  message('Loading plot scale data from cache ', path_data_plot_scale_out)
  pnts = readRDS(path_data_plot_scale_out)
}

##############################################################################
# Focal patch scale and home range neighborhood scale covariates (limited to the extent of the species home range)

# TODO: DEBUG IN PROGRESS MANGO

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
  homeranges = data.frame(
    common_name = c("min", "median", "mean", "max"),
    home_range_radius_m = c(100, buff_median, buff_mean, buff_max)
  )
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
