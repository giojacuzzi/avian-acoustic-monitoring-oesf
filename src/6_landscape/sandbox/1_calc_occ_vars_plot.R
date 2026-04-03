source("src/3_gis/3_calc_occ_vars_config.R")

source("src/3_gis/1_preprocess_gis_data.R")
paved_roads = st_make_valid(roads %>% filter(road_usgs1 %in% c("Primary Highway", "Light-Duty Road")))

##############################################################################
# Local plot scale covariates

if (overwrite_data_plot_scale_cache) {
  
  message("Generating plot scale data cache (current time ", time_start <- Sys.time(), ")")
  
  # Plot-level spatial scale buffer
  plot_buffer = 100 # 100 meters
  
  # Elevation [m]
  rast_elevation = load_raster("data/environment/elevation/elevation.tif")
  pnts$elevation = extract(rast_elevation, vect(pnts))[,2]
  
  # Age (mean and cv) [#]
  # TODO: do this in script #2
  pnts$age_point = terra::extract(rast_age, vect(pnts))[, 2]
  stats_age      = compute_raster_buffer_value_func(rast_age, pnts, plot_buffer, summary_stats)
  pnts$age_mean = as.numeric(stats_age[,2])
  pnts$age_cv   = as.numeric(stats_age[,4])
  hist(pnts$age_point, breaks = seq(0, max(pnts$age_point, na.rm = TRUE) + 10, by = 10))
  hist(pnts$age_mean, breaks = seq(0, max(pnts$age_mean, na.rm = TRUE) + 10, by = 10))
  table(pnts$age_point)
  table(round(pnts$age_mean))
  
  # Distance to edge
  rast_edges = boundaries(rast_cover, inner = TRUE, classes = TRUE)
  rast_edges[rast_edges == 0] = NA
  # Compute a distance raster (m)
  rast_dist_edge = distance(rast_edges)
  # Get distance at each site
  pnts = pnts %>% mutate(dist_edge = terra::extract(rast_dist_edge, vect(pnts))[, 2])
  
  # Distance to roads
  dist_road_paved = st_distance(pnts, paved_roads)
  dist_road_paved = apply(dist_road_paved, 1, min)
  pnts$dist_road_paved = dist_road_paved
  # mapview(pnts, zcol = "dist_road_paved") + mapview(paved_roads)
  
  dist_road_all = st_distance(pnts, roads)
  dist_road_all = apply(dist_road_all, 1, min)
  pnts$dist_road_all = dist_road_all
  # mapview(pnts, zcol = "dist_road_all") + mapview(roads, zcol = "road_usgs1")
  
  # Distance to water (type 1-3 and all types) [m]
  dist_watercourse_major = st_distance(pnts, boundary_watercourses)
  dist_watercourse_major = apply(dist_watercourse_major, 1, min)
  pnts$dist_watercourse_major = dist_watercourse_major
  # mapview(pnts, zcol = "dist_watercourse_major") + mapview(boundary_watercourses)
  
  dist_watercourse_all = st_distance(pnts, watercourses)
  dist_watercourse_all = apply(dist_watercourse_all, 1, min)
  pnts$dist_watercourse_all = dist_watercourse_all
  # mapview(pnts, zcol = "dist_watercourse_all") + mapview(watercourses)
  
  # Distance to nearest edge [m]
  # For each site, find its patch, then find the distance to the nearest non-patch cell
  # hr_patch_ids = patches(rast_cover, directions=8, values=TRUE)
  # pnts$dist_nearest_edge = NA
  # for (i in 1:nrow(pnts)) {
  #   d = pnts[i,]
  #   cover_class = terra::extract(rast_cover, vect(d))[,2]
  #   pid = terra::extract(hr_patch_ids, vect(d))[,2]
  #   m = trim(classify(hr_patch_ids, cbind(pid, 1), others = NA))
  #   na_cells = which(is.na(values(m)))
  #   na_coords = xyFromCell(m, na_cells)
  #   d_coords = st_coordinates(d)
  #   distances = sqrt((na_coords[,1] - d_coords[1])^2 + (na_coords[,2] - d_coords[2])^2)
  #   min_dist = min(distances)
  #   pnts[i, 'dist_nearest_edge'] = min_dist
  # mapview(d) + mapview(m) + mapview(st_buffer(d, 100), col.regions = 'transparent', lwd = 2)
  # }
  # mapview(rast_cover) + mapview(pnts, label = pnts$site, zcol = 'dist_nearest_edge') + mapview(st_buffer(pnts, 100), col.regions = 'transparent', lwd = 2)
  
  if (pnts_name == "sites") {
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
    intersect(pnts$site, data_plot$site)
    length(intersect(pnts$site, data_plot$site))
    setdiff(data_plot$site, pnts$site) # This should be 0
    pnts$hs = FALSE
    pnts[pnts$site %in% data_plot$site, 'hs'] = TRUE
    mapview(pnts, zcol = 'hs')
    
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
    
    ggplot(pnts, aes(x = stratum_4, y = plot_tree_all_diversity)) +
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
    dps_long = pnts %>% select(site, stratum_4, starts_with(metric)) %>%
      pivot_longer(cols = starts_with(metric), names_to = "species", values_to = "density") %>%
      mutate(species = gsub(metric, "", species))
    ggplot(dps_long, aes(x = species, y = density, fill = species)) +
      geom_boxplot() +
      facet_wrap(~ stratum_4, scales = "free_y") +
      coord_cartesian(ylim = c(0, 1250)) +
      labs(title = "Tree species density by stage") +
      theme_minimal()
  }
  
  dir.create(dirname(path_data_plot_scale_out), recursive = TRUE, showWarnings = FALSE)
  saveRDS(pnts, path_data_plot_scale_out)
  message(crayon::green("Cached plot scale data to", path_data_plot_scale_out, "(", round(as.numeric(difftime(Sys.time(), time_start, units = 'mins')), 2), "min )"))
  
} else { # overwrite_data_plot_scale_cache is FALSE
  message('Loading plot scale data from cache ', path_data_plot_scale_out)
  pnts = readRDS(path_data_plot_scale_out)
}