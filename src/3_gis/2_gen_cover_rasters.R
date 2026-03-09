## 2_gen_cover_rasters.R #########################################################################################
# Generate and inspect raster(s) for landscape cover
#
## OUTPUT:
path_out = "data/cache/3_gis/2_gen_cover_rasters"
#
## INPUT:
# Base RS-FRIS data (0.1 acre resolution, i.e. ~404m2 or 20.10836 * 20.10836 m grain, roughly 1% of the area of a 100m radius circle)
# RS-FRIS 4.0 uses a combination of 2019 and 2020 photogrammetry.
# RS-FRIS 5.0 uses a combination of 2021 and 2022 photogrammetry. 
# TODO: Calculate on a yearly basis!
path_rsfris = "data/environment/rsfris_study_area"
path_trait_data = "data/cache/trait_data/trait_data.csv"
###########################################################################################################

source("src/global.R")

if (!dir.exists(path_out)) dir.create(path_out, recursive = TRUE)

species_trait_data = read_csv(path_trait_data, show_col_types = FALSE)

options(mapview.maxpixels = 2117676)

##############################################################################
# Study area, sites, boundaries, and helper functions

# ARU locations (sites)
# Used to define study area
mapview(aru_sites, label = aru_sites$name)

source("src/3_gis/1_preprocess_gis_data.R")

for (t in rsfris_version_years$year) {
  
  year_baseline = t
  version = rsfris_version_years %>% filter(year == t) %>% pull(version)

  # Store site-specific covariate data
  data_plot_scale = aru_sites
  
  # Derive land cover stages/strata ---------------------------------------------------------------
  
  path_rsfris_t = paste0(path_rsfris, "/", version)
  
  # WADNR delineated patch polygons. Manually determined from a combination of aerial imagery and remote sensing.
  wadnr_patch_sf = st_read('data/environment/GIS Data/Forest Development Strata/AgeStrataFixedDIS_RSFRIS20200130.shp') %>%
    st_transform(crs_m) %>% clean_names() %>% mutate(stratum = str_to_lower(stratum)) %>% rename(wadnr_patch_stratum = stratum)
  wadnr_patch_sf$wadnr_patch_stratum[wadnr_patch_sf$wadnr_patch_stratum == "initiation"] = stages_3[['class']][1]
  wadnr_patch_sf$wadnr_patch_stratum[wadnr_patch_sf$wadnr_patch_stratum == "canclose"]   = stages_3[['class']][1] # redefine canclose to initiation
  wadnr_patch_sf$wadnr_patch_stratum[wadnr_patch_sf$wadnr_patch_stratum == "stemex"]     = stages_3[['class']][2]
  wadnr_patch_sf$wadnr_patch_stratum[wadnr_patch_sf$wadnr_patch_stratum == "mature"]     = stages_3[['class']][3]
  data_plot_scale = st_join(data_plot_scale, wadnr_patch_sf)
  
  # Origin year is calculated from multiple data sources using the following logic: Origin year is reported at the pixel-level (1/10th ac scale). The default value is the predicted origin year, based on RS-FRIS models. These models rely on remotely sensed data (LiDAR and DAP) and are constructed primarily from height-to-age relationships. The default RS-FRIS predicted ages are overwritten using the following data, if available:
  # - Tree core data from DNR's historic inventory (FRIS) is used for stands whose origin year is 1900 or earlier.
  # - FRIS tree core data from younger stands (post-1900) or FRIS data based on stratified samples is used only if RS-FRIS data is not available.
  # - Information from the Land Resource Manager system (LRM) is used for completed harvests ('TEMP_RET_REM', 'TEMP_RET_1ST', 'VRH', 'SEEDTREE_INT', 'CLEAR_CUT', 'PATCH_REGEN', 'LANDUSE_CONV', 'SEEDTREE_REM', 'SHELTER_INT', 'SHELTER_REM'). STAND_ORIGIN_DT was used if populated, otherwise FMA_DT (<forest management activity date?>).
  # Data sources are mapped in the "Combined Origin Year Data Source" raster.
  rast_origin = round(load_raster(paste0(path_rsfris_t, '/ORIGIN_YEAR.tif')))
  summary((terra::values(rast_origin)))
  hist(na.omit(values(rast_origin)))
  
  # Flag patches of missing origin year data (these were harvested after baseline year)
  rast_origin_missing = rast_origin
  rast_origin_missing[rast_origin_missing < (year_baseline + 1)] = NA
  summary((terra::values(rast_origin_missing)))
  
  # Impute patches
  # (see Powell vegetation stages white paper; also consult ESRI World Imagery Wayback for visual inspection)
  rast_origin[rast_origin >= (year_baseline + 1)] = NA
  # mapview(wadnr_patch_sf) + mapview(rast_origin_missing)
  rast_wadnr_patch = rasterize(
    vect(wadnr_patch_sf), rast_origin_missing,
    field = "wadnr_patch_stratum"
  )
  
  rast_missing_masked = mask(rast_wadnr_patch, rast_origin_missing)
  rast_missing_imputed = catalyze(rast_missing_masked)
  rast_missing_imputed[rast_missing_imputed == 1] = year_baseline - stages_3 %>% filter(class == "compex") %>% pull(age_min)
  rast_missing_imputed[rast_missing_imputed == 2] = year_baseline - stages_3 %>% filter(class == "mature") %>% pull(age_min)
  rast_missing_imputed[rast_missing_imputed == 3] = year_baseline - stages_3 %>% filter(class == "standinit") %>% pull(age_min)
  
  rast_origin_imputed = cover(rast_missing_imputed, rast_origin)
  
  # Calculate age
  rast_age = round(year_baseline - rast_origin_imputed)
  
  # Determine thinning treatment
  # TODO: Check thinning class validity for sites: ca263i (mature), ap022i (standinit), az041i (standinit), 
  poly_thinning_treatment = st_read('data/environment/GIS Data/Forest Development Strata/ThinAfter94NoHarvSinceClipByInitBuf3.shp') %>% 
    st_transform(crs_m) %>% select(TECHNIQUE_, FMA_DT, FMA_STATUS) %>% janitor::clean_names() %>%
    mutate(technique = technique %>% str_to_lower(), fma_status = fma_status %>% str_to_lower()) %>%
    rename(thinning_treatment = technique, thinning_status = fma_status, thinning_date = fma_dt)
  stopifnot(sum(poly_thinning_treatment$thinning_status == "completed") == nrow(poly_thinning_treatment)) # check that all treatments are completed
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
  
  rast_thin = rasterize(vect(poly_thinning_treatment), rast_age, field = "thinning_status")
  rast_thin = as.numeric(rast_thin)
  
  ## Classify cover by age breaks
  
  rast_stage_3 = classify(rast_age, as.matrix(stages_3 %>% select(-class)), include.lowest = TRUE, right = FALSE)
  levels(rast_stage_3) = data.frame(ID = stages_3$idx, stage = stages_3$class)
  data_plot_scale$stage_3 = terra::extract(rast_stage_3, vect(data_plot_scale))[, 2]
  table(data_plot_scale$stage_3)
  
  rast_stage_4 = classify(rast_age, as.matrix(stages_4 %>% select(-class)), include.lowest = TRUE, right = FALSE)
  levels(rast_stage_4) = data.frame(ID = stages_4$idx, stage = stages_4$class)
  data_plot_scale$stage_4 = terra::extract(rast_stage_4, vect(data_plot_scale))[, 2]
  table(data_plot_scale$stage_4)
  
  rast_strata_4 = classify(rast_age, as.matrix(strata_4 %>% filter(class != "thin") %>% select(-class)), include.lowest = TRUE, right = FALSE)
  rast_strata_4 = cover(rast_thin, rast_strata_4)
  levels(rast_strata_4) = data.frame(ID = strata_4$idx, stage = strata_4$class)
  data_plot_scale$stratum_4 = terra::extract(rast_strata_4, vect(data_plot_scale))[, 2]
  table(data_plot_scale$stratum_4)
  
  rast_strata_5 = classify(rast_age, as.matrix(strata_5 %>% filter(class != "thin") %>% select(-class)), include.lowest = TRUE, right = FALSE)
  rast_strata_5 = cover(rast_thin, rast_strata_5)
  levels(rast_strata_5) = data.frame(ID = strata_5$idx, stage = strata_5$class)
  data_plot_scale$stratum_5 = terra::extract(rast_strata_5, vect(data_plot_scale))[, 2]
  table(data_plot_scale$stratum_5)
  
  # Assemble cover class raster ---------------------------------------------------------------
  
  # Loop over all strata/stages
  all_raster_list = list(
    "stage_3" = rast_stage_3,
    "stage_4" = rast_stage_4,
    "strata_4" = rast_strata_4,
    "strata_5" = rast_strata_5
  )
  rast_cover = list()
  for (rast_name in names(all_raster_list)) {
    message("Assembling cover class raster for: ", rast_name)
    r = all_raster_list[[rast_name]]
    print(unique(r))
  
    # Further delineate patch boundaries by roads, waterbodies, and watercourses
    template = rast(ext(r), resolution = res(r), crs = crs(r))
    
    message("Assembling water cover")
    
    watercourse_half_width = res(r)[1] / 2
    watercourses_buffered = (st_buffer(boundary_watercourses, dist = watercourse_half_width))
    rast_watercourses = as.factor(rasterize(
      vect(watercourses_buffered), template, field = "sl_wtrty_c"
    ))
    i = nrow(unique(r))+1
    rast_watercourses[rast_watercourses == 1] = i
    levels(rast_watercourses) = data.frame(ID = i, stage = "water")
    r = cover(rast_watercourses, r)
    
    rast_waterbodies = rasterize(
      vect(boundary_waterbodies), template, field = "objectid"
    )
    rast_waterbodies[!is.na(rast_waterbodies)] = i
    levels(rast_waterbodies) = data.frame(ID = i, stage = "water")
    r = cover(rast_waterbodies, r)
    
    # Fill missing water bodies off coast
    ocean_polyfill = st_polygon(list(rbind(
      c(-180, -90), c(-124.42, -90), c(-124.42, 47.63), c(-180, 47.63), c(-180, -90)
    ))) |> st_sfc(crs = 4326) |> st_transform(crs(rast_waterbodies))
    p = terra::project(vect(ocean_polyfill), crs(rast_waterbodies))
    mask_r = rasterize(p, rast_waterbodies, field=1, background=NA, touches=TRUE)
    rast_waterbodies[!is.na(mask_r[])] = i
    r = cover(rast_waterbodies, r)
    
    message("Assembling road cover")
    
    # "The Hoh-Clearwater Mainline is our only double-lane mainline, and it is about 26 feet wide for the asphalt surface. If we’re looking at right of way widths, we could easily assume a minimum width of about 50 feet for a 12-foot road, probably a good 60-80 feet for a 14-20 foot wide road, and about 100 feet for the Hoh Mainline. 100 feet might be good for US 101 as well for right of way width, and maybe about 30 feet for actual road surface width."
    paved_primary_roads = st_make_valid(roads %>% filter(road_usgs1 %in% c("Primary Highway", "Light-Duty Road")))
    road_half_width_primary = max(100 / 2 * conv_ft_to_m, res(r)[1] / 2)
    # "Our minimum road surface width is going to be 12 feet. That would cover most of our roads. Main arterials/single lane mainlines will have a minimum surface width of 14 feet, with a few up to 20 feet."
    paved_secondary_roads = st_make_valid(roads %>% filter(road_usgs1 %in% c("Light-Duty Road")))
    road_half_width_secondary = max(30 / 2 * conv_ft_to_m, res(r)[1] / 2)
    
    paved_roads_primary_buffered = (st_buffer(paved_primary_roads, dist = road_half_width_primary))
    paved_roads_secondary_buffered = (st_buffer(paved_secondary_roads, dist = road_half_width_secondary))
    
    rast_road_paved_primary = rasterize(
      vect(paved_roads_primary_buffered), template
    )
    i = nrow(unique(r))+1
    rast_road_paved_primary[rast_road_paved_primary == 1] = i
    levels(rast_road_paved_primary) = data.frame(ID = i, stage = "road_paved")
    r = cover(rast_road_paved_primary, r)
    
    rast_road_paved_secondary = rasterize(
      vect(paved_roads_secondary_buffered), template
    )
    rast_road_paved_secondary[rast_road_paved_secondary == 1] = i
    levels(rast_road_paved_secondary) = data.frame(ID = i, stage = "road_paved")
    r = cover(rast_road_paved_secondary, r)
    
    message("Cropping to study area")
    r = crop(r, ext(vect(study_area)))
    
    message("Imputing any remaining missing patches with cover class 'other'")
    rast_other = r
    rast_other[] = NA
    i = nrow(unique(r))+1
    rast_other[is.na(r[])] = i
    levels(rast_other) = data.frame(ID = i, stage = "other")
    r = cover(rast_other, r)
    
    rast_cover[[rast_name]] = r
    
    path_out_rast = paste0(path_out, "/rast_cover_", t, "_", rast_name, ".tif")
    writeRaster(r, path_out_rast, overwrite=TRUE)
    message(crayon::green("Cached raster cover data to:", path_out_rast))
  }
  
  # mapview(rast_cover[["strata_5"]]) +
  #   mapview(aru_sites, label = aru_sites$site) +
  #   mapview(st_buffer(aru_sites, 100), col.regions = 'transparent', lwd = 2)
  
  # Clean cover class raster ---------------------------------------------------------------
  
  min_species_homerange_area_m2 = round(pi * min(species_trait_data$home_range_radius_m)^2,0)
  rast_cover_clean = list()
  for (rast_name in names(all_raster_list)) {
  
    # Clean raster by generalizing minimum patch area and width
      
    message("Cleaning cover raster '", rast_name, "' (current time ", time_start <- Sys.time(), ")")
    
    r = rast_cover[[rast_name]]
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
    pb = progress_bar$new(format = "[:bar] :percent :elapsedfull (ETA :eta)", total = length(cls), clear = FALSE)
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
      patch_list[[as.character(cl)]] <- r_patches
      pb$tick()
    }
    
    # Merge patch ids for each raster type into one raster
    patch_ids = Reduce(cover, patch_list)
    
    # Create a raster stack (cover class and patch id)
    patch_stack = c(patch_ids, r)
    names(patch_stack) = c("patch_id", "cover_class")
    
    # mapview(patch_stack[["patch_id"]], col.regions = viridis) + 
    #   mapview(patch_stack[["cover_class"]])
    
    # Identify patches with maximum width less than minimum width ##################################################
    message("Identifying patch candidates with negligible width")
    narrow_patches = c()
    pids = unique(na.omit(as.vector(values(patch_stack[['patch_id']]))))
    pb = progress_bar$new(format = "[:bar] :percent :elapsedfull (ETA :eta)", total = length(pids), clear = FALSE)
    for (pid in pids) {
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
      pb$tick()
    }
    
    # Identify patch ids with negligible patch width
    patch_stack[["patch_narrow"]] = classify(patch_stack[["patch_id"]], rcl = cbind(narrow_patches, narrow_patches), others = NA)
    patch_narrow_ids = unique(patch_stack[["patch_narrow"]])[,1]
    
    # Reclassify those cells to value 0
    values(patch_stack[["cover_class"]])[which(values(!is.na(patch_stack[['patch_narrow']])))] = 0
    # mapview(patch_stack[['cover_class']])
    
    # Identify patches with area less than minimum area ##################################################
    message("Identifying patch candidates with negligible area")
    
    # Identify patches with negligible area
    freq_table = freq(patch_stack[["patch_id"]])
    small_patches = freq_table[freq_table$count < min_cells_area, "value"]
    patch_stack[["patch_small"]] = classify(patch_stack[["patch_id"]], rcl = cbind(small_patches, small_patches), others = NA)
    patch_small_ids = unique(patch_stack[["patch_small"]])[,1]
    
    # Reclassify those cells to value 0
    values(patch_stack[["cover_class"]])[which(values(!is.na(patch_stack[['patch_small']])))] = 0
    # mapview(patch_stack[['cover_class']])
    
    r_zero_alt = patch_stack[['cover_class']]
    values(r_zero_alt)[values(r_zero_alt) != 0] = NA
    # mapview(r_zero_alt)
    p = patches(r_zero_alt, directions = 8, values=TRUE)
    
    # Reclassify (generalize) negligible patches as the mode of adjacent nearest neighbors
    message("Reclassifying negligible patches as the mode of adjacent nearest neighbors")
    r_clean = patch_stack[["cover_class"]]
    i = 1
    patch_ids_to_generalize = unique(na.omit(values(p)))
    total = length(patch_ids_to_generalize)
    pb = progress_bar$new(format = "[:bar] :percent :elapsedfull (ETA :eta)", total = total, clear = FALSE)
    for (small_patch_id in patch_ids_to_generalize) {
      # print(round(i / total, 3))
      
      # mapview(trim(classify(patch_stack[["patch_id"]], cbind(small_patch_id, 1), others = NA)))
      
      patch_cells = which(values(p) == small_patch_id)
      
      # Find adjacent cells to the patch
      adj_cells = adjacent(r, cells = patch_cells, directions = 4, pairs = TRUE)
      
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
      values(r_clean)[patch_cells] <- majority_class
      
      i = i + 1
      pb$tick()
    }
    r_clean[is.nan(r_clean)] = NA
    # mapview(r_clean) + mapview(patch_stack[['cover_class']])
    
    r_clean_fact = r_clean
    names(r_clean_fact) = names(r)
    levels(r_clean_fact) = levels(r)
    
    rast_cover_clean[[rast_name]] = r_clean_fact
    
    # Cache clean cover raster
    path_out_rast_clean = paste0(path_out, "/rast_cover_", t, "_clean_", rast_name, ".tif")
    writeRaster(r_clean_fact, path_out_rast_clean, overwrite=TRUE)   
    message(crayon::green("Cached raster cover clean data to:", path_out_rast_clean))
  }
  
  # Overwrite site management strata data with manual logs
  site_key = read_csv(path_site_key, show_col_types = FALSE) %>%
    mutate(site = str_to_lower(site)) %>% mutate(site_agg = str_to_lower(site_agg)) %>%
    select(site, stratum) %>% filter(stratum != "HARVESTED") %>% distinct() %>%
    mutate(stratum = dplyr::recode(stratum,
                            "STAND INIT" = "standinit",
                            "COMP EXCL"  = "compex",
                            "THINNED"    = "thin",
                            "MATURE"     = "mature"))
  
  site_cover_class_sf = data_plot_scale %>%
    left_join(site_key, by = "site") %>%  mutate(stratum_4 = coalesce(stratum, stratum_4)) %>% select(-stratum)
  site_cover_class_sf = site_cover_class_sf %>%
    left_join(site_key, by = "site") %>%
    mutate(
      stratum_5 = if_else(
        stratum_5 %in% c("standinit", "thin", "compex") & !is.na(stratum),
        stratum,
        stratum_5
      )
    ) %>%
    select(-stratum)
  
  # Cache site cover class data
  path_out_site_cover_class_sf = paste0(path_out, "/site_cover_class_", t, "_sf.rds")
  write_rds(site_cover_class_sf, path_out_site_cover_class_sf)
  message(crayon::green("Cached site cover class sf data to:", path_out_site_cover_class_sf))
  
  message(crayon::green("Finished generating cover cache (", round(as.numeric(difftime(Sys.time(), time_start, units = 'mins')), 2), " minutes)"))
  
  # Compare original vs clean
  plot(rast_cover[["stage_3"]], main = "stage_3 (original)")
  plot(rast_cover_clean[["stage_3"]], main = "stage_3 (clean)")
  
  # View each raster
  plot(rast_cover[["stage_3"]], main = "Developmental stage (3-class)", col = c(
    "standinit"  = "#d8c18a",
    "compex"     = "#3c8273",
    "mature"     = "#9b652b",
    "water"      = "#6495ed",
    "road_paved" = "gray50",
    "other"      = "gray80"
  ))
  plot(rast_cover[["stage_4"]], main = "Developmental stage (4-class)", col = c(
    "standinit"  = "#d8c18a",
    "compex"     = "#3c8273",
    "underdev"   = "#9b652b",
    "old"        = "#5C4033",
    "water"      = "#6495ed",
    "road_paved" = "gray50",
    "other"      = "gray80"
  ))
  plot(rast_cover[["strata_4"]], main = "Management strata (4-class)", col = c(
    "thin"       = "#b2675e",
    "standinit"  = "#d8c18a",
    "compex"     = "#3c8273",
    "mature"     = "#9b652b",
    "water"      = "#6495ed",
    "road_paved" = "gray50",
    "other"      = "gray80"
  ))
  plot(rast_cover[["strata_5"]], main = "Managment strata (5-class)", col = c(
    "thin"       = "#b2675e",
    "standinit"  = "#d8c18a",
    "compex"     = "#3c8273",
    "underdev"   = "#9b652b",
    "old"        = "#5C4033",
    "water"      = "#6495ed",
    "road_paved" = "gray50",
    "other"      = "gray80"
  ))
  
  # Inspect dynamically
  mapview(rast_cover_clean[["strata_5"]],
          alpha.regions = 1.0,
          col.regions = c(
            "thin"       = "#b2675e",
            "standinit"  = "#d8c18a",
            "compex"     = "#3c8273",
            "underdev"   = "#9b652b",
            "old"        = "#5C4033",
            "water"      = "#6495ed",
            "road_paved" = "gray50",
            "other"      = "gray80"
          )) +
    mapview(wadnr_patch_sf) +
    mapview(poly_thinning_treatment) +
    mapview(data_plot_scale, label = data_plot_scale$site) +
    mapview(st_buffer(data_plot_scale, 100), col.regions = 'transparent', lwd = 2)

}
