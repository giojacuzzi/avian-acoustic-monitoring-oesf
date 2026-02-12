# 1_assemble_detection_arrays.R ####################################################################################
# Derive MSOM-ready data arrays:
# - Survey date (yday) array (site x survey x season)
# - Putative detection array (site x survey x species x season) from both manually validated and thresholded prediction data.
#
# CONFIG:
min_sites_detected = 1 # Minimum number of sites present to retain a species for analysis
#
# OUTPUT:
out_cache_dir  = "data/cache/4_msom/1_assemble_detection_arrays"
path_out_y     = paste0(out_cache_dir, "/y.rds")
path_out_xyday = paste0(out_cache_dir, "/xyday.rds")
#
# INPUT:
path_community_array_predictions = "data/cache/1_pam/3_derive_observation_data/community_array_predictions.rds" # TODO: rename cache directory and re-run script 3
path_community_array_surveydates = "data/cache/1_pam/3_derive_observation_data/community_array_surveydates.rds"
path_calibration_results         = "data/cache/1_pam/1_classifier_calibration/calibration_results_raw.csv"
path_annotations                 = "data/cache/1_pam/1_classifier_calibration/annotations_clean.csv"
##################################################################################################################

source("src/global.R")

if (!dir.exists(out_cache_dir)) dir.create(out_cache_dir, recursive = TRUE)

# Load dependencies ----------------------------------------------------------------------------------------------

message("Loading species-specific thresholds from ", path_calibration_results)
calibration_results = read_csv(path_calibration_results, show_col_types = FALSE)

message("Loading site confirmation annotations from ", path_annotations)
annotations = read_csv(path_annotations, show_col_types = FALSE)

message("Loading community survey date array from ", path_community_array_surveydates)
community_array_surveydates = readRDS(path_community_array_surveydates)

message("Loading community prediction array from ", path_community_array_predictions)
community_array_predictions = readRDS(path_community_array_predictions)

sites   = dimnames(community_array_predictions)[["site"]]
surveys = dimnames(community_array_predictions)[["survey"]]
seasons = dimnames(community_array_predictions)[["season"]]
message("Found ", length(sites), " sites, ", length(seasons), " seasons, and up to ", length(surveys), " surveys per site")

# Derive survey date (yday) arrays ------------------------------------------------------------------------

message("Deriving survey date (yday) arrays x_yday")
x_yday = setNames(
  lapply(seq_along(seasons), function(t) {
    matrix(
      unlist(lapply(community_array_surveydates[, , t, 1], function(x) {
        if (!is.null(x)) yday(x$survey_date) else NA
      })),
      nrow = dim(community_array_surveydates)[1],
      ncol = dim(community_array_surveydates)[2],
      dimnames = dimnames(community_array_surveydates)[1:2]
    )
  }),
  seasons
)
# Inspect examples
x_yday[["2020"]]

# Derive putative (uncertain) detection-nondetection arrays for each species with optimized thresholds -----------------------

message("Deriving putative observation arrays for each species with optimized thresholds (detection '1' nondetection '0')")
species = dimnames(community_array_predictions)[["common_name"]]
species = species[!str_starts(species, "abiotic")] # Manually exclude non-avian classes
species = species[!str_starts(species, "biotic")]
ylist   = setNames(lapply(species, function(x) {
  setNames(vector("list", length(seasons)), seasons)
}), species)

class_discrepancies = sort(c(setdiff(species, calibration_results$common_name), setdiff(calibration_results$common_name, species)))
if (length(class_discrepancies) > 0) {
  message(crayon::yellow("WARNING:", length(class_discrepancies), "class discrepancies"))
  message(crayon::yellow(paste(class_discrepancies, collapse = ", ")))
}

species_thresholds = calibration_results %>% rename(species = common_name) %>% select(species, model, threshold)

# Populate putative detection-nondetection matrices via species-specific score thresholds
for (t in seasons) {
  message("Season ", t)
  pb = progress_bar$new(format = "[:bar] :percent :elapsedfull (ETA :eta)", total = length(species), clear = FALSE)
  for (i in species) {
    n_row        = dim(community_array_predictions)[1]
    n_col        = dim(community_array_predictions)[2]
    dim_names    = dimnames(community_array_predictions)[1:2]
    species_data = community_array_predictions[, , t, i]
    
    if (i %in% species_thresholds$species) {
      sp_threshdata = species_thresholds %>% filter(species == i)
      model     = sp_threshdata %>% pull(model)
      threshold = sp_threshdata %>% pull(threshold)
      
      confidence_model = switch(model, source = "confidence_source", target = "confidence_target", stop("Incompatible model type ", model))
      
      mat_obs = matrix(
        unlist(lapply(species_data, function(x) {
          if (!is.null(x)) {
            as.integer(any(x[[confidence_model]] > threshold, na.rm = TRUE))
          } else {
            NA
          }
        })),
        nrow = n_row, ncol = n_col, dimnames = dim_names
      )
    } else {
      message(crayon::yellow("WARNING: No threshold found for", i))
      # There is no threshold for this species
      mat_obs = matrix(
        unlist(lapply(species_data, function(x) if (!is.null(x)) 0.0 else NA)),
        nrow = n_row, ncol = n_col, dimnames = dim_names)
    }
    ylist[[i]][[t]] = mat_obs # Store the resulting putative observations
    pb$tick()
  }
}

# Overwrite manually reviewed site confirmations ---------------------------------------------------------

# Incorporate site confirmations, including manual observations as value 2
message("Overwriting manually reviewed site confirmations with value '2'")

for (i in species) {
  species_confirmations = annotations %>%
    select(all_of(i), serialno, season, yday, site, site_agg, file) %>%
    rename(confirmation = 1) %>% filter(confirmation == "1") %>%
    distinct()
  
  for (r in seq_len(nrow(species_confirmations))) {
    
    file = species_confirmations$file[r]
    j    = species_confirmations$site[r]
    t    = species_confirmations$season[r] %>% as.character()
    yday = species_confirmations$yday[r]
    
    if (is.na(j)) {
      message(yellow("WARNING: No matching site for", i, "confirmation in file", file))
      next
    }
    
    j_ydays = x_yday[[t]][j, ]
    k = which(!is.na(j_ydays) & j_ydays == yday)
    if (length(k) != 0) {
      ylist[[i]][[t]][j,k] <- 2 # Overwrite with confirmed detection
    }
  }
}

# Inspect examples
ylist[["american robin"]][["2020"]]
any(ylist[["bewick's wren"]][["2020"]] == 1, na.rm = TRUE)

# Exclude species that were not detected -----------------------------------------------------------

# Determine naive species site occurence per year and in total
naive_occurrence_per_year = lapply(names(ylist[[1]]), function(season) {
  sapply(ylist, function(species_list) {
    mat <- species_list[[season]]
    sum(apply(mat, 1, function(x) any(x >= 1, na.rm = TRUE)))
  })
})
total_occurrences = Reduce(`+`, naive_occurrence_per_year)
total_occurrences = tibble(
  species = names(total_occurrences),
  total_occurrences = as.numeric(total_occurrences)
)
print(total_occurrences, n = Inf)

message("Excluding species with insufficient detections (fewer than ", min_sites_detected, " sites):")
species_to_remove = total_occurrences %>% filter(total_occurrences < min_sites_detected) %>% pull(species)
print(species_to_remove)

ylist[species_to_remove] = NULL
species = names(ylist)

# Observed naive species occurrence
total_occurrences = total_occurrences %>% filter(species %in% names(ylist))
message(length(species), " species detected")
message("Most commonly observed species:")
print(total_occurrences %>% arrange(desc(total_occurrences)), n = 15)
message("Least commonly observed species:")
print(total_occurrences %>% arrange(total_occurrences), n = 15)

# Exclude any sites that were not surveyed -----------------------------------------------------------

surveys_per_season = lapply(names(ylist[[1]]), function(season) {
  season_counts = sapply(ylist, function(species_list) {
    rowSums(!is.na(species_list[[season]]))
  })
  as.data.frame(season_counts)
})
site_survey_counts = lapply(surveys_per_season, function(df) {
  rowSums(df, na.rm = TRUE)
})
site_survey_counts = rowSums(do.call(cbind, site_survey_counts)) # combine across years per site
sites_not_surveyed = unlist(unique(lapply(site_survey_counts, function(site_counts) {
  names(site_counts)[site_counts == 0]
})))
if (length(sites_not_surveyed) > 0) {
  message("Excluding ", length(sites_not_surveyed), " sites with no surveys")
  ylist = lapply(ylist, function(species_mat_list) {
    lapply(species_mat_list, function(mat) { mat[!(rownames(mat) %in% sites_not_surveyed), , drop = FALSE] })
  })
  x_yday = lapply(x_yday, function(mat) { mat[!(rownames(mat) %in% sites_not_surveyed), , drop = FALSE] })
}

# Convert to an MSOM-ready numeric array --------------------------------------------------------------------

message("Formatting observation data for MSOM")

species = names(ylist)
seasons = names(ylist[[1]])
sites   = dimnames(ylist[[1]][[1]])$site
surveys = dimnames(ylist[[1]][[1]])$survey

# Format ylist as 4D array (site × survey × season x species)
y = array(NA, dim = c(length(sites), length(surveys), length(seasons), length(species)),
          dimnames = list(site = sites, survey = surveys, season = seasons, species = species))
for (t in seasons) {
  for (i in species) {
    y[ , , t, i] = as.matrix(ylist[[i]][[t]])
  }
}
# access e.g. y[ , , "2020", "common raven"]

# Format x_yday as 3D array (site x survey x season)
x_yday = simplify2array(x_yday)
names(dimnames(x_yday)) = c("site", "survey", "season")

# Left-align data (moving any missing NA surveys to the right) to allow for direct indexing by number of surveys per site
left_align_row = function(x) {
  non_na = x[!is.na(x)]
  c(non_na, rep(NA, length(x) - length(non_na)))
}

# Align y and yday data
y_unaligned = y
x_yday_unaligned = x_yday
for (t in dimnames(y)[['season']]) {
  for (i in dimnames(y)[['species']]) {
    sp_season_mat               = y[, , t, i] # Extract site × survey matrix for this season × species
    sp_season_aligned           = t(apply(sp_season_mat, 1, left_align_row)) # Left-align across surveys for each site
    dimnames(sp_season_aligned) = dimnames(sp_season_mat) # Restore dimnames
    y[, , t, i]                 = sp_season_aligned # Put back into aligned array
  }
  x_season_mat               = x_yday[, , t] # Extract site × survey matrix for this season
  x_season_aligned           = t(apply(x_season_mat, 1, left_align_row)) # Left-align across surveys for each site
  dimnames(x_season_aligned) = dimnames(x_season_mat) # Restore dimnames
  x_yday[, , t]              = x_season_aligned # Put back into aligned array
  
}

# Cache results -------------------------------------------------------------------------------------

saveRDS(y, path_out_y)
message(crayon::green("Cached y to", path_out_y))

saveRDS(x_yday, path_out_xyday)
message(crayon::green("Cached xyday to", path_out_xyday))
