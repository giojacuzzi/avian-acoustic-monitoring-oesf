############################################################################################################################
# A standard single-season multi-species occupancy model assuming no false positives
#
model_name = "multiseason"
generate_diagnostic_plots = FALSE
############################################################################################################################

path_community_survey_data = "data/cache/1_derive_community_survey_data/community_survey_data_2025-08-15.rds"
path_species_thresholds    = "data/cache/1_calculate_species_thresholds/species_thresholds.csv"
path_plot_scale_data       = "data/cache/occurrence_covariates/data_plot_scale.rds"
path_homerange_scale_data  = "data/cache/occurrence_covariates/data_homerange_scale.rds"
path_out = paste0("data/cache/models/", model_name, "_", format(Sys.Date(), "%Y-%m-%d"), ".rds")

pacman::p_load(progress, car, tidyverse, jagsUI, MCMCvis, glue, ggplot2, ggrepel)
theme_set(theme_classic())

############################################################################################################################
# Load occurrence covariate data

# local plot scale
occ_data_plot_shared = readRDS(path_plot_scale_data) %>% sf::st_drop_geometry() %>% arrange(site) %>% mutate(site = tolower(site)) %>% select(
  site, elev
)
occ_data_plot_rs = readRDS(path_homerange_scale_data)[['plot']] %>% arrange(site) %>% mutate(site = tolower(site)) %>% select(
  site, homerange_downvol_mean, homerange_htmax_cv, homerange_qmd_mean, homerange_treeden_all_mean, homerange_treeden_gt4in_dbh_mean, homerange_ba_mean, homerange_snagden_gt15dbh_mean
)
occ_data_plot = occ_data_plot_shared %>% full_join(occ_data_plot_rs, by = "site") %>% select(
  site, elev, homerange_treeden_all_mean, homerange_qmd_mean, homerange_htmax_cv
)

# species-specific homerange scales
occ_data_homerange = readRDS(path_homerange_scale_data)
occ_data_homerange = map(occ_data_homerange, ~
                           .x %>%
                           arrange(site) %>%
                           mutate(site = tolower(site)) %>%
                           select(
                             buffer_radius_m,
                             site,
                             cover_forest_diversity,
                             density_edge_cw,
                             density_roads_paved,
                             density_streams_major,
                             focalpatch_area_homeange_pcnt,
                             prop_abund_standinit,
                             prop_abund_lsog,
                             prop_abund_comthin,
                             shape_idx
                           )
)
names(occ_data_homerange) = tolower(names(occ_data_homerange))

############################################################################################################################
# Load species-specific thresholds
message("Loading species-specific thresholds")
species_thresholds = read.csv(path_species_thresholds)
species_thresholds_source = species_thresholds %>% filter(model == 'source')
species_thresholds_target = species_thresholds %>% filter(model == 'target')
species_only_in_source = species_thresholds_source %>% filter(!species %in% species_thresholds_target$species)

# Manually choose bioacoustic classifier model for specific species
target_species = species_thresholds_target %>% filter(species %in% c(
  "marbled murrelet",
  "western screech-owl",
  "sooty grouse",
  "northern pygmy-owl",
  "western wood-pewee",
  "red-breasted nuthatch",
  "northern saw-whet owl",
  "white-crowned sparrow",
  "townsend's warbler",
  "dark-eyed junco",
  "hermit thrush",
  "golden-crowned kinglet",
  "song sparrow",
  "band-tailed pigeon",
  "pileated woodpecker",
  "rufous hummingbird",
  "red crossbill"
))
source_species = species_thresholds_source %>% filter(model == 'source') %>% filter(species %in% c(
  "ruby-crowned kinglet",
  "violet-green swallow",
  "american robin",
  "wilson's warbler",
  "spotted towhee",
  "purple finch",
  "olive-sided flycatcher",
  "western tanager",
  "hutton's vireo",
  "black-throated gray warbler",
  "varied thrush",
  "pacific-slope flycatcher",
  "pacific wren",
  "swainson's thrush",
  "barred owl",
  "belted kingfisher",
  "hairy woodpecker",
  "northern flicker",
  "hammond's flycatcher",
  "common raven"
))
species_thresholds_manual_selection = rbind(species_only_in_source, source_species)
species_thresholds_manual_selection = rbind(species_thresholds_manual_selection, target_species)
stopifnot(nrow(species_thresholds_manual_selection) == length(species_thresholds_source$species))

# Set minimum threshold to 0.5
species_thresholds_manual_selection = species_thresholds_manual_selection %>%
  mutate(
    threshold = ifelse(n_pos == 0, NA, ifelse(t_conf_tp >= 0.5, t_conf_tp, 0.5)),
    precision = ifelse(t_conf_tp >= 0.5, precision_tp, precision_0.5),
    recall    = ifelse(t_conf_tp >= 0.5, recall_tp, recall_0.5)
  ) %>% select(species, model, threshold, precision, recall, auc_pr, auc_roc, n_pos, n_neg) %>% arrange(species)
species_thresholds = species_thresholds_manual_selection

# TODO: Finalize species-specific thresholds to prevent need for manual override here
species_thresholds = species_thresholds %>%
  mutate(threshold = if_else(species == "vaux's swift", 0.5, threshold))

############################################################################################################################
# Derive putative observation and survey date matricies
message("Loading community survey data")
community_survey_data = readRDS(path_community_survey_data)
dimnames(community_survey_data)

message("Deriving detection-nondetection and yday matricies for each species")
seasons = dimnames(community_survey_data)[["season"]]
species = dimnames(community_survey_data)[["common_name"]]
ylist   = setNames(lapply(species, function(x) {
    setNames(vector("list", length(seasons)), seasons)
  }), species)
xlist_yday = setNames(lapply(species, function(x) {
  setNames(vector("list", length(seasons)), seasons)
}), species)

species_discrepancies = sort(c(setdiff(species, species_thresholds$species), setdiff(species_thresholds$species, species)))
if (length(species_discrepancies) > 0) {
  message(crayon::yellow("WARNING:", length(species_discrepancies), "species discrepancies"))
  message(crayon::yellow(paste(species_discrepancies, collapse = ", ")))
}

# Populate putative observation matrices (detection-nondetection via thresholded confidence scores)
for (t in seasons) {
  message("Populating putative observation matrices for season: ", t)
  pb = progress_bar$new(format = "[:bar] :percent :elapsedfull (ETA :eta)", total = length(species), clear = FALSE)
  for (i in species) {
    n_row = dim(community_survey_data)[1]
    n_col = dim(community_survey_data)[2]
    dim_names = dimnames(community_survey_data)[1:2]
    species_data = community_survey_data[, , t, i]
    
    # Use probablistic thresholding
    if (i %in% species_thresholds$species) {
      sp_threshdata = species_thresholds %>% filter(species == i)
      model     = sp_threshdata %>% pull(model)
      threshold = sp_threshdata %>% pull(threshold)
      if (model == "source") {
        mat_obs = matrix(
          unlist(lapply(species_data, function(x) if (!is.null(x)) as.integer(any(x$confidence_source >= threshold, na.rm = TRUE)) else NA)),
          nrow = n_row, ncol = n_col, dimnames = dim_names)
      } else if (model == "target") {
        mat_obs = matrix(
          unlist(lapply(species_data, function(x) if (!is.null(x)) as.integer(any(x$confidence_target >= threshold, na.rm = TRUE)) else NA)),
          nrow = n_row, ncol = n_col, dimnames = dim_names)
      }
    } else {
      # There is no threshold for this species
      mat_obs = matrix(
        unlist(lapply(species_data, function(x) if (!is.null(x)) 0.0 else NA)),
        nrow = n_row, ncol = n_col, dimnames = dim_names)
    }
    ylist[[i]][[t]] = mat_obs # Store the resulting putative observations
    pb$tick()
  }
}

# TODO: Incorporate site confirmations

# Survey date matrix (day of year)
x_yday = setNames(
  lapply(seq_along(seasons), function(t) {
    matrix(
      unlist(lapply(community_survey_data[, , t, 1], function(x) {
        if (!is.null(x)) yday(x$survey_date) else NA
      })),
      nrow = dim(community_survey_data)[1],
      ncol = dim(community_survey_data)[2],
      dimnames = dimnames(community_survey_data)[1:2]
    )
  }),
  seasons
)

# Discard sites with no environmental data
sites_missing_environmental_data = setdiff(dimnames(community_survey_data)$site, occ_data_plot$site)
if (length(sites_missing_environmental_data) > 0) {
  message("Discarding ", length(sites_missing_environmental_data), " sites with missing environmental data")
  ylist = lapply(ylist, function(species_mat_list) {
    lapply(species_mat_list, function(mat) { mat[!(rownames(mat) %in% sites_missing_environmental_data), , drop = FALSE] })
  })
  x_yday = lapply(x_yday, function(mat) { mat[!(rownames(mat) %in% sites_missing_environmental_data), , drop = FALSE] })
}
# Discard sites with no observations
sites_with_environmental_data_missing_observations = setdiff(occ_data_plot$site, dimnames(community_survey_data)$site)
if (length(sites_with_environmental_data_missing_observations) > 0) {
  message("Discarding ", length(sites_with_environmental_data_missing_observations), " sites with missing observations")
  ylist = lapply(ylist, function(species_mat_list) {
    lapply(species_mat_list, function(mat) { mat[!(rownames(mat) %in% sites_with_environmental_data_missing_observations), , drop = FALSE] })
  })
  x_yday = lapply(x_yday, function(mat) { mat[!(rownames(mat) %in% sites_with_environmental_data_missing_observations), , drop = FALSE] })
}

# TODO: Discard sites with no survey observations and surveys with no site observations?
# surveys_per_season = lapply(names(ylist[[1]]), function(season) {
#   season_counts = sapply(ylist, function(species_list) {
#     rowSums(!is.na(species_list[[season]]))
#   })
#   as.data.frame(season_counts)
# })
# site_survey_counts = lapply(surveys_per_season, function(df) {
#   rowSums(df, na.rm = TRUE)
# })
# sites_not_surveyed = lapply(site_survey_counts, function(site_counts) {
#   names(site_counts)[site_counts == 0]
# })
# 
# survey_site_counts = lapply(surveys_per_season, function(df) {
#   colSums(df, na.rm = TRUE)
# })
# sites_not_surveyed = lapply(site_survey_counts, function(site_counts) {
#   names(site_counts)[site_counts == 0]
# })
# 
# 
# survey_site_counts = lapply(ylist, function(x) { colSums(!is.na(x))})
# sites_per_survey = as.data.frame(t(do.call(rbind, survey_site_counts)))
# surveys_not_conducted = rownames(sites_per_survey)[rowSums(sites_per_survey) == 0]
# if (length(surveys_not_conducted) > 0) {
#   message("Discarding ", length(surveys_not_conducted), " surveys with no site observations")
#   ylist = lapply(ylist, function(mat) { mat[, !(colnames(mat) %in% surveys_not_conducted), drop = FALSE] })
#   x_yday = x_yday[, !colnames(x_yday) %in% surveys_not_conducted]
# }

sites   = dimnames(ylist[[1]][[1]])$site
surveys = dimnames(ylist[[1]][[1]])$survey

# Inspect the detection history and covariate data
message("Total number of sites: ", length(sites))
# lapply(ylist, head)
# head(x_yday)
sites_per_season = sapply(names(ylist[[1]]), function(t) {
  mat <- ylist[[1]][[t]]
  sum(rowSums(!is.na(mat)) > 0)
})
print(sites_per_season)

for (t in seasons) {
  n_surveys_per_site = apply(!is.na(ylist[[1]][[t]]), 1, sum)
  message("Season ", t, ": ", sum(n_surveys_per_site), " total sampling periods (surveys) conducted across ", sites_per_season[t], " sampling units (sites)")
  message("Sampling periods (surveys) conducted per site: median ", median(n_surveys_per_site), ", range ", min(n_surveys_per_site), "–", max(n_surveys_per_site))
  print(table(n_surveys_per_site))
}

############################################################################################################################
# Inspect naive detections, occurrence, and richness

# Obtain naive occurrence stats and exclude species detected at fewer than a minimum number of sites
# naive_occurrence = sapply(ylist, function(mat) { sum(apply(mat, 1, function(x) any(x == 1, na.rm = TRUE))) })
# naive_occurrence = data.frame(species = names(naive_occurrence), nsites = naive_occurrence) %>%
#   arrange(desc(nsites)) %>% mutate(species = factor(species, levels = rev(species))) %>% mutate(prob = nsites / length(sites))

naive_occurrence_per_year = lapply(names(ylist[[1]]), function(season) {
  sapply(ylist, function(species_list) {
    mat <- species_list[[season]]
    sum(apply(mat, 1, function(x) any(x == 1, na.rm = TRUE)))
  })
})
total_occurrences = Reduce(`+`, naive_occurrence_per_year)
naive_occurrence_per_year_df = do.call(rbind, lapply(seq_along(naive_occurrence_per_year), function(i) {
  data.frame(
    season = names(sites_per_season)[i],      # use the season names from sites_per_season
    species = names(naive_occurrence_per_year[[i]]),
    count = naive_occurrence_per_year[[i]],
    stringsAsFactors = TRUE
  )
})) %>% mutate(prob = count / sites_per_season[season])

species_metadata = read.csv("data/traits/species_metadata(included).csv", nrows = 107) %>% select(common_name, scientific_name, home_range_radius_m, residency, habitat_association, habitat_ebird) %>% mutate(species = tolower(common_name))
ggplot(left_join(naive_occurrence_per_year_df, species_metadata, by = "species") %>% filter(season != "2022"),
       aes(x = season, y = prob, color = habitat_association, group = species)) +
  geom_line() +
  geom_point() +
  geom_text_repel(
    aes(label = species),
    nudge_x = 0.2,          # adjust horizontal position
    direction = "y",         # repel only vertically
    hjust = 0,
    segment.size = 0.2
  ) +
  theme_bw() +
  labs(
    x = "Season/Year",
    y = "Naive occurrence probability",
    title = "Species occurrence over time"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"  # hide legend when using labels
  )

min_sites_detected = 1
message("Excluding species that were detected at fewer than the minimum number of sites (", min_sites_detected, "):")
species_to_remove = names(total_occurrences[total_occurrences < min_sites_detected])
# TODO: incorporate these additional species determined to be absent via manual review above into the naive statistics
# TODO: do not remove species for which there are manual validations
print(species_to_remove)
ylist[species_to_remove] = NULL
species = names(ylist)

# Naive species detections
message(length(species), " species detected")
message("Most commonly detected species:")
# print(naive_occurrence %>% slice_max(prob, n=1) %>% pull(species) %>% as.character())
which(total_occurrences == max(total_occurrences))
message("Least commonly detected species:")
# print(naive_occurrence %>% slice_min(prob, n=1) %>% pull(species) %>% as.character())
which(total_occurrences == min(total_occurrences[total_occurrences > 0]))
# p = ggplot(naive_occurrence, aes(x = nsites, y = species)) +
#   geom_bar(stat = "identity") +
#   geom_vline(xintercept = mean(naive_occurrence$nsites), color = "blue") +
#   labs(title = "Species occurrence across sites", x = "Number of sites detected as occurrence", y = ""); print(p)

# Naive species occurrence
# naive_occurrence = naive_occurrence %>% filter(species %in% names(ylist))
# message("Naive occurrence rate: mean ", round(mean(naive_occurrence$prob),2),
#         ", min ", round(min(naive_occurrence$prob),2), ", max ", round(max(naive_occurrence$prob),2))
# p = ggplot(naive_occurrence, aes(x = prob, y = species)) +
#   geom_point(stat = "identity") +
#   geom_vline(xintercept = mean(naive_occurrence$prob), color = "blue") +
#   labs(title = "Naive species occurrence", x = "Proportion of sites detected as occurrence", y = ""); print(p)

# Naive species richness per site
# species_per_site = setNames(rep(0, length(sites)), sites)
# for (sp in ylist) {
#   detected = rowSums(sp, na.rm = TRUE) > 0
#   species_per_site[detected] = species_per_site[detected] + 1
# }
# species_per_site = data.frame(site = names(species_per_site), species_detected = as.vector(species_per_site))
# mean_species_detected = mean(species_per_site$species_detected)
# message("Naive species richness per site (assuming no false positive detections): mean ", round(mean(species_per_site$species_detected),2), ", range ", min(species_per_site$species_detected), "–", max(species_per_site$species_detected))
# p = ggplot(species_per_site, aes(x = species_detected)) +
#   geom_histogram(binwidth = 1) +
#   geom_vline(xintercept = mean_species_detected, color = "blue") +
#   labs(title = "Naive species richness per site", x = "Number of species detected", y = "Number of sites"); print(p)

############################################################################################################################
# Format observation detection-nondetection and covariate data for modeling as 3D arrays (site × survey × species)

# Observed (uncertain) detection-nondetection data
y = array(NA, dim = c(length(sites), length(surveys), length(seasons), length(species)),
          dimnames = list(site = sites, survey = surveys, season = seasons, species = species))
for (t in seasons) {
  for (i in species) {
    y[ , , t, i] = as.matrix(ylist[[i]][[t]])
  }
}

# Left-align data (moving any missing NA surveys to the right) to allow for direct indexing by number of surveys per site
y_unaligned = y
left_align_row = function(x) {
  non_na = x[!is.na(x)]
  c(non_na, rep(NA, length(x) - length(non_na)))
}
for (t in dimnames(y)[['season']]) {
  for (i in dimnames(y)[['species']]) {
    # Extract site × survey matrix for this season × species
    sp_season_mat <- y[, , t, i]
    
    # Left-align across surveys for each site
    sp_season_aligned <- t(apply(sp_season_mat, 1, left_align_row))
    
    # Restore dimnames
    dimnames(sp_season_aligned) <- dimnames(sp_season_mat)
    
    # Put back into aligned array
    y[, , t, i] <- sp_season_aligned
  }
}
n_surveys_per_site = list()
for (t in dimnames(y)[["season"]]) {
  # Count number of non-NA surveys per site
  n_surveys_per_site[[t]] = apply(!is.na(y[, , t, 1]), 1, sum)
}

x_yday_unaligned = x_yday
for (t in names(x_yday)) {
  season_mat <- x_yday[[t]]
  
  # Left-align across surveys for each site
  season_mat_aligned <- t(apply(season_mat, 1, left_align_row))
  
  # Restore dimnames
  dimnames(season_mat_aligned) <- dimnames(season_mat)
  
  # Put back into aligned list
  x_yday[[t]] <- season_mat_aligned
}

# Observed (certain) detection data, i.e. "site confirmation design"
if (model_name == "fp_Miller") {
  path_in = paste0("/Users/giojacuzzi/repos/few-shot-transfer-learning-bioacoustics/data/cache", "/predictions_", model, ".parquet")
  predictions = arrow::read_parquet(path_in) %>% arrange(file, label_predicted)
  predictions = predictions %>%
    mutate(
      serialno = str_extract(file, "SMA\\d+"),
      season = str_extract(file, "SMA\\d+_(\\d{4})") %>% str_remove("^SMA\\d+_"),
      yday = str_extract(file, "SMA\\d+_(\\d{8})") %>% str_remove("^SMA\\d+_") %>% as.Date(format = "%Y%m%d") %>% yday()
    )
  path_site_key = "data/sites/site_key_long.csv"
  site_key = read.csv(path_site_key) %>% mutate(site = tolower(site), site_agg = tolower(site_agg), date = as.Date(date, format = "%m/%d/%y"), season = as.character(year(date)))
  site_key$yday = yday(site_key$date)
  predictions = predictions %>%
    left_join(
      site_key %>% select(serialno, yday, site, site_agg, season),
      by = c("serialno", "yday", "season")
    )
  predictions = predictions %>% select(label_truth, season, serialno, yday, site, site_agg)
  message("Number of certain detections per species:")
  print(table(predictions$label_truth))
  message("By season:")
  print(table(predictions$season, predictions$label_truth))
  
  ## TODO: Miller model incorporation!
  
  # --- 1. Map site names to row indices in y ---
  site_idx <- match(predictions$site_agg, dimnames(y)[["site"]])
  
  # --- 2. Map species names to species dimension indices in y ---
  species_idx <- match(predictions$label_truth, dimnames(y)[["species"]])
  
  # --- 3. Map seasons to season indices in y ---
  season_idx <- match(predictions$season, dimnames(y)[["season"]])
  
  # --- 3. Map yday to survey column indices using x_yday (season-aware) ---
  survey_idx <- mapply(function(site_name, yday_val, season_val) {
    # Skip if site not in this season's matrix
    if (!site_name %in% rownames(x_yday[[season_val]])) return(NA)
    
    # Find survey columns matching this yday
    cols <- which(x_yday[[season_val]][site_name, ] == yday_val)
    if (length(cols) == 0) return(NA)  # no match
    cols[1]  # take first if multiple matches
  },
  predictions$site_agg,
  predictions$yday,
  predictions$season)
  
  # --- 4. Keep only valid rows ---
  keep <- !is.na(site_idx) & !is.na(survey_idx) &
          !is.na(species_idx) & !is.na(season_idx) &
          predictions$label_truth != "0"

  # --- 5. Vectorized assignment into 4D y ---
  y[cbind(site_idx[keep], survey_idx[keep], season_idx[keep], species_idx[keep])] <- 2
  
  counts_df <- data.frame(
    species = species,
    nondetection = integer(length(species)),
    unconfirmed = integer(length(species)),
    confirmed = integer(length(species))
  )
  
  
  # # Loop over species and count
  # for (i in seq_along(species)) {
  #   vals <- as.vector(y[ , , 1, i])       # flatten site x survey for this species
  #   tab <- table(factor(vals, levels = 0:2))  # ensure all levels 0,1,2 are counted
  #   counts_df[i, c("nondetection", "unconfirmed", "confirmed")] <- as.integer(tab)
  # }
  # counts_df
  
  for (i in seq_along(dimnames(y)[["species"]])) {
    # flatten over sites × surveys × seasons
    vals <- as.vector(y[ , , , i])
    # ensure all levels (0,1,2) are represented
    tab <- table(factor(vals, levels = 0:2))
    counts_df[i, c("nondetection", "unconfirmed", "confirmed")] <- as.integer(tab)
  }
  counts_df
  
}

# Get detection covariate data
# detection_data = readRDS("data/cache/detection_covariates/data_detection.rds") %>% filter(year == "2020")

# x_yday_df = as.data.frame(x_yday[["2020"]])
# x_yday_df$site = rownames(x_yday_df)
# x_yday_long = tidyr::pivot_longer(
#   x_yday_df,
#   cols = -site,
#   names_to = "survey",
#   values_to = "yday"
# )
# 
# get_var_matrix = function(variable) {
#   detection_data_long = x_yday_long %>%
#     left_join(detection_data, by = c("site", "yday")) %>%
#     select(site, survey, !!sym(variable))
#   x = tidyr::pivot_wider(
#     detection_data_long,
#     names_from = survey,
#     values_from = !!sym(variable)
#   )
#   x = as.data.frame(x)
#   rownames(x) = x$site
#   x$site = NULL
#   return(as.matrix(x))
# }
# x_tmax = get_var_matrix("tmax_deg_c")
# x_prcp = get_var_matrix("prcp_mm_day")

# Get detection covariate data
detection_data = readRDS("data/cache/detection_covariates/data_detection.rds") %>% mutate(year = as.character(year))

# Convert x_yday list into a long dataframe with season info
x_yday_long_all <- lapply(names(x_yday), function(season) {
  x_df <- as.data.frame(x_yday[[season]])
  x_df$site <- rownames(x_df)
  tidyr::pivot_longer(
    x_df,
    cols = -site,
    names_to = "survey",
    values_to = "yday"
  ) %>% 
    mutate(
      season = as.character(season),
      survey = as.integer(survey)  # ensures numeric ordering
    )
}) %>% bind_rows()

# Ensure detection_data$year is character to match season
detection_data <- detection_data %>%
  mutate(year = as.character(year))

# Function to create site × survey × season array for a given covariate
get_var_array <- function(variable) {
  # Join covariate values to x_yday_long_all
  cov_long <- x_yday_long_all %>%
    left_join(
      detection_data,
      by = c("site" = "site", "yday" = "yday", "season" = "year")
    ) %>%
    select(site, survey, season, !!sym(variable))
  
  # Define dimensions
  sites   <- sort(unique(cov_long$site))
  surveys <- sort(unique(cov_long$survey))
  seasons <- sort(unique(cov_long$season))
  
  # Initialize array
  cov_array <- array(
    NA,
    dim = c(length(sites), length(surveys), length(seasons)),
    dimnames = list(site = sites, survey = surveys, season = seasons)
  )
  
  # Fill array
  for (row in seq_len(nrow(cov_long))) {
    s  <- cov_long$site[row]
    sv <- cov_long$survey[row]
    se <- cov_long$season[row]
    cov_array[s, sv, se] <- cov_long[[variable]][row]
  }
  
  return(cov_array)
}

# Example: create arrays for temperature and precipitation
x_tmax <- get_var_array("tmax_deg_c")
x_prcp <- get_var_array("prcp_mm_day")

############################################################################################################################
# Assemble occurrence covariate data

occ_data_plot = occ_data_plot %>% filter(site %in% dimnames(y)$site) # discard data for irrelevant sites
stopifnot(dimnames(y)[["site"]] == occ_data_plot$site) # check that covariate data are aligned with observation matrix by site

param_alpha_names = c(
  "elev",
  "homerange_htmax_cv",
  "homerange_qmd_mean",
  "homerange_treeden_all_mean"
)
param_delta_names = c(
  # "cover_forest_diversity",
  "density_roads_paved",
  "density_streams_major",
  "focalpatch_area_homeange_pcnt",
  "prop_abund_comthin",
  "prop_abund_lsog",
  "prop_abund_standinit",
  "shape_idx"
)

# Store alpha parameter ID, variable name, and standardize data to have mean 0, standard deviation 1
param_alpha_data = tibble(param = paste0("alpha", 1:length(param_alpha_names)), name  = param_alpha_names)
param_alpha_data = param_alpha_data %>% rowwise() %>% mutate(scaled = list(scale(occ_data_plot[[name]]))) %>% ungroup()
n_alpha_params = nrow(param_alpha_data)

# Standardize species names of occ homerange data
names(occ_data_homerange)[names(occ_data_homerange) == "western flycatcher"] = "pacific-slope flycatcher"

# Store minimum (i.e. floor) plot homerange data for those species who have smaller estimated home range sizes
occ_data_homerange_floor = occ_data_homerange[['plot']] %>% filter(site %in% dimnames(y)$site)
stopifnot(dimnames(y)$site == occ_data_homerange_floor$site)

# Store delta parameter ID, variable name, and standardized data in species-specific matricies
delta_data = setNames(vector('list', length(param_delta_names)), param_delta_names)
for (param in param_delta_names) {
  message(param)
  
  species_delta_data = setNames(vector('list', length(species)), species)
  
  for (i in 1:length(species)) {
    species_name = species[i]
    # message(i, " ", species_name)
    # discard data for irrelevant sites
    species_occ_data = occ_data_homerange[[species_name]] %>% filter(site %in% dimnames(y)$site)
    # check that data are aligned with observation matrix by site
    stopifnot(dimnames(y)$site == species_occ_data$site)
    if (unique(species_occ_data %>% pull(buffer_radius_m)) >= 100) {
      param_data = species_occ_data %>% select(site, all_of(param))
    } else {
      param_data = occ_data_homerange_floor %>% select(site, all_of(param))
    }
    species_delta_data[[species_name]] = scale(param_data[[param]])
  }
  delta_data[[param]] = species_delta_data
}
param_delta_data = tibble(param = paste0("delta", 1:length(param_delta_names)), name  = param_delta_names)
param_delta_data = param_delta_data %>% rowwise() %>% mutate(data = list(delta_data[[name]])) %>% ungroup()
n_delta_params = nrow(param_delta_data)

# Assemble detection covariate data
detect_data = list(
  yday        = x_yday,
  prcp_mm_day = x_prcp,
  tmax_deg_c  = x_tmax
)
# Store beta parameter ID, variable name, and standardize data to have mean 0, standard deviation 1
param_beta_data = tibble(param = paste0("beta", seq_along(detect_data)), name = names(detect_data))

# param_beta_data = param_beta_data %>% rowwise() %>% mutate(scaled = list(scale(as.vector(detect_data[[name]])))) %>% ungroup()
param_beta_data = param_beta_data %>% rowwise() %>% mutate(scaled = list(scale(unlist(detect_data[[name]], use.names = FALSE)))) %>% ungroup()

n_beta_params = nrow(param_beta_data)

param_season_data = tibble(param = "season", name = "season", scaled = list(scale(1:length(seasons))))

############################################################################################################################
# Prepare all data for the model

# Initialize latent occupancy state z[i] as 1 if a detection occurred at site i, and 0 otherwise
z = array(NA, dim = c(length(sites), length(seasons), length(species)), dimnames = list(sites, seasons, species))
for (j in seq_along(sites)) {
  for (t in seq_along(seasons)) {
    for (i in seq_along(species)) {
      z[j, t, i] = (sum(y[j, , t, i], na.rm = TRUE) > 0) * 1
    }
  }
}

# Model data constants and covariates
I = length(species)
J = length(sites)
K = as.matrix(as.data.frame(lapply(n_surveys_per_site, as.vector)))
T = length(seasons)

msom_data = list(
  y = y, # observed (detection-nondetection) data matrix
  J = J, # number of sites sampled
  K = K, # number of secondary sampling periods (surveys) per site per season (site x season)
  T = T, # number of primary sampling periods (seasons)
  I = I, # number of species observed
  eps = 1e-6
)
for (a in seq_len(n_alpha_params)) { # Add alpha covariates
  msom_data[[paste0("x_", param_alpha_data$param[a])]] <- as.vector(param_alpha_data$scaled[[a]])
}
for (d in seq_len(n_delta_params)) { # Add delta covariates
  param_data_new = do.call(cbind, param_delta_data %>% filter(name == param_delta_data$name[d]) %>% pull(data) %>% .[[1]])
  msom_data[[paste0("x_", param_delta_data$param[d])]] = param_data_new
}
msom_data[["x_season"]] = as.vector(param_season_data$scaled[[1]])

to_3d_array <- function(x) {
  if (is.list(x)) {
    # yday case: list of season-specific matrices
    arr <- abind::abind(x, along = 3)
    dimnames(arr)[[3]] <- names(x)  # keep season names
    return(arr)
  } else {
    # already an array (prcp_mm_day, tmax_deg_c)
    return(x)
  }
}

# Loop over beta parameters
for (b in seq_len(n_beta_params)) {
  arr <- to_3d_array(detect_data[[param_beta_data$name[b]]])  # now always 3D
  vec <- as.vector(param_beta_data$scaled[[b]])               # flatten to 1D
  
  stopifnot(length(vec) == prod(dim(arr)))  # safety check
  
  msom_data[[paste0("x_", param_beta_data$param[b])]] <- array(
    vec,
    dim = dim(arr),
    dimnames = dimnames(arr)
  )
}

str(msom_data)

############################################################################################################################
# Following:
# - https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.1365-2664.2009.01664.x
# - https://esajournals.onlinelibrary.wiley.com/doi/10.1002/eap.2293
# - https://www.sciencedirect.com/science/article/abs/pii/S0006320709004819
model_template = paste0("
model{

  ## Community level hyperpriors

  # Occurrence
  psi.mean ~ dunif(0,1)                   # probability scale
  mu.u <- log(psi.mean) - log(1-psi.mean) # logit scale
  sigma.u ~ dunif(0,5)                    # standard deviation
  tau.u <- pow(sigma.u,-2)                # precision
  
  # Covariate effects on occurrence
  mu.alpha1 ~ dnorm(0,0.01)
  sigma.alpha1 ~ dunif(0,5)
  tau.alpha1 <- pow(sigma.alpha1,-2)
  mu.alpha2 ~ dnorm(0,0.01)
  sigma.alpha2 ~ dunif(0,5)
  tau.alpha2 <- pow(sigma.alpha2,-2)
  mu.alpha3 ~ dnorm(0,0.01)
  sigma.alpha3 ~ dunif(0,5)
  tau.alpha3 <- pow(sigma.alpha3,-2)
  mu.alpha4 ~ dnorm(0,0.01)
  sigma.alpha4 ~ dunif(0,5)
  tau.alpha4 <- pow(sigma.alpha4,-2)
  mu.delta1 ~ dnorm(0,0.01)
  sigma.delta1 ~ dunif(0,5)
  tau.delta1 <- pow(sigma.delta1,-2)
  mu.delta2 ~ dnorm(0,0.01)
  sigma.delta2 ~ dunif(0,5)
  tau.delta2 <- pow(sigma.delta2,-2)
  mu.delta3 ~ dnorm(0,0.01)
  sigma.delta3 ~ dunif(0,5)
  tau.delta3 <- pow(sigma.delta3,-2)
  mu.delta4 ~ dnorm(0,0.01)
  sigma.delta4 ~ dunif(0,5)
  tau.delta4 <- pow(sigma.delta4,-2)
  mu.delta5 ~ dnorm(0,0.01)
  sigma.delta5 ~ dunif(0,5)
  tau.delta5 <- pow(sigma.delta5,-2)
  mu.delta6 ~ dnorm(0,0.01)
  sigma.delta6 ~ dunif(0,5)
  tau.delta6 <- pow(sigma.delta6,-2)
  mu.delta7 ~ dnorm(0,0.01)
  sigma.delta7 ~ dunif(0,5)
  tau.delta7 <- pow(sigma.delta7,-2)
  mu.season ~ dnorm(0,0.01)
  sigma.season ~ dunif(0,5)
  tau.season <- pow(sigma.season,-2)
  
  # Detection (unconfirmed true-positive)
  v.mean  ~ dunif(0,1)
  mu.v  <- log(v.mean) - log(1-v.mean)
  sigma.v ~ dunif(0,5)
  tau.v <- pow(sigma.v,-2)

  # Covariate effects on detection (unconfirmed true-positive)
  mu.beta1    ~ dnorm(0,0.01)
  sigma.beta1 ~ dunif(0,5)
  tau.beta1  <- pow(sigma.beta1,-2)
  mu.beta2    ~ dnorm(0,0.01)
  sigma.beta2 ~ dunif(0,5)
  tau.beta2  <- pow(sigma.beta2,-2)
  mu.beta3    ~ dnorm(0,0.01)
  sigma.beta3 ~ dunif(0,5)
  tau.beta3  <- pow(sigma.beta3,-2)
  
  # False-positive unconfirmed detection (Royle and Link 2006, Miller et al. 2011)
  # mu.w ~ dnorm(0,0.01)
  # sigma.w ~ dunif(0,5)
  # tau.w <- pow(sigma.w,-2)
  
  # Confirmed true positive detection (Miller et al. 2011)
  # mu.b ~ dnorm(0,0.01)
  # sigma.b ~ dunif(0,5)
  # tau.b <- pow(sigma.b,-2)

  for (i in 1:I) { # for each species
  
      ## Species level priors
      
      # Occupancy
      u[i]      ~ dnorm(mu.u, tau.u)
      alpha1[i] ~ dnorm(mu.alpha1,tau.alpha1)
      alpha2[i] ~ dnorm(mu.alpha2,tau.alpha2)
      alpha3[i] ~ dnorm(mu.alpha3,tau.alpha3)
      alpha4[i] ~ dnorm(mu.alpha4,tau.alpha4)
      delta1[i] ~ dnorm(mu.delta1,tau.delta1)
      delta2[i] ~ dnorm(mu.delta2,tau.delta2)
      delta3[i] ~ dnorm(mu.delta3,tau.delta3)
      delta4[i] ~ dnorm(mu.delta4,tau.delta4)
      delta5[i] ~ dnorm(mu.delta5,tau.delta5)
      delta6[i] ~ dnorm(mu.delta6,tau.delta6)
      delta7[i] ~ dnorm(mu.delta7,tau.delta7)
      season[i] ~ dnorm(mu.season,tau.season)
  
      # Unconfirmed true positive detection
      v[i]     ~ dnorm(mu.v, tau.v)
      beta1[i] ~ dnorm(mu.beta1,tau.beta1)
      beta2[i] ~ dnorm(mu.beta2,tau.beta2)
      beta3[i] ~ dnorm(mu.beta3,tau.beta3)
      
      # Unconfirmed false positive detection
      # w[i] ~ dnorm(mu.w,tau.w)
      
      # Confirmed true positive detection
      # gamma[i] ~ dnorm(mu.b, tau.b)
      # b[i] <- 1 / (1 + exp(-gamma[i]))
  
      for (j in 1:J) { # for each site j
      
          for (t in 1:T) { # for each season t
          
              # Ecological process model for latent occurrence z
              logit(psi[j,t,i]) <- u[i] + alpha1[i]*x_alpha1[j] + alpha2[i]*x_alpha2[j] + alpha3[i]*x_alpha3[j] + alpha4[i]*x_alpha4[j] + delta1[i]*x_delta1[j,i] + delta2[i]*x_delta2[j,i] + delta3[i]*x_delta3[j,i] + delta4[i]*x_delta4[j,i] + delta5[i]*x_delta5[j,i] + delta6[i]*x_delta6[j,i] + delta7[i]*x_delta7[j,i] + season[i]*x_season[t]
              z[j,t,i] ~ dbern(psi[j,t,i])
              
              for (k in 1:K[j,t]) { # for each sampling period (survey) k at site j during season t
    
                  ## Observation model, assumming no false positives (e.g. Zipkin et al. 2009)
                  
                  # p11 is true-positive detection probability given z = 1
                  logit(p11[j,k,t,i]) <- v[i] + beta1[i]*x_beta1[j,k,t] + beta2[i]*x_beta2[j,k,t] + beta3[i]*x_beta3[j,k,t]
                
                  p[j,k,t,i] <- z[j,t,i] * p11[j,k,t,i]
                  
                  # Observed outcome
                  y[j,k,t,i] ~ dbern(p[j,k,t,i])
                  
                  # Simulated replicate for posterior predictive checks
                  y.sim[j,k,t,i] ~ dbern(p[j,k,t,i])

                  # Deviance (Broms et al. 2016)
                  d.obs.dev[j,k,t,i] <- -2 * ( y[j,k,t,i] * log(p[j,k,t,i]) + (1 - y[j,k,t,i]) * log(1 - p[j,k,t,i]) )
                  d.sim.dev[j,k,t,i] <- -2 * ( y.sim[j,k,t,i] * log(p[j,k,t,i]) + (1 - y.sim[j,k,t,i]) * log(1 - p[j,k,t,i]) )
      
              } # K surveys
              
              ## Sums per site/species for each posterior predictive check
              # Bernoulli deviance contribution
              d.obs.dev.sum[j,i] <- sum(d.obs.dev[j,1:K[j, ], ,i])
              d.sim.dev.sum[j,i] <- sum(d.sim.dev[j,1:K[j, ], ,i])
          
          } # T seasons

      } # J sites
  } # I species
  
  ## Derived quantities
  
  # Discrepancy measure between observed and simulated data is defined as mean(D.obs > D.sim)
  D.obs.dev <- sum(d.obs.dev.sum[1:J,1:I])
  D.sim.dev <- sum(d.sim.dev.sum[1:J,1:I])
  bayes.p.dev <- step(D.sim.dev - D.obs.dev)

  # Estimated number of occuring sites per species per season (among the sampled population of sites)
  for (i in 1:I) {
    for (t in 1:T) {
      Nocc[i,t] <- sum(z[ , t, i])
    }
  }

  # Estimated number of occuring species per site per season (among the species that were detected anywhere)
  for (j in 1:J) {
    for (t in 1:T) {
      Nsite[j,t] <- sum(z[j, t, ])
    }
  }
}
")
model_spec = model_template
cat(strsplit(model_spec, "\n")[[1]], sep = "\n") # print model specification to console
model_file = tempfile()
writeLines(model_spec, con = model_file)

############################################################################################################################
# Run JAGS

message("\n", "System CPU: "); print(as.data.frame(t(benchmarkme::get_cpu())))
message("System RAM: "); print(benchmarkme::get_ram())

message("Running JAGS (current time ", time_start <- Sys.time(), ")")

msom = jags(data = msom_data,
            inits = function() { list( # initial values to avoid data/model conflicts
              z = z,
              v = rep(logit(0.70), length(species)),
              w = rep(logit(0.05), length(species))
            ) },
            parameters.to.save = c( # monitored parameters
              "psi", "z",
              "mu.u", "sigma.u", "u",
              "mu.v", "sigma.v", "v",
              # "mu.w", "sigma.w", "w",
              # "mu.b", "sigma.b", "b",
              paste0("mu.alpha", 1:n_alpha_params), paste0("sigma.alpha", 1:n_alpha_params), paste0("alpha", 1:n_alpha_params),
              paste0("mu.delta", 1:n_delta_params), paste0("sigma.delta", 1:n_delta_params), paste0("delta", 1:n_delta_params),
              paste0("mu.beta",  1:n_beta_params),  paste0("sigma.beta",  1:n_beta_params),  paste0("beta",  1:n_beta_params),
              "mu.season", "sigma.season", "season",
              "d.obs.dev", "d.sim.dev", # "bayes.p.se",
              "Nsite", "Nocc"
            ),
            model.file = model_file,
            n.chains = 2, n.adapt = 100, n.iter = 2000, n.burnin = 1000, n.thin = 1,
            parallel = TRUE, DIC = FALSE)

message("Finished running JAGS (", round(as.numeric(difftime(Sys.time(), time_start, units = 'mins')), 2), " minutes)")

############################################################################################################################
# Retrieve summary data and investigate goodness-of-fit

msom_summary = summary(msom)
msom_summary = msom_summary %>% as_tibble() %>%
  mutate(param = rownames(summary(msom)), overlap0 = as.factor(overlap0)) %>% relocate(param, .before = 1) %>%
  mutate(prob = plogis(mean), prob_lower95 = plogis(`2.5%`), prob_upper95 = plogis(`97.5%`))
rhat_threshold = 1.1
suspected_nonconvergence = msom_summary %>% filter(Rhat >= rhat_threshold) %>% filter(!str_starts(param, "z\\[") & !str_starts(param, "psi\\["))
suspected_nonconvergence = suspected_nonconvergence %>% mutate(
  index = str_extract(param, "(?<=\\[)\\d+(?=\\])"),
  index = as.integer(index),
  species = ifelse(!is.na(index), species[index], NA)
)
if (nrow(suspected_nonconvergence) > 1) {
  message("The following ", nrow(suspected_nonconvergence), " parameters may not have converged:")
  print(suspected_nonconvergence)
} else {
  message("All parameters appear to have converged (rhat < ", rhat_threshold, ")")
}

## Posterior predictive checks
# "If the observed data are consistent with the model in question, then the Bayesian p-value should be close to 0.5. In practice, a p-value close to 0 or 1 indicates that the model is inadequate in some way -- close to 0 suggests a lack of fit and close to 1 suggests that the model over-fits the data, which may occur when it is too complex." (MacKenzie et al. 2018)

# Squared error (Zipkin et al. 2009)
# Does the model predict observed detections/non-detections about as well as it predicts data simulated under its own assumptions?
# Summing these gives a measure of average squared prediction error across all sites, surveys, and species.
# This Bayesian p-value is the probability (proportion of iterations) that the discrepancy for simulated data exceeds that for observed data.
message("Baysian p-value (squared error):")
print(mean(msom$sims.list$bayes.p.se))
# Bernoulli deviance contribution (Broms et al. 2016)
# Is the overall likelihood of the observed detection histories under the fitted model about the same as the likelihood of new data generated from that model?
# This Bayesian p-value is the probability (proportion of iterations) that the simulated deviance is greater than the observed deviance.
message("Baysian p-value (deviance residuals):")
print(mean(msom$sims.list$bayes.p.dev))
{
  # Extract the posterior samples
  d.obs.dev <- msom$sims.list$d.obs.dev
  d.sim.dev <- msom$sims.list$d.sim.dev
  
  # Get the dimensions
  # Typically: sims x J x Kmax x T x I
  dim(d.obs.dev)  # [n.sims, J, Kmax, T, I]
  
  n.sims <- dim(d.obs.dev)[1]
  
  # Suppose you have K[j,t] giving the number of surveys at site j, season t
  # Compute D.obs.dev and D.sim.dev for each posterior sample
  D.obs.dev <- numeric(n.sims)
  D.sim.dev <- numeric(n.sims)
  
  for (s in 1:n.sims) {
    total.obs <- 0
    total.sim <- 0
    for (j in 1:J) {
      for (t in 1:T) {
        Kjt <- K[j,t]
        for (i in 1:I) {
          total.obs <- total.obs + sum(d.obs.dev[s, j, 1:Kjt, t, i])
          total.sim <- total.sim + sum(d.sim.dev[s, j, 1:Kjt, t, i])
        }
      }
    }
    D.obs.dev[s] <- total.obs
    D.sim.dev[s] <- total.sim
  }
  
  # Compute Bayesian p-value
  # > 0.5, model underfit (simulated deviance larger than observed)
  # < 0.5, model overfit  (observed deviance larger than simulated)
  bayes.p.dev <- mean(D.sim.dev > D.obs.dev)
  bayes.p.dev
}


# Get posterior samples (i.e. MCMC simulated draws) for estimated occurrence psi and latent state z
psi_samples = msom$sims.list$psi
z_samples   = msom$sims.list$z
# These should be of dimension: samples x sites x seasons x species
n_samples = dim(psi_samples)[1]
for (s in list(psi_samples, z_samples)) {
  stopifnot(identical(dim(s)[2], J))
  stopifnot(identical(dim(s)[3], T))
  stopifnot(identical(dim(s)[4], I))
}

auc_combined = rep(NA, n_samples)
pb = progress_bar$new(format = "[:bar] :percent :elapsedfull (ETA :eta)", total = n_samples, clear = FALSE)
for (s in 1:n_samples) {
  psi_all = as.vector(psi_samples[s, , , ]) # combine all species into a single vector
  z_all   = as.vector(z_samples[s, , , ])
  
  if (length(unique(z_all)) > 1) {
    auc_combined[s] = as.numeric(pROC::auc(pROC::roc(z_all, psi_all, quiet=TRUE)))
  } else {
    stop("Cannot calculate ROC") # latent z states are identical for this sample, something is wrong
  }
  pb$tick()
}
auc_combined_mean = round(mean(auc_combined, na.rm=TRUE), 3)
auc_combined_bci  = round(quantile(auc_combined, probs=c(0.025, 0.975), na.rm=TRUE), 3)
auc_combined = data.frame(
  mean_auc = auc_combined_mean,
  lower95_auc = auc_combined_bci[[1]],
  upper95_auc = auc_combined_bci[[2]]
)

message("Combined mean AUC: ", auc_combined_mean, " (95% BCI ", auc_combined_bci[1], "–", auc_combined_bci[2], ")")

auc_species = matrix(NA, nrow=n_samples, ncol=I)
pb = progress_bar$new(format = "[:bar] :percent :elapsedfull (ETA :eta)", total = n_samples, clear = FALSE)
for (s in 1:n_samples) {
  for (i in 1:I) {
    psi_vec = as.vector(psi_samples[s, , , i])
    z_vec   = as.vector(z_samples[s, , , i])
    
    if (length(unique(z_vec)) > 1) {
      auc_species[s, i] = as.numeric(pROC::auc(pROC::roc(z_vec, psi_vec, quiet=TRUE)))
    } else {
      # species[i] latent z state is identical at all sites for this sample, cannot calculate ROC
    }
  }
  pb$tick()
}
auc_species_mean = round(apply(auc_species, 2, mean, na.rm=TRUE), 3)
auc_species_bci  = round(apply(auc_species, 2, quantile, probs=c(0.025, 0.975), na.rm=TRUE), 3)
# nocc_samples = msom$sims.list$Nocc # posterior mean number of occupied sites (z=1) per species (i.e. effective positive sample size)
# nocc_mean = round(apply(nocc_samples, 2, mean),1)
# nocc_bci = apply(nocc_samples, 2, quantile, probs = c(0.025, 0.975))
auc_species = data.frame(
  species = species[1:I],
  mean_auc = auc_species_mean,
  lower95_auc = auc_species_bci["2.5%", ],
  upper95_auc = auc_species_bci["97.5%", ]
  # mean_nocc = nocc_mean,
  # lower95_nocc = nocc_bci["2.5%", ],
  # upper95_nocc = nocc_bci["97.5%", ]
)

message("Species-specific mean AUC: ", round(mean(auc_species$mean_auc, na.rm = TRUE),3),
        " (range ",  round(min(auc_species$mean_auc, na.rm = TRUE),3), "–", round(max(auc_species$mean_auc, na.rm = TRUE),3), ")")
print(auc_species)

# Write results to cache
msom_results = list(
  model_spec = model_spec,
  msom_summary = msom_summary,
  p.se        = mean(msom$sims.list$bayes.p.se),
  p.dev = mean(msom$sims.list$bayes.p.dev),
  auc_combined = auc_combined,
  auc_species  = auc_species,
  param_alpha_data = param_alpha_data,
  param_delta_data = param_delta_data,
  param_season_data = param_season_data,
  param_beta_data = param_beta_data,
  sites = sites,
  species = species
)
if (!dir.exists(dirname(path_out))) dir.create(dirname(path_out), recursive = TRUE)
saveRDS(msom_results, file = path_out)
message(crayon::green("Cached model and results to", path_out))
