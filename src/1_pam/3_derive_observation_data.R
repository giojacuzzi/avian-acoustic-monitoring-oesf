# 3_derive_observation_data.R ####################################################################################
# Derive a multidimensional observation array (site, survey, species, season) from aggregated prediction data
#
# CONFIG:
min_prediction_count = 1 # Minimum number of predictions in a survey (i.e. secondary sampling period) to retain
#
# OUTPUT:
# A multidimensional data structure with dimensions [site × survey × season × species], with each
# element containing a named list of observational data (including survey date and confidence scores)
path_out_community_array_predictions = "data/cache/1_pam/3_derive_observation_data/community_array_predictions.rds"
path_out_community_array_surveydates = "data/cache/1_pam/3_derive_observation_data/community_array_surveydates.rds"
#
# INPUT:
# Cached dataframe of all predictions
path_prediction_data    = "data/cache/1_pam/2_agg_raw_predictions/prediction_data.feather"
# Cached dataframe of prediction file counts (i.e. recordings) per site-survey
path_survey_file_counts = "data/cache/1_pam/2_agg_raw_predictions/survey_file_counts.feather"
##################################################################################################################

# TODO: REVISE THE CODE BELOW

source("src/global.R")

message("Deriving community survey data (current time ", time_start <- format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ")")

# Load survey file counts
message("Loading survey file counts from ", path_survey_file_counts)
survey_file_counts = arrow::read_feather(path_survey_file_counts) %>%
  group_by(season, unit) %>%
  mutate(survey = as.integer(survey_date - min(survey_date)) + 1) %>%
  ungroup()

# Load prediction data
message("Loading aggregated prediction data from ", path_prediction_data)
prediction_data = arrow::read_feather(path_prediction_data) %>% mutate(survey_date = as.Date(time))

## Clean data
message("Cleaning data")
prediction_data_filtered = prediction_data %>% mutate(
  common_name       = tolower(common_name),
  site              = as.factor(tolower(as.character(unit))),
  unit_agg          = as.factor(tolower(as.character(unit_agg))),
  confidence_source = replace_na(confidence_source, 0.0),
  confidence_target = replace_na(confidence_target, 0.0)
)
survey_file_counts_filtered = survey_file_counts %>% mutate(
  site = as.factor(tolower(as.character(unit)))
)

## Manually exclude specific survey(s)
# Incomplete and redone in a later deployment
survey_file_counts_filtered = survey_file_counts_filtered %>% filter(!(season == 2020 & deploy == 2 & serialno == "SMA00403"))
prediction_data_filtered    = prediction_data_filtered    %>% filter(!(season == 2020 & deploy == 2 & serialno == "SMA00403"))
# Conducted as a replacement for an earlier deployment that was at first unsuccessfully retrieved and presumed lost
survey_file_counts_filtered = survey_file_counts_filtered %>% filter(!(season == 2022 & deploy == 3 & serialno == "SMA08215"))
prediction_data_filtered    = prediction_data_filtered    %>% filter(!(season == 2022 & deploy == 3 & serialno == "SMA08215"))
# SM2 model ARUs have a substantially different construction than SM4 and Mini models, and were used for a negligible number of surveys
survey_file_counts_filtered = survey_file_counts_filtered %>% mutate(serialno = as.character(serialno)) %>% filter(!startsWith(serialno, "SM2"))
prediction_data_filtered    = prediction_data_filtered    %>% mutate(serialno = as.character(serialno)) %>% filter(!startsWith(serialno, "SM2"))

# Remove predictions that are not associated with a sampling site
predicitions_with_missing_units = prediction_data_filtered %>% filter(is.na(site))
if (nrow(predicitions_with_missing_units) > 0) {
  warning(crayon::yellow("Discarding predictions with missing sites:"))
  print(crayon::yellow(predicitions_with_missing_units %>% distinct(site, season, deploy, serialno)))
  prediction_data_filtered = prediction_data_filtered %>% filter(!is.na(site))
}

# Group predictions by unique site-survey combinations, filtering for species with a minimum number of predictions
prediction_data_filtered = prediction_data_filtered %>%
  group_by(site, survey_date, common_name) %>%
  filter(n() >= min_prediction_count) %>%
  ungroup()

# Filter out incomplete surveys (i.e. site-survey combinations that did not record all 24h)
complete_surveys = survey_file_counts_filtered %>% filter(n_prediction_files == 24)
message("Discarding predictions from ", nrow(survey_file_counts_filtered) - nrow(complete_surveys), " incomplete surveys (i.e. < 24 hr)")
prediction_data_filtered_complete = prediction_data_filtered %>% semi_join(complete_surveys, by = c("site", "survey_date"))
survey_file_counts_filtered = complete_surveys

# # Determine maximum number of surveys conducted at a site for each season
# max_surveys_per_season = prediction_data_filtered_complete %>%
#   distinct(season, site, survey_date) %>%  # Unique surveys per site per season
#   count(season, site) %>%                  # Count surveys per site within each season
#   group_by(season) %>%
#   summarise(max_surveys = max(n))

message("Calculating survey period and numbers")

# Determine maximum survey period at a site for each season (i.e. days from the first to last survey date, inclusive)
# This range will form the columns in the observation matrix for a given species-season
max_survey_period_per_season = survey_file_counts_filtered %>%
  distinct(season, site, survey_date) %>%      # Unique survey dates per site per season
  group_by(season, site) %>%
  summarise(
    survey_period = as.integer(max(survey_date) - min(survey_date)) + 1, .groups = "drop"
  ) %>%
  group_by(season) %>%
  summarise(max_survey_period = max(survey_period))  # Max period per season

# Join survey numbers with prediction data
prediction_data_with_survey = prediction_data_filtered_complete %>%
  left_join(
    survey_file_counts_filtered %>%
      select(season, deploy, serialno, site, survey_date, survey),
    by = c("season", "deploy", "serialno", "site", "survey_date")
  )

message("Deriving species predictions and survey dates per season, site, and survey")

# Initialize community observation array
species = class_labels$common_name
seasons = sort(unique(survey_file_counts_filtered$season))
sites   = sort(unique(survey_file_counts_filtered$site)) # Matrix rows
surveys = 1:max(setNames(max_survey_period_per_season$max_survey_period, as.character(max_survey_period_per_season$season))) # Matrix columns

community_array_predictions = vector("list", length(sites) * length(surveys) * length(seasons) * length(species))
dim(community_array_predictions) = c(length(sites), length(surveys), length(seasons), length(species))
dimnames(community_array_predictions) = list(
  site = sites,
  survey = as.character(surveys),
  season = seasons,
  common_name = species
)
community_array_surveydates = community_array_predictions


# Initialize all site-survey elements that were actually surveyed
# (i.e. present in survey_file_counts_filtered) with a list
# containing survey date and 0.0 confidence
for (i in seq_len(nrow(survey_file_counts_filtered))) {
  row      = survey_file_counts_filtered[i, ]
  unit_i   = as.character(row$site)
  survey_i = as.character(row$survey)
  season_i = as.character(row$season)
  
  for (species_i in dimnames(community_array_predictions)$common_name) {
    community_array_predictions[[unit_i, survey_i, season_i, species_i]] = list(
      confidence_source = c(0.0),
      confidence_target = c(0.0)
    )
    community_array_surveydates[[unit_i, survey_i, season_i, species_i]] = list(
      survey_date = row$survey_date,
      yday = yday(row$survey_date)
    )
  }
}

# Example: Get site-survey matrix for a given season and species
(slice = community_array_surveydates[, , "2022", "barred owl", drop = FALSE])
options(max.print = 1e6)
matrix( # Survey date
  unlist(lapply(slice, function(x) if (!is.null(x)) x$survey_date else NA)),
  nrow = dim(slice)[1],
  ncol = dim(slice)[2],
  dimnames = dimnames(slice)[1:2]
)
matrix( # yday
  unlist(lapply(slice, function(x) if (!is.null(x)) x$yday else NA)),
  nrow = dim(slice)[1],
  ncol = dim(slice)[2],
  dimnames = dimnames(slice)[1:2]
)
(slice = community_array_predictions[, , "2022", "barred owl", drop = FALSE])
matrix( # Max source confidence
  unlist(lapply(slice, function(x) if (!is.null(x)) max(x$confidence_source, na.rm = TRUE) else NA)),
  nrow = dim(slice)[1],
  ncol = dim(slice)[2],
  dimnames = dimnames(slice)[1:2]
)

# Ensure common_name is a character to match species names in `community_array_predictions`
prediction_data_with_survey_edit = prediction_data_with_survey %>%
  mutate(common_name = as.character(common_name),
         site = as.character(site),
         season = as.character(season),
         survey = as.character(survey))

# Populate `community_array_predictions` with predicted classifier scores
prediction_data_with_survey_edit %>%
  group_by(site, survey, season, common_name) %>%
  summarise(conf_vector_source = list(confidence_source), conf_vector_target = list(confidence_target), .groups = "drop") %>% # Create a list of confidence values for each site-survey-season-common_name
  rowwise() %>% # For each site-survey-season-common_name with its corresponding list of confidence values...
  mutate(
    updated = { # Update the corresponding community_array_predictions list element
      element <- community_array_predictions[[site, survey, season, common_name]]  # Get the current list element
      element$confidence_source <- conf_vector_source # Overwrite "confidence" vectors
      element$confidence_target <- conf_vector_target                             
      community_array_predictions[[site, survey, season, common_name]] <<- element # Overwrite the list element
      TRUE # Populate the `updated` column
    }
  )

# Example data retrieval: Get observation matrix
slice_surveydates = community_array_surveydates[, , "2022", "barred owl", drop = FALSE]
slice_predictions = community_array_predictions[, , "2022", "barred owl", drop = FALSE]
options(max.print = 1e6)
matrix( # Survey date
  unlist(lapply(slice_surveydates, function(x) if (!is.null(x)) x$survey_date else NA)),
  nrow = dim(slice_surveydates)[1],
  ncol = dim(slice_surveydates)[2],
  dimnames = dimnames(slice_surveydates)[1:2]
)
matrix( # yday
  unlist(lapply(slice_surveydates, function(x) if (!is.null(x)) x$yday else NA)),
  nrow = dim(slice_surveydates)[1],
  ncol = dim(slice_surveydates)[2],
  dimnames = dimnames(slice_surveydates)[1:2]
)
matrix( # Max confidence (source)
  unlist(lapply(slice_predictions, function(x) if (!is.null(x)) max(x$confidence_source, na.rm = TRUE) else NA)),
  nrow = dim(slice_predictions)[1],
  ncol = dim(slice_predictions)[2],
  dimnames = dimnames(slice_predictions)[1:2]
)
matrix( # Max confidence (target)
  unlist(lapply(slice_predictions, function(x) if (!is.null(x)) max(x$confidence_target, na.rm = TRUE) else NA)),
  nrow = dim(slice_predictions)[1],
  ncol = dim(slice_predictions)[2],
  dimnames = dimnames(slice_predictions)[1:2]
)
matrix( # e.g. thresholded
  unlist(lapply(slice_after, function(x) if (!is.null(x)) as.integer(any(x$confidence_source > 0.9, na.rm = TRUE)) else NA)),
  nrow = dim(slice_after)[1],
  ncol = dim(slice_after)[2],
  dimnames = dimnames(slice_after)[1:2]
)

# Examine the data structure
print(dimnames(community_array_predictions))

# Write results to cache
if (!dir.exists(dirname(path_out_community_array_predictions))) dir.create(dirname(path_out_community_array_predictions), recursive = TRUE)
saveRDS(community_array_predictions, path_out_community_array_predictions)
message(crayon::green("Cached community predictions array to", path_out_community_array_predictions))
saveRDS(community_array_surveydates, path_out_community_array_surveydates)
message(crayon::green("Cached community survey dates array to", path_out_community_array_surveydates))

message(crayon::green("Finished aggregating community survey arrays (", round(as.numeric(difftime(Sys.time(), time_start, units = 'mins')), 2), " min )"))
