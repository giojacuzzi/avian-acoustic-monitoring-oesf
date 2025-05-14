####################################################################################
# Derive a multidimensional observation array (unit, survey, species, season) from aggregated prediction data
#
# INPUT:
# Species list used to produce raw prediction data
path_species_list = "models/ensemble/ensemble_species_list.txt"
# Cached dataframe of prediction file counts (i.e. recordings) per unit-survey
path_out_survey_file_counts = "data/cache/0_aggregate_raw_prediction_data/survey_file_counts.feather"
# Cached dataframe of all predictions
path_prediction_data = "data/cache/0_aggregate_raw_prediction_data/prediction_data.feather"
# Minimum number of predictions in a survey (i.e. secondary sampling period) to retain
min_prediction_count = 1
#
# OUTPUT:
# A multidimensional data structure with dimensions [unit × survey × season × species], with each
# element containing a named list of observational data (including survey date and confidence scores)
path_out_community_survey_data = "data/cache/1_derive_community_array/community_survey_data.rds"
####################################################################################

library(tidyverse)

# Load survey file counts
message("Loading survey file counts from ", path_out_survey_file_counts)
survey_file_counts = arrow::read_feather(path_out_survey_file_counts) %>%
  group_by(season, unit) %>%
  mutate(survey = as.integer(survey_date - min(survey_date)) + 1) %>%
  ungroup()

# Load prediction data
message("Loading aggregated prediction data from ", path_prediction_data)
prediction_data = arrow::read_feather(path_prediction_data) %>% mutate(survey_date = as.Date(time))

## Clean data
message("Cleaning data")
prediction_data_filtered = prediction_data
survey_file_counts_filtered = survey_file_counts

## Manually exclude specific survey(s)
# Incomplete and redone in a later deployment
survey_file_counts_filtered = survey_file_counts_filtered %>% filter(!(season == 2020 & deploy == 2 & serialno == "SMA00403"))
prediction_data_filtered = prediction_data_filtered %>% filter(!(season == 2020 & deploy == 2 & serialno == "SMA00403"))
# Conducted as a replacement for an earlier deployment that was at first unsuccessfully retrieved and presumed lost
survey_file_counts_filtered = survey_file_counts_filtered %>% filter(!(season == 2022 & deploy == 3 & serialno == "SMA08215"))
prediction_data_filtered = prediction_data_filtered %>% filter(!(season == 2022 & deploy == 3 & serialno == "SMA08215"))
# SM2 model ARUs have a substantially different construction than SM4 and Mini models, and were used for a negligible number of surveys
survey_file_counts_filtered = survey_file_counts_filtered %>% mutate(serialno = as.character(serialno)) %>% filter(!startsWith(serialno, "SM2"))
prediction_data_filtered = prediction_data_filtered %>% mutate(serialno = as.character(serialno)) %>% filter(!startsWith(serialno, "SM2"))

# TODO: If needed, filter for predictions only within a specific timeframe

# Remove predictions that are not associated with a sampling unit
predicitions_with_missing_units = prediction_data_filtered %>% filter(is.na(unit))
if (nrow(predicitions_with_missing_units) > 0) {
  warning("Discarding predictions with missing units:")
  print(predicitions_with_missing_units %>% distinct(unit, season, deploy, serialno))
  prediction_data_filtered = prediction_data_filtered %>% filter(!is.na(unit))
}

# Group predictions by unique unit-survey combinations, filtering for species with a minimum number of predictions
prediction_data_filtered = prediction_data_filtered %>%
  group_by(unit, survey_date, common_name) %>%
  filter(n() >= min_prediction_count) %>%
  ungroup()

# Filter out incomplete surveys (i.e. unit-survey combinations that did not record all 24h)
complete_surveys = survey_file_counts_filtered %>%
  filter(n_prediction_files == 24)
message("Discarding predictions from ", nrow(survey_file_counts_filtered) - nrow(complete_surveys), " incomplete surveys (i.e. < 24 hr)")
prediction_data_filtered_complete = prediction_data_filtered %>%
  semi_join(complete_surveys, by = c("unit", "survey_date"))
survey_file_counts_filtered = complete_surveys

# # Determine maximum number of surveys conducted at a unit for each season
# max_surveys_per_season = prediction_data_filtered_complete %>%
#   distinct(season, unit, survey_date) %>%  # Unique surveys per unit per season
#   count(season, unit) %>%                  # Count surveys per unit within each season
#   group_by(season) %>%
#   summarise(max_surveys = max(n))

message("Calculating survey period and numbers")

# Determine maximum survey period at a unit for each season (i.e. days from the first to last survey date, inclusive)
# This range will form the columns in the observation matrix for a given species-season
max_survey_period_per_season = survey_file_counts_filtered %>%
  distinct(season, unit, survey_date) %>%      # Unique survey dates per unit per season
  group_by(season, unit) %>%
  summarise(
    survey_period = as.integer(max(survey_date) - min(survey_date)) + 1, .groups = "drop"
  ) %>%
  group_by(season) %>%
  summarise(max_survey_period = max(survey_period))  # Max period per season

# Join survey numbers with prediction data
prediction_data_with_survey = prediction_data_filtered_complete %>%
  left_join(
    survey_file_counts_filtered %>%
      select(season, deploy, serialno, unit, survey_date, survey),
    by = c("season", "deploy", "serialno", "unit", "survey_date")
  )

# TODO: Obtain weighted confidence scores from classification ensemble here

message("Deriving species observations per season, unit, and survey")

# Initialize community observation array
species = sort(sapply(strsplit(readLines(path_species_list), "_"), `[`, 2))
seasons = sort(unique(survey_file_counts_filtered$season))
units   = sort(unique(survey_file_counts_filtered$unit)) # Matrix rows
surveys = 1:max(setNames(max_survey_period_per_season$max_survey_period, as.character(max_survey_period_per_season$season))) # Matrix columns

community_survey_data = vector("list", length(units) * length(surveys) * length(seasons) * length(species))
dim(community_survey_data) = c(length(units), length(surveys), length(seasons), length(species))
dimnames(community_survey_data) = list(
  unit = units,
  survey = as.character(surveys),
  season = seasons,
  common_name = species
)

# Initialize all unit-survey elements that were actually surveyed
# (i.e. present in survey_file_counts_filtered) with a list
# containing survey date and 0.0 confidence
for (i in seq_len(nrow(survey_file_counts_filtered))) {
  row      = survey_file_counts_filtered[i, ]
  unit_i   = as.character(row$unit)
  survey_i = as.character(row$survey)
  season_i = as.character(row$season)

  for (species_i in dimnames(community_survey_data)$common_name) {
    community_survey_data[[unit_i, survey_i, season_i, species_i]] = list(
      survey_date = row$survey_date,
      confidence = c(0.0)
    )
  }
}

# Example: Get unit-survey matrix for a given season and species
(slice = community_survey_data[, , "2022", "Barred Owl", drop = FALSE])
options(max.print = 1e6)
matrix( # Survey date
  unlist(lapply(slice, function(x) if (!is.null(x)) x$survey_date else NA)),
  nrow = dim(slice)[1],
  ncol = dim(slice)[2],
  dimnames = dimnames(slice)[1:2]
)
matrix( # Max confidence
  unlist(lapply(slice, function(x) if (!is.null(x)) max(x$confidence, na.rm = TRUE) else NA)),
  nrow = dim(slice)[1],
  ncol = dim(slice)[2],
  dimnames = dimnames(slice)[1:2]
)

# Ensure common_name is a character to match species names in `community_survey_data`
prediction_data_with_survey_edit = prediction_data_with_survey %>%
  mutate(common_name = as.character(common_name),
         unit = as.character(unit),
         season = as.character(season),
         survey = as.character(survey))

# Populate `community_survey_data` with confidence scores
prediction_data_with_survey_edit %>%
  group_by(unit, survey, season, common_name) %>%
  summarise(conf_vector = list(confidence), .groups = "drop") %>% # Create a list of confidence values for each unit-survey-season-common_name
  rowwise() %>% # For each unit-survey-season-common_name with its corresponding list of confidence values...
  mutate(
    updated = { # Update the corresponding community_survey_data list element
      element <- community_survey_data[[unit, survey, season, common_name]]  # Get the current list element
      element$confidence <- conf_vector                                      # Overwrite "confidence" vector 
      community_survey_data[[unit, survey, season, common_name]] <<- element # Overwrite the list element
      TRUE # Populate the `updated` column
    }
  )

# Example data retrieval: Get observation matrix
slice_after = community_survey_data[, , "2022", "Barred Owl", drop = FALSE]
options(max.print = 1e6)
matrix( # Survey date
  unlist(lapply(slice_after, function(x) if (!is.null(x)) x$survey_date else NA)),
  nrow = dim(slice_after)[1],
  ncol = dim(slice_after)[2],
  dimnames = dimnames(slice_after)[1:2]
)
matrix( # Max confidence
  unlist(lapply(slice_after, function(x) if (!is.null(x)) max(x$confidence, na.rm = TRUE) else NA)),
  nrow = dim(slice_after)[1],
  ncol = dim(slice_after)[2],
  dimnames = dimnames(slice_after)[1:2]
)
matrix( # Thresholded
  unlist(lapply(slice_after, function(x) if (!is.null(x)) as.integer(any(x$confidence > 0.9, na.rm = TRUE)) else NA)),
  nrow = dim(slice_after)[1],
  ncol = dim(slice_after)[2],
  dimnames = dimnames(slice_after)[1:2]
)

# Write results to cache
message("Caching community survey data to ", path_out_community_survey_data)
if (!dir.exists(dirname(path_out_community_survey_data))) dir.create(dirname(path_out_community_survey_data), recursive = TRUE)
saveRDS(community_survey_data, path_out_community_survey_data)
