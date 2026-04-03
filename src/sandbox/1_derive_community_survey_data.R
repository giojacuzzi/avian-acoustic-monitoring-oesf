####################################################################################
# Derive a multidimensional observation array (site, survey, species, season) from aggregated prediction data
#
# INPUT:
# Species list used to produce raw prediction data
path_species_list = "models/ensemble/ensemble_species_list.txt"
# Cached dataframe of prediction file counts (i.e. recordings) per site-survey
path_survey_file_counts = "data/cache/0_aggregate_raw_prediction_data/survey_file_counts.feather"
# Cached dataframe of all predictions
path_prediction_data = "data/cache/0_aggregate_raw_prediction_data/prediction_data.feather"
# Minimum number of predictions in a survey (i.e. secondary sampling period) to retain
min_prediction_count = 1
#
# OUTPUT:
# A multidimensional data structure with dimensions [site × survey × season × species], with each
# element containing a named list of observational data (including survey date and confidence scores)
path_out_community_survey_data = paste0("data/cache/1_derive_community_survey_data/community_survey_data_", format(Sys.Date(), "%Y-%m-%d"), ".rds")
####################################################################################

message("Deriving community survey data (current time ", time_start <- format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ")")

library(tidyverse)

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
  common_name = tolower(common_name),
  site = as.factor(tolower(as.character(unit))),
  unit_agg = as.factor(tolower(as.character(unit_agg))),
  confidence_source = replace_na(confidence_source, 0.0),
  confidence_target = replace_na(confidence_target, 0.0)
)
survey_file_counts_filtered = survey_file_counts %>% mutate(
  site = as.factor(tolower(as.character(unit)))
)

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
complete_surveys = survey_file_counts_filtered %>%
  filter(n_prediction_files == 24)
message("Discarding predictions from ", nrow(survey_file_counts_filtered) - nrow(complete_surveys), " incomplete surveys (i.e. < 24 hr)")
prediction_data_filtered_complete = prediction_data_filtered %>%
  semi_join(complete_surveys, by = c("site", "survey_date"))
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

message("Deriving species observations per season, site, and survey")

# Initialize community observation array
species = tolower(sort(sapply(strsplit(readLines(path_species_list), "_"), `[`, 2)))
seasons = sort(unique(survey_file_counts_filtered$season))
sites   = sort(unique(survey_file_counts_filtered$site)) # Matrix rows
surveys = 1:max(setNames(max_survey_period_per_season$max_survey_period, as.character(max_survey_period_per_season$season))) # Matrix columns

community_survey_data = vector("list", length(sites) * length(surveys) * length(seasons) * length(species))
dim(community_survey_data) = c(length(sites), length(surveys), length(seasons), length(species))
dimnames(community_survey_data) = list(
  site = sites,
  survey = as.character(surveys),
  season = seasons,
  common_name = species
)

# Initialize all site-survey elements that were actually surveyed
# (i.e. present in survey_file_counts_filtered) with a list
# containing survey date and 0.0 confidence
for (i in seq_len(nrow(survey_file_counts_filtered))) {
  row      = survey_file_counts_filtered[i, ]
  unit_i   = as.character(row$site)
  survey_i = as.character(row$survey)
  season_i = as.character(row$season)

  for (species_i in dimnames(community_survey_data)$common_name) {
    community_survey_data[[unit_i, survey_i, season_i, species_i]] = list(
      survey_date = row$survey_date,
      confidence = c(0.0),
      confidence_source = c(0.0),
      confidence_target = c(0.0)
    )
  }
}

# Example: Get site-survey matrix for a given season and species
(slice = community_survey_data[, , "2022", "barred owl", drop = FALSE])
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
         site = as.character(site),
         season = as.character(season),
         survey = as.character(survey))

# Populate `community_survey_data` with confidence scores
prediction_data_with_survey_edit %>%
  group_by(site, survey, season, common_name) %>%
  summarise(conf_vector = list(confidence), .groups = "drop") %>% # Create a list of confidence values for each site-survey-season-common_name
  rowwise() %>% # For each site-survey-season-common_name with its corresponding list of confidence values...
  mutate(
    updated = { # Update the corresponding community_survey_data list element
      element <- community_survey_data[[site, survey, season, common_name]]  # Get the current list element
      element$confidence <- conf_vector                                      # Overwrite "confidence" vector 
      community_survey_data[[site, survey, season, common_name]] <<- element # Overwrite the list element
      TRUE # Populate the `updated` column
    }
  )
prediction_data_with_survey_edit %>%
  group_by(site, survey, season, common_name) %>%
  summarise(conf_vector = list(confidence_source), .groups = "drop") %>% # Create a list of confidence values for each site-survey-season-common_name
  rowwise() %>% # For each site-survey-season-common_name with its corresponding list of confidence values...
  mutate(
    updated = { # Update the corresponding community_survey_data list element
      element <- community_survey_data[[site, survey, season, common_name]]  # Get the current list element
      element$confidence_source <- conf_vector                              # Overwrite "confidence_source" vector
      community_survey_data[[site, survey, season, common_name]] <<- element # Overwrite the list element
      TRUE # Populate the `updated` column
    }
  )
prediction_data_with_survey_edit %>%
  group_by(site, survey, season, common_name) %>%
  summarise(conf_vector = list(confidence_target), .groups = "drop") %>% # Create a list of confidence values for each site-survey-season-common_name
  rowwise() %>% # For each site-survey-season-common_name with its corresponding list of confidence values...
  mutate(
    updated = { # Update the corresponding community_survey_data list element
      element <- community_survey_data[[site, survey, season, common_name]]  # Get the current list element
      element$confidence_target <- conf_vector                              # Overwrite "confidence_source" vector
      community_survey_data[[site, survey, season, common_name]] <<- element # Overwrite the list element
      TRUE # Populate the `updated` column
    }
  )

# TODO: Overwrite confidence scores for manually verified data (i.e. 100% confidence)

# Example data retrieval: Get observation matrix
slice_after = community_survey_data[, , "2022", "barred owl", drop = FALSE]
options(max.print = 1e6)
matrix( # Survey date
  unlist(lapply(slice_after, function(x) if (!is.null(x)) x$survey_date else NA)),
  nrow = dim(slice_after)[1],
  ncol = dim(slice_after)[2],
  dimnames = dimnames(slice_after)[1:2]
)
matrix( # Max confidence (combined)
  unlist(lapply(slice_after, function(x) if (!is.null(x)) max(x$confidence, na.rm = TRUE) else NA)),
  nrow = dim(slice_after)[1],
  ncol = dim(slice_after)[2],
  dimnames = dimnames(slice_after)[1:2]
)
matrix( # Max confidence (source)
  unlist(lapply(slice_after, function(x) if (!is.null(x)) max(x$confidence_source, na.rm = TRUE) else NA)),
  nrow = dim(slice_after)[1],
  ncol = dim(slice_after)[2],
  dimnames = dimnames(slice_after)[1:2]
)
matrix( # Max confidence (target)
  unlist(lapply(slice_after, function(x) if (!is.null(x)) max(x$confidence_target, na.rm = TRUE) else NA)),
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

# Examine the data structure
print(dimnames(community_survey_data))

# Write results to cache
if (!dir.exists(dirname(path_out_community_survey_data))) dir.create(dirname(path_out_community_survey_data), recursive = TRUE)
saveRDS(community_survey_data, path_out_community_survey_data)
message(crayon::green("Cached community survey data to ", path_out_community_survey_data))

message(crayon::green("Finished deriving community survey data (", round(as.numeric(difftime(Sys.time(), time_start, units = 'mins')), 2), " min )"))
