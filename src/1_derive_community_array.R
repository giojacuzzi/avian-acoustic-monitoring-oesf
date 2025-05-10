####################################################################################
# Derive a 5-dimensional observation matrix (unit, survey, species, season) from aggregated prediction data
#
# Input:
# - Species list used to produce raw prediction data
# - Cached dataframe of recording counts (i.e. prediction files) per unit-date
# - Cached dataframe of all predictions
path_species_list = "models/ensemble/ensemble_species_list.txt"
path_recording_counts = "data/cache/0_aggregate_raw_prediction_data/recording_counts.feather"
path_prediction_data = "data/cache/0_aggregate_raw_prediction_data/prediction_data.feather"
#
# Output:
# - A 4-dimensional array with dimensions [unit × survey × season × species], with each element representing an observation (e.g. max confidence)
path_out_community_array = "data/cache/1_derive_community_array/community_array.rds"
####################################################################################

library(tidyverse)
library(arrow)
library(lubridate)

# Load species list of common names
species_list = sapply(strsplit(readLines(path_species_list), "_"), `[`, 2)

# Load prediction data
message("Loading aggregated prediction data from ", path_prediction_data, "...")
prediction_data = read_feather(path_prediction_data)

# Inspect prediction data
str(prediction_data)

# Clean prediction data
message("Cleaning prediction data...")
prediction_data_filtered = prediction_data

# Manually exclude specific survey(s)
prediction_data_filtered = prediction_data_filtered %>% filter(!(season == 2020 & deploy == 2 & serialno == "SMA00403"))  # Incomplete and redone in a later deployment
prediction_data_filtered = prediction_data_filtered %>% filter(!(season == 2022 & deploy == 3 & serialno == "SMA08215")) # Conducted as a replacement for an earlier deployment that was at first unsuccessfully retrieved and presumed lost
prediction_data_filtered = prediction_data_filtered %>% mutate(serialno = as.character(serialno)) %>% filter(!startsWith(serialno, "SM2")) # SM2 model ARUs have a substantially different construction than SM4 and Mini models, and were used for a negligible number of surveys

predicitions_with_missing_units = prediction_data_filtered %>% filter(is.na(unit))
if (nrow(predicitions_with_missing_units) > 0) {
  warning("Predictions with missing units:")
  print(predicitions_with_missing_units %>% distinct(unit, season, deploy, serialno), n = 1000)
}

# Remove predictions that are not associated with a sampling unit
prediction_data_filtered = prediction_data_filtered %>% filter(!is.na(unit))

message("Deriving species observations per season, unit, and survey...")

# Define primary and secondary sampling periods
prediction_data_filtered = prediction_data_filtered %>%
  mutate(
    season = factor(year(time)), # Define a "survey" (secondary period) as a continuous 24h recording period (from 00:00:00 to 23:59:59)
    survey_date = as.Date(time)  # Define a "season" (primary period) by year
  )

# Group predictions by unique unit-survey combinations, filtering for species with a minimum number of predictions
prediction_data_filtered = prediction_data_filtered %>%
  group_by(unit, survey_date, common_name) %>%
  filter(n() > 1) %>%
  ungroup()

# Filter out incomplete surveys (i.e. unit-survey combinations that did not record all 24h)
recording_counts = read_feather(path_recording_counts)
prediction_data_filtered = prediction_data_filtered %>%
  left_join(recording_counts, by = c("serialno", "survey_date")) %>%
  filter(recordings == 24) %>%
  select(-recordings)

# TODO: Obtain weighted confidence scores from classification ensemble here

# Summarize all predictions per unit-survey to a single observation per species using the maximum confidence value
unit_survey_max_conf = prediction_data_filtered %>%
  mutate(conf = coalesce(confidence, confidence_source, confidence_target)) %>%
  group_by(unit, season, survey_date, common_name) %>%
  summarise(max_conf = max(conf, na.rm = TRUE), .groups = "drop")

# TODO: Do any thresholding or probability calculation here

# Add any missing species per unit-survey
unit_survey_max_conf = unit_survey_max_conf %>%
  group_by(unit, survey_date) %>%
  complete(common_name = species_list, fill = list(max_conf = 0.0)) %>% # NOTE: As the original CNN output was limited to scores >= 0.5 to conserve memory, we are manually setting any predictions with score < 0.5 to 0.0 here
  ungroup() %>%
  mutate(season = as.factor(format(survey_date, "%Y")))

# Count surveys (i.e. "visits") for each unit per season, accounting for gaps in dates
unit_survey_max_conf = unit_survey_max_conf %>%
  group_by(unit, season) %>%
  mutate(
    min_date = min(survey_date),
    survey = as.integer(difftime(survey_date, min_date, units = "days")) + 1
  ) %>%
  ungroup() %>%
  select(-min_date)
summary(unit_survey_max_conf$survey)

# Standardize number of surveys according to the maximum, adding any missing surveys per species with confidence NA
complete_surveys = unit_survey_max_conf %>% distinct(unit, season, common_name) %>% crossing(survey = 1:max(unit_survey_max_conf$survey))
unit_survey_max_conf = complete_surveys %>% left_join(unit_survey_max_conf, by = c("unit", "season", "common_name", "survey"))

message("Deriving community observation array...")
# Initialize the 4D community observation array
units   = sort(unique(unit_survey_max_conf$unit))
surveys = sort(unique(unit_survey_max_conf$survey))
seasons = sort(unique(unit_survey_max_conf$season))
species = sort(unique(unit_survey_max_conf$common_name))
community_array = array(
  NA,
  dim = c(length(units), length(surveys), length(seasons), length(species)),
  dimnames = list(unit = units, survey = surveys, season = seasons, common_name = species)
)
# Populate the array
for (season in seasons) {
  for (name in species) {
    season_name_data = unit_survey_max_conf[unit_survey_max_conf$season == season & unit_survey_max_conf$common_name == name, ]
    season_name_matrix = tapply(season_name_data$max_conf, list(season_name_data$unit, season_name_data$survey), FUN = identity)
    community_array[, , season, name] = as.matrix(season_name_matrix[units, as.character(surveys)])
  }
}

message("Unit surveys per season:")
print(apply(community_array, c(1, 3), function(x) { sum(apply(!is.na(x), 1, any)) }))

message("Units per season:")
print(apply(community_array, 3, function(x) { sum(apply(x, 1, function(u) any(!is.na(u)))) }))

# Example data access
# View(community_array[, , "2020", "Barred Owl"])

# Write results to cache
message("Caching community array...")
if (!dir.exists(dirname(path_out_community_array))) dir.create(dirname(path_out_community_array), recursive = TRUE)
saveRDS(community_array, path_out_community_array)
