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
# Output:
# - A 4-dimensional array with dimensions [unit × survey × season × species], with each element representing an observation (e.g. max confidence)
path_out_community_array = "data/cache/1_derive_community_array/community_array.rds"
####################################################################################

library(tidyverse)
library(arrow)
library(lubridate)

# Load survey file counts
message("Loading survey file counts from ", path_out_survey_file_counts, "...")
survey_file_counts = read_feather(path_out_survey_file_counts) %>%
  group_by(season, unit) %>%
  mutate(survey = as.integer(survey_date - min(survey_date)) + 1) %>%
  ungroup()

# Load prediction data
message("Loading aggregated prediction data from ", path_prediction_data, "...")
prediction_data = read_feather(path_prediction_data) %>% mutate(survey_date = as.Date(time))

## Clean data
message("Cleaning data...")
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

# TODO: If needed, filter for predictions only in a specific timeframe

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

message("Deriving species observations per season, unit, and survey...")

# Initialize community observation array
species = sort(sapply(strsplit(readLines(path_species_list), "_"), `[`, 2))
seasons = sort(unique(survey_file_counts_filtered$season))
units   = sort(unique(survey_file_counts_filtered$unit)) # Matrix rows
surveys = 1:max(setNames(max_survey_period_per_season$max_survey_period, as.character(max_survey_period_per_season$season))) # Matrix columns

community_array <- vector("list", length(units) * length(surveys) * length(seasons) * length(species))
dim(community_array) <- c(length(units), length(surveys), length(seasons), length(species))
dimnames(community_array) <- list(
  unit = units,
  survey = as.character(surveys),
  season = seasons,
  common_name = species
)

# Loop over each row in the tibble and populate matching entries in the list-array
for (i in seq_len(nrow(survey_file_counts_filtered))) {
  row <- survey_file_counts_filtered[i, ]
  unit_i <- as.character(row$unit)
  survey_i <- as.character(row$survey)
  season_i <- as.character(row$season)

  for (species_i in dimnames(community_array)$common_name) { # Initialize all surveyed sites with 0.0 for each species
    community_array[[unit_i, survey_i, season_i, species_i]] <- list(
      survey_date = row$survey_date,
      confidence = c(0.0)
    )
  }
}

# EX: Get observation matrix
slice = community_array[, , "2022", "Barred Owl", drop = FALSE]
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

# Make sure `common_name` is character to match species names in community_array
prediction_data_with_survey_edit <- prediction_data_with_survey %>%
  mutate(common_name = as.character(common_name),
         unit = as.character(unit),
         season = as.character(season),
         survey = as.character(survey))

# Iterate over each row group by unique unit/survey/season/species
prediction_data_with_survey_edit %>%
  group_by(unit, survey, season, common_name) %>%
  summarise(conf_vector = list(confidence), .groups = "drop") %>%
  rowwise() %>%
  mutate(
    # For each group, update the corresponding community_array cell
    updated = {
      # Access current list at array position
      entry <- community_array[[unit, survey, season, common_name]]
      
      # If entry is NULL, create a new list
      if (is.null(entry)) {
        entry <- list()
      }
      
      # Overwrite or add the "confidence" element
      entry$confidence <- conf_vector
      
      # Assign back to array
      community_array[[unit, survey, season, common_name]] <<- entry
      
      TRUE  # placeholder to make mutate work
    }
  )

# EX: Get observation matrix
slice = community_array[, , "2022", "Barred Owl", drop = FALSE]
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
matrix( # Thresholded
  unlist(lapply(slice, function(x) if (!is.null(x)) as.integer(any(x$confidence > 0.9, na.rm = TRUE)) else NA)),
  nrow = dim(slice)[1],
  ncol = dim(slice)[2],
  dimnames = dimnames(slice)[1:2]
)








# community_array = array(
#   NA,
#   dim = c(length(units), length(surveys), length(seasons), length(species)),
#   dimnames = list(unit = units, survey = as.character(surveys), season = seasons, common_name = species)
# )

# data_structure <- list()
# for (sp in species) {
#   print(sp)
#   for (sn in seasons) {
#     print(sn)
#     for (un in units) {
#       print(un)
#       for (sv in surveys) {
#         # TODO: For surveys that were conducted (i.e. present in survey_file_counts_filtered) but had no predictions,
#         # Give a default confidence history of 0.0
#         if (survey_file_counts_filtered %>% filter(season == sn, unit == un, survey == sv) %>% nrow() > 0) {
#           data_structure[[sp]][[sn]][[un]][[sv]] = list(
#             survey_date = (survey_file_counts_filtered %>% filter(season == sn, unit == un, survey == sv))[["survey_date"]],
#             confidence = c(0)
#           )
#         } else {
#           data_structure[[sp]][[sn]][[un]][[sv]] = list()
#         }
#         # if (prediction_data_with_survey %>% filter(common_name == common_name, season == season, unit == unit, survey == survey) %>% nrow() > 0) {
#         #   (data_structure[[common_name]][[season]][[unit]][[survey]])[[confidence]] = (prediction_data_with_survey %>% filter(common_name == common_name, season == season, unit == unit, survey == survey))[["confidence"]]
#         # }
#       }
#     }
#     # community_array[, , season, common_name] =
#       # season_name_data = unit_survey_max_conf[unit_survey_max_conf$season == season & unit_survey_max_conf$common_name == name, ]
#       # season_name_matrix = tapply(season_name_data$max_conf, list(season_name_data$unit, season_name_data$survey), FUN = identity)
#       # community_array[, , season, name] = as.matrix(season_name_matrix[units, as.character(surveys)])
#   }
# }
# 
# 
# 
# survey_matrix <- prediction_data_with_survey %>%
#   filter(common_name == "Barred Owl", season == "2020") %>%
#   select(unit, survey, survey_date, confidence) %>%
#   distinct() %>%
#   mutate(
#     unit = factor(unit, levels = units),
#     survey = as.integer(survey)
#   ) %>%
#   group_by(unit, survey) %>%
#   summarise(
#     value = list(list(
#       survey_date = unique(survey_date)[1],
#       confidence = confidence
#     )),
#     .groups = "drop"
#   ) %>%
#   pivot_wider(
#     names_from = survey,
#     values_from = value,
#     values_fill = list(list()), # Empty list for units that were surveyed, but have no predictions
#     names_sort = TRUE
#   ) %>%
#   complete(unit = units, fill = list(NULL)) %>%  # NULL values for all units that were not surveyed
#   arrange(unit)
# 
# missing_surveys <- setdiff(as.character(surveys), names(survey_matrix))
# for (s in missing_surveys) {
#   survey_matrix[[s]] <- rep(list(NULL), nrow(survey_matrix))
# }
# 
# survey_matrix %>% select(-unit) %>%
#   mutate(across(everything(), ~map(.x, ~if (is.list(.x)) .x$survey_date else NA))) %>%
#   mutate(across(everything(), ~ ifelse(is.logical(.), NA, as.numeric(.))))
# 
# 
# ##############################################################
# 
# 
# 
# message("Deriving community observation array...")
# # Initialize the 4D community observation array
# 
# community_array = array(
#   NA,
#   dim = c(length(units), length(surveys), length(seasons), length(species)),
#   dimnames = list(unit = units, survey = as.character(surveys), season = seasons, common_name = species)
# )
# # Populate the array
# for (season in seasons) {
#   for (name in species) {
#     community_array[, , season, name] =
#     # season_name_data = unit_survey_max_conf[unit_survey_max_conf$season == season & unit_survey_max_conf$common_name == name, ]
#     # season_name_matrix = tapply(season_name_data$max_conf, list(season_name_data$unit, season_name_data$survey), FUN = identity)
#     # community_array[, , season, name] = as.matrix(season_name_matrix[units, as.character(surveys)])
#   }
# }
# 
# 
# # filtered_data <- prediction_data_filtered %>%
# #   filter(common_name == "Barred Owl", season == "2022") %>%
# #   select(unit, survey, survey_date) %>%
# #   distinct()  # remove duplicates, if any
# # survey_matrix <- filtered_data %>%
# #   pivot_wider(names_from = survey, values_from = survey_date, names_sort = T)
# 
# # survey_matrix <- prediction_data_filtered %>%
# #   filter(common_name == "Barred Owl", season == "2020") %>%
# #   select(unit, survey, survey_date) %>%
# #   distinct() %>%
# #   mutate(
# #     unit = factor(unit, levels = units),  # enforce unit row order
# #   ) %>%
# #   pivot_wider(
# #     names_from = survey,
# #     values_from = survey_date,
# #     values_fn = list,         # <- store values in a list
# #     values_fill = NA,   # <- fill empty cells with NA in a list
# #     names_sort = TRUE
# #   ) %>%
# #   arrange(unit)
# 
# survey_matrix <- prediction_data_filtered %>%
#   filter(common_name == "Barred Owl", season == "2022") %>%
#   select(unit, survey, survey_date, confidence) %>%
#   distinct() %>%
#   mutate(
#     unit = factor(unit, levels = units),
#     survey = as.integer(survey)
#   ) %>%
#   group_by(unit, survey) %>%
#   summarise(
#     value = list(list(
#       survey_date = unique(survey_date)[1],
#       confidence = confidence
#     )),
#     .groups = "drop"
#   ) %>%
#   pivot_wider(
#     names_from = survey,
#     values_from = value,
#     values_fill = list(),
#     names_sort = TRUE
#   ) %>%
#   complete(unit = units, fill = list(value = list())) %>%  # Ensure all units appear
#   arrange(unit)
# 
# # Get survey_dates
# fdsa = survey_matrix %>% select(-unit) %>%
#   mutate(across(everything(), ~map(.x, ~if (is.list(.x)) .x$survey_date else NA))) %>%
#   mutate(across(everything(), ~ ifelse(is.logical(.), NA, as.numeric(.))))
# fdsa_mat = as.matrix(fdsa)
# 
# test = survey_matrix %>% select(-unit) %>%
#   mutate(across(everything(), ~ map_chr(., ~ ifelse(length(.x) > 0, as.character(.x[1]), NA))))
# 
# test %>%
#   mutate(across(everything(), ~ ifelse(is.na(.), NA, as.Date(as.character(.), origin = "1970-01-01"))))
# 
# 
# 
# 
# 
# 
# nested_data <- prediction_data_filtered %>%
#   mutate(survey = as.character(survey)) %>%
#   group_by(unit, survey, season, common_name) %>%
#   summarise(
#     survey_date = list(survey_date),
#     confidence = list(confidence),
#     .groups = "drop"
#   )
# 
# # Create full set of dimension combinations
# full_grid <- expand.grid(
#   unit = unique(prediction_data_filtered$unit),
#   survey = as.character(unique(prediction_data_filtered$survey)),
#   season = unique(prediction_data_filtered$season),
#   common_name = unique(prediction_data_filtered$common_name),
#   stringsAsFactors = FALSE
# )
# 
# # Join with nested data and replace NAs with empty lists
# full_joined <- full_grid %>%
#   left_join(nested_data, by = c("unit", "survey", "season", "common_name")) %>%
#   mutate(
#     survey_date = map(survey_date, ~ .x %||% as.POSIXct(character())),
#     confidence = map(confidence, ~ .x %||% numeric()),
#     combined = map2(survey_date, confidence, ~ list(survey_date = .x, confidence = .y))
#   )
# 
# # Turn into a 4D array
# community_array <- array(
#   full_joined$combined,
#   dim = c(
#     length(unique(full_grid$unit)),
#     length(unique(full_grid$survey)),
#     length(unique(full_grid$season)),
#     length(unique(full_grid$common_name))
#   ),
#   dimnames = list(
#     unit = unique(full_grid$unit),
#     survey = unique(full_grid$survey),
#     season = unique(full_grid$season),
#     common_name = unique(full_grid$common_name)
#   )
# )
