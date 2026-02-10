# 2_agg_raw_predictions.R ####################################################################################
# Aggregate raw prediction (.csv) data and associated metadata into a single table
# 
# OUTPUT:
path_out_dir = "data/cache/1_pam/2_agg_raw_predictions"
path_out_prediction_data    = paste0(path_out_dir, '/prediction_data.feather')    # All predictions
path_out_survey_file_counts = paste0(path_out_dir, '/survey_file_counts.feather') # Prediction file counts (i.e. recordings) per site-survey
#
# INPUT:
# Raw .csv prediction data with metadata contained in directory structure
path_in_dir = "/Volumes/gioj/OESF/predictions"
##############################################################################################################

source("src/global.R")

# Parse metadata from prediction .csv path
parse_metadata = function(file_path) {
  # Split the file path into metadata components
  metadata_components = strsplit(normalizePath(file_path), .Platform$file.sep)[[1]]
  
  # Extract metadata
  n = length(metadata_components)
  season = metadata_components[n - 4]
  deploy = metadata_components[n - 3]
  deploy = strsplit(deploy, 'Deployment')[[1]][2]
  serialno_date_time = metadata_components[n]
  serialno_date_time = strsplit(serialno_date_time, ".BirdNET.results.csv")[[1]]
  serialno = strsplit(serialno_date_time, '_')[[1]][1]
  date = strsplit(serialno_date_time, '_')[[1]][2]
  time = strsplit(serialno_date_time, '_')[[1]][3]
  rec_time_start = as.POSIXct(paste(date,time), format = "%Y%m%d %H%M%S", tz = "UTC")
  
  return(list(
    season = season,
    deploy = deploy,
    serialno = serialno,
    rec_time_start = rec_time_start
  ))
}


message("Locating all prediction files (i.e. recordings)")
survey_files = list.files(path = path_in_dir, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
message("Located ", length(survey_files), " survey files")

message("Parsing metadata for each prediction file")
metadata_list = lapply(seq_along(survey_files), function(i) {
  c(parse_metadata(survey_files[[i]]), file_path = survey_files[[i]])
})
predictions_metadata = data.table::rbindlist(metadata_list, fill = TRUE)

message("Matching sampling unit IDs")
unit_key = read.csv(path_unit_key) %>%
  mutate(season = as.character(season), deploy = as.character(deploy))
predictions_metadata_matched = predictions_metadata %>%
  left_join(
    unit_key %>% select(season, deploy, serialno, unit, unit_agg),
    by = c("season", "deploy", "serialno")
  )

# Drop prediction files that were not matched to a unit
# TODO: Add missing site ID and stage for SMA00403 2020 deployment 2 to "unit_key.csv"
predictions_metadata_no_match = predictions_metadata_matched %>% filter(is.na(unit) | unit == "")
predictions_metadata_no_match = predictions_metadata_no_match %>%
  distinct(season, deploy, serialno)
if (nrow(predictions_metadata_no_match) > 0) {
  message("Prediction files for the following surveys could not be matched to a unit and are discarded:")
  print(as.data.frame(predictions_metadata_no_match))
  predictions_metadata_matched = predictions_metadata_matched %>% filter(!is.na(unit) & unit != "")
}
predictions_metadata_matched = as_tibble(predictions_metadata_matched)

message("Located predictions for ", length(unique(predictions_metadata_matched$unit)), " unique sampling units across ", length(unique(predictions_metadata_matched$season)), " seasons")

# Count the number of prediction files (recordings) per unit-survey pairing
message("Caching prediction file (i.e. recording) counts per unit survey to ", path_out_survey_file_counts)
survey_file_counts = predictions_metadata_matched %>% mutate(survey_date = as.Date(rec_time_start)) %>%
  group_by(season, deploy, serialno, unit, survey_date) %>%
  summarise(n_prediction_files = n_distinct(file_path), .groups = "drop") %>%
  mutate(season = factor(season), unit = factor(unit)) %>%
  arrange(season, unit, survey_date)
if (!dir.exists(path_out_dir)) dir.create(path_out_dir, recursive = TRUE)
write_feather(survey_file_counts, path_out_survey_file_counts)

# Read and aggregate all predictions
message('Aggregating prediction data (this may take some time)')
n = nrow(predictions_metadata_matched)
all_data_list = vector("list", n)  # preallocate list

for (i in seq_len(n)) { # ETA ~45 min
  if (i %% ceiling(n/100) == 0 || i == n) cat(sprintf("%.0f%% (%d/%d)\n", 100 * i / n, i, n))
  file_path = predictions_metadata_matched[[i, 'file_path']]
  rec_time_start = predictions_metadata_matched[[i, 'rec_time_start']]
  # Read the prediction history data
  prediction_data = read_csv(file_path, show_col_types = FALSE) %>% mutate(across(contains("Confidence"), ~ as.numeric(.)))
  if (nrow(prediction_data) > 0) {
    # Compute prediction time, add metadata, and drop irrelevant columns
    prediction_data = prediction_data %>% mutate(
      time     = as.POSIXct(rec_time_start + `Start (s)`),
      season   = predictions_metadata_matched[[i, 'season']],
      deploy   = predictions_metadata_matched[[i, 'deploy']],
      serialno = predictions_metadata_matched[[i, 'serialno']],
      unit     = predictions_metadata_matched[[i, 'unit']],
      unit_agg = predictions_metadata_matched[[i, 'unit_agg']]
    )
  }
  prediction_data = prediction_data %>% select(-`Start (s)`, -`End (s)`, -`Label`, -`Scientific name`)
  all_data_list[[i]] = prediction_data
}
prediction_data = bind_rows(all_data_list)
prediction_data = prediction_data %>% clean_names()

message(nrow(prediction_data), " predictions aggregated")

# Format data
message("Formatting prediction data")
prediction_data = prediction_data %>%
  mutate(across(c(common_name, season, deploy, serialno, unit, unit_agg), as.factor)) %>%
  mutate(common_name = str_to_lower(common_name))

# Write results to cache
write_feather(prediction_data, path_out_prediction_data)
message(crayon::green("Cached prediction data to", path_out_prediction_data))
