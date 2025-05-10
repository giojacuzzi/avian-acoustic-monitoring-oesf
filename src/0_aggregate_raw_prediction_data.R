####################################################################################
# Aggregate raw prediction (.csv) data and associated metadata
#
# Directory structure is ".../season/deployment/serialno/date/serialno_date_time.csv"
# For example: "predictions/2020/Deployment1/S4A04271_20200412_Data/S4A04271_20200411/S4A04271_20200411_235938.BirdNET.results.csv"
#
# Assumes that the entire active survey period is represented with individual .csv files (i.e. a recording with no detections should still be represented by an empty .csv file)
#
# Input:
# - Raw .csv prediction data with metadata contained in directory structure (see above)
# - "unit_key.csv" table linking unique sampling unit IDs with season/serialno/deploy combinations
path_in_dir = "/Volumes/gioj/OESF_processed/predictions"
path_unit_key = "data/unit_key.csv"
#
# Output:
# - Cached dataframe of recording counts (i.e. prediction files) per unit-date
# - Cached dataframe of all predictions
path_out_dir = "data/cache/0_aggregate_raw_prediction_data"
####################################################################################

path_out_recording_counts = paste0(path_out_dir, '/recording_counts.feather')
path_out_prediction_data = paste0(path_out_dir, '/prediction_data.feather')

library(tidyverse)
library(janitor)
library(arrow)

# List all individual prediction files
message("Locating all prediction files (i.e. recordings)...")
prediction_files = list.files(path = path_in_dir, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
message("Located ", length(prediction_files), " files")

# Count the number of prediction files (recordings) per serialno-date pairing
message("Caching recording counts...")
serialno_date = sapply(prediction_files, function(f) { basename(dirname(f)) })
recording_counts = as.data.frame(table(serialno_date))
recording_counts$serialno = sapply(recording_counts$serialno_date, function(s) {strsplit(as.character(s), "_")[[1]][1]})
recording_counts$survey_date = sapply(recording_counts$serialno_date, function(s) {strsplit(as.character(s), "_")[[1]][2]})
recording_counts$survey_date = as.Date(recording_counts$survey_date, format = "%Y%m%d", tz = "UTC")
recording_counts = recording_counts %>% select(
    serialno,
    survey_date,
    recordings = Freq
)
if (!dir.exists(path_out_dir)) dir.create(path_out_dir, recursive = TRUE)
write_feather(recording_counts, path_out_recording_counts)

# Read prediction history .csv and parse metadata from path
read_file_and_metadata = function(file_path) {
  
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

  timestamp = as.POSIXct(paste(date,time), format = "%Y%m%d %H%M%S", tz = "UTC")
  
  # Read the prediction history data
  data = read_csv(file_path, show_col_types = FALSE) %>% mutate(across(contains("Confidence"), ~ as.numeric(.)))
  if (nrow(data) > 0) {
    # Add metadata and drop irrelevant columns
    data = data %>% mutate(
        deploy = as.numeric(deploy),
        serialno = serialno,
        time = timestamp + `Start (s)`
    )
    # message(paste(serialno, timestamp))
  }
  data = data %>% select(-`Start (s)`, -`End (s)`)
  return(data)
}

# Read and aggregate all predictions
message('Aggregating prediction data...')
n = length(prediction_files)
all_data_list = vector("list", n)  # preallocate list

for (i in seq_len(n)) {
  if (i %% ceiling(n/100) == 0 || i == n) {
    cat(sprintf("%.0f%% (%d/%d)\n", 100 * i / n, i, n))
  }
  file = prediction_files[i]
  all_data_list[[i]] = read_file_and_metadata(file)
}
prediction_data = bind_rows(all_data_list)
prediction_data = prediction_data %>% clean_names()

# Match season/deploy/serialno combinations to unique sampling unit IDs
message("Matching sampling unit IDs...")
unit_key = read_csv(path_unit_key)
prediction_data_matched = prediction_data %>%
  mutate(season = year(time)) %>%
  left_join(
    unit_key %>% select(season, deploy, serialno, unit, unit_agg),
    by = c("season", "deploy", "serialno")
  )

# TODO: Check for NA unit/unit_agg values

message(nrow(prediction_data_matched), " predictions aggregated for ", length(unique(prediction_data_matched$unit)), " unique sampling units across ", length(unique(prediction_data_matched$season)), " seasons")

# Format data
message("Formatting prediction data...")
prediction_data_matched = prediction_data_matched %>% mutate(across(c(label, common_name, scientific_name, serialno, season, unit, unit_agg), as.factor))

# Write results to cache
message("Caching prediction data...")
write_feather(prediction_data_matched, path_out_prediction_data)
