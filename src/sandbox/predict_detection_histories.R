# Generate binary detection history predictions from classifier confidence scores using an arbitrary threshold.

# INPUT
# Input data directory structure is assumed to be ".../site/visit/" (e.g. ".../Bp235i/20200618"), such as:
#
# directory
#   site1
#     visit1
#       recording1.csv
#       recording2.csv
#       ...
#     visit2
#     visit3
#   site2
#     visit1
#     visit2
#     ...
#   ...
#
# `site` is a subdirectory that indicates the unique site identifier, and `visit` is a subdirectory of `site`
# that indicates the unique visit identifier (i.e. date) and contains one or more recording .csv files with
# the raw confidence scores for that visit.
path_input_dir = "/Users/giojacuzzi/Desktop/detection_history_demo"
path_species_list = "models/ensemble/ensemble_species_list.txt"
threshold = 0.9

# OUTPUT
path_output_dir = "data/processed/predicted_detection_histories"

############################################################################################################

library(tidyverse)

# Read species list from file
species = sapply(strsplit(readLines(path_species_list), "_"), `[`, 2)
species = tibble(species = species) %>% arrange(species)

# Function to generate binary presence/absence detection history from
# a single .csv of detection confidence scores using a threshold `t`
generate_detection_history_from_predictions = function(path_scores, t) {
  # Read detections from file
  detections = read_csv(path_scores, show_col_types = FALSE)
  # Extract the species common name and confidence scores and rename
  detections = detections[, c('Common name', 'Confidence')]
  # Retain only the most confident detection for each unique species
  detections_max_score <- detections %>%
    group_by(`Common name`) %>%
    summarize(Confidence = max(Confidence, na.rm = TRUE), .groups = "drop")
  # Convert confidence scores to binary presence/absence prediction
  detections_presence <- detections_max_score %>%
    mutate(Confidence = case_when(
      is.na(Confidence) ~ NA_real_,  # preserve NAs
      Confidence >= t ~ 1,
      TRUE ~ 0
    ))
  detections_presence <- detections_presence %>% rename(species = `Common name`, presence = Confidence)
  # Merge with all species not detected
  detection_history_prediction <- species %>%
    left_join(detections_presence, by = "species") %>%
    mutate(presence = coalesce(presence, 0))
  return(detection_history_prediction)
}

# For each site in the input directory, generate a site-specific detection history
paths_site = list.dirs(path_input_dir, recursive = FALSE)
for (path_site in paths_site) {
  site_id = basename(path_site)
  message('Site ', site_id)
  site_detection_history = species
  
  # For each visit associated with the site, generate a visit-specific detection history
  paths_visit = list.dirs(path_site, recursive = FALSE)
  for (path_visit in paths_visit) {
    visit_id = basename(path_visit)
    message('Generating detection history predictions for visit ', visit_id, '...')
    visit_detections = species %>% mutate(presence = 0)
    
    # For each file associated with the visit, generate a file-specific detection history
    files_visit = list.files(path_visit, full.names = TRUE)
    for (file in files_visit) {
      message('File ', basename(file))
      # Generate detection history from predictions in the file
      file_detections = generate_detection_history_from_predictions(file, threshold)
      # Update any species for the visit that was predicted present in the file
      visit_detections = visit_detections %>% mutate(presence = pmax(presence, file_detections$presence == 1))
    }
    message(nrow(visit_detections[visit_detections$presence == 1,]), ' species detected')
    visit_detections = visit_detections %>% rename(!!visit_id := presence)
    # View(visit_detections)
    site_detection_history = site_detection_history %>% full_join(visit_detections, by = "species")
  }
  # View(site_detection_history)
  path_output_file = paste0(path_output_dir, '/', paste0(site_id, '.csv'))
  if (!dir.exists(path_output_dir)) dir.create(path_output_dir, recursive = TRUE)
  write.csv(as.data.frame(site_detection_history), path_output_file, row.names = FALSE)
  message('Saved site detection history to file ', path_output_file)
}
