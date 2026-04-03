# Get predictions above a confidence threshold for the target submodel for a specific class
# These can be used with the audio data to derive species- and model-specific segments

focal_class = "marbled murrelet"
threshold = 0.5

out_dir = paste0("data/debug/helper_get_class_predictions/", focal_class)

pred_data = read_feather("data/cache/1_pam/2_agg_raw_predictions/NEW_prediction_data.feather")
class_data = pred_data %>% filter(common_name == focal_class, confidence_target > threshold) %>% arrange(desc(confidence_target))

class_data = class_data %>% rename(Confidence = confidence_target, `Start (s)` = start_s, `End (s)` = end_s, Label = label, `Scientific name` = scientific_name, `Common name` = common_name) 

class_data$dir  = dirname(class_data$file_path)
class_data$file = sub("\\.BirdNET\\.results\\.csv$", ".wav", basename(class_data$file_path))

# Save .csv files matching directory structure
library(dplyr)
library(readr)
library(purrr)

class_data %>%
  group_split(file_path) %>%
  walk(function(df) {
    
    rel_path <- unique(df$file_path)
    
    # Build full output path under out_dir
    full_path <- file.path(out_dir, rel_path)
    
    # Create nested directories safely
    dir.create(dirname(full_path), recursive = TRUE, showWarnings = FALSE)
    
    # Write CSV
    write_csv(df, full_path)
  })
