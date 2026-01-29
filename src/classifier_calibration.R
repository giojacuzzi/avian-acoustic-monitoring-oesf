## Calculate species-specific probabilistic score thresholds using Platt scaling (Platt 2000, Wood and Kahl 2024).

path_annotations = "/Users/giojacuzzi/Library/CloudStorage/GoogleDrive-giojacuzzi@gmail.com/My Drive/Research/Projects/C4 - OESF avian communities/data/calibration/annotations"
path_annotations_Jacuzzi_and_Olden_2025 = "/Users/giojacuzzi/Library/CloudStorage/GoogleDrive-giojacuzzi@gmail.com/My Drive/Research/Projects/C4 - OESF avian communities/data/calibration/annotations_Jacuzzi_Olden_2025/test_data_annotations.csv"

library(tidyverse)

# Load class labels
class_labels = readLines("models/ensemble/ensemble_class_labels.txt") %>% tolower() %>% tibble(label = .) %>%
  separate(label, into = c("scientific_name", "common_name"), sep = "_", extra = "merge", fill  = "right", remove = FALSE) %>%
  select(label, common_name, scientific_name)

# Load validation annotations from Jacuzzi and Olden 2025 ---------------------------------------------------------------------------------
message("Loading Jacuzzi and Olden 2025 validation annotations from ", path_annotations_Jacuzzi_and_Olden_2025)
annotations_Jacuzzi_Olden_2025 = read_csv(path_annotations_Jacuzzi_and_Olden_2025, show_col_types = FALSE) %>%
  mutate(labels = str_to_lower(labels), focal_class = str_to_lower(focal_class)) %>%
  mutate(not_focal = str_detect(labels, "\\bnot_focal\\b")) %>%
  mutate(unknown = str_detect(labels, "\\bunknown\\b")) %>%
  separate_rows(labels, sep = ", ") %>%
  rename(label = labels)

# Get all true positive files for each class
files_tp = map(class_labels$label, function(lbl) {
  annotations_Jacuzzi_Olden_2025 %>%
    # The file contains the label
    filter(label == lbl) %>%
    pull(file) %>% unique()
})
names(files_tp) = class_labels$label

message("Number of TP samples:")
n_tp = map_int(files_tp, length)
n_tp = tibble(class = names(n_tp), n_tp = n_tp)

# Get all true negative files for each class
# NOTE: Annotations from Jacuzzi and Olden 2025 are exhaustive and pooled across focal species
files_tn = map(class_labels$label, function(lbl) {
  annotations_Jacuzzi_Olden_2025 %>%
    # The file does not contain a true positive
    group_by(file) %>% filter(!any(label == lbl)) %>% ungroup() %>%
    # The file contains "unknown" AND NOT the label, OR the file contains "not_focal" AND focal_class EQUALS the label
    filter( (unknown & label != lbl) | (not_focal & focal_class == lbl) ) %>%
    pull(file) %>% unique()
})
names(files_tn) = class_labels$label

message("Number of TN samples:")
n_tn = map_int(files_tn, length)
n_tn = tibble(class = names(n_tn), n_tn = n_tn)

message("Number of TP and TN samples per class from Jacuzzi and Olden 2025:")
n_tptn = full_join(n_tp, n_tn, by = "class") %>% separate(
  class,
  into = c("scientific_name", "common_name"),
  sep = "_",
  fill = "right",
  extra = "merge",
  remove = FALSE
)
print(n_tptn, n = Inf)

# Load all other validation annotations ---------------------------------------------------------------------------------------------

# Calculate species-specific probabilistic score thresholds using Platt scaling (Platt 2000, Wood and Kahl 2024). Takes as input a parquet file containing the validation dataset following the format below:

# label_predicted       confidence  file                                  label_truth 
# <chr>                 <dbl>       <chr>                                 <chr>       
# townsend's warbler    0.996       101_20210508_090002_1293.0s_1296.0s   townsend's warbler
# varied thrush         0.0003      101_20210508_090002_1293.0s_1296.0s   0            
# violet-green swallow  0.0001      101_20210508_090002_1293.0s_1296.0s   0   

message("TODO: Load all other validation annotations from ", path_validations)

