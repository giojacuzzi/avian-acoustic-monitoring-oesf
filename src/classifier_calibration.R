## Calculate species-specific probabilistic score thresholds using Platt scaling (Platt 2000, Wood and Kahl 2024).

overwrite_prediction_cache = FALSE
overwrite_annotation_cache = FALSE
visualize_calibration_regressions = TRUE

# path_predictions = "/Users/giojacuzzi/Library/CloudStorage/GoogleDrive-giojacuzzi@gmail.com/My Drive/Research/Projects/C4 - OESF avian communities/data/calibration/predictions"
path_jo_predictions_raw = "/Users/giojacuzzi/Library/CloudStorage/GoogleDrive-giojacuzzi@gmail.com/My Drive/Research/Projects/C4 - OESF avian communities/data/calibration/predictions/Jacuzzi_Olden_2025"

path_wadnr_annotations = "/Users/giojacuzzi/Library/CloudStorage/GoogleDrive-giojacuzzi@gmail.com/My Drive/Research/Projects/C4 - OESF avian communities/data/calibration/annotations"
path_jo_annotations_raw = "/Users/giojacuzzi/Library/CloudStorage/GoogleDrive-giojacuzzi@gmail.com/My Drive/Research/Projects/C4 - OESF avian communities/data/calibration/annotations_Jacuzzi_Olden_2025/test_data_annotations.csv"

out_cache_dir = "data/cache/classifier_calibration"
if (!dir.exists(out_cache_dir)) dir.create(out_cache_dir, recursive = TRUE)
path_jo_predictions_cache = paste0(out_cache_dir, "/jo_predictions.rds")
path_jo_annotations_cache = paste0(out_cache_dir, "/jo_annotations.rds")

# Naive thresholds
threshold_min_classifier_score = 0.5 # Naive classifier minimum confidence score threshold to assume binary presence/absence. # "For most false-positive models in our study, using a mid-range threshold of 0.50 or above generally yielded stable estimates." (Katsis et al. 2025)
threshold_min_detected_days = 2 # Minimum number of unique days detected to retain species detections at a site
# Classifier calibration (Platt scaling)
tp_min_prob = 0.95

library(tidyverse)
library(janitor)
library(PRROC)
library(progress)

theme_set(theme_minimal())

conf_to_logit = function(c) {
  c = min(max(c, 0.00001), 1.0 - 0.00001) # prevent undefined logit for extreme scores due to rounding error
  return(log(c / (1 - c)))
}

logit_to_conf = function(l) {
  return(1 / (1 + exp(-l)))
}

# Load class labels
class_labels = readLines("data/models/ensemble/ensemble_class_labels.txt") %>% tolower() %>% tibble(label = .) %>%
  separate(label, into = c("scientific_name", "common_name"), sep = "_", extra = "merge", fill  = "right", remove = FALSE) %>%
  select(label, common_name, scientific_name)

# Aggregate classifier predictions ---------------------------------------------------------------------------------------------

# Aggregate J&O classifier predictions
if (!overwrite_prediction_cache) {
  
  message("Loading cached classifier predictions from ", path_jo_predictions_cache)
  jo_predictions = readRDS(path_jo_predictions_cache)
  
} else {
  
  message("Aggregating raw classifier predictions from ", path_jo_predictions_raw)
  jo_predictions = list.files(path_jo_predictions_raw, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE) %>%
    map_dfr(~ read_csv(.x, show_col_types = FALSE) %>% mutate(file = basename(.x)))
  
  message("Formatting predictions")
  jo_predictions = jo_predictions %>% clean_names() %>%
    mutate(label = str_to_lower(label), common_name = str_to_lower((common_name)), scientific_name = str_to_lower(scientific_name)) %>%
    select(file, label, confidence_source, confidence_target) %>%
    rename(label_predicted = label) %>%
    mutate(file = str_remove(file, "\\.BirdNET\\.results\\.csv$"))
  
  saveRDS(jo_predictions, path_jo_predictions_cache)
  message(crayon::green("Cached", path_jo_predictions_cache))

}

# Aggregate annotations from Jacuzzi and Olden 2025 -------------------------------------------------------------

if (!overwrite_annotation_cache) {

  message("Loading cached JO annotations from ", path_jo_annotations_cache)
  jo_annotations = readRDS(path_jo_annotations_cache)
  
} else {
  
  # Load Jacuzzi and Olden 2025 annotations
  message("Loading Jacuzzi and Olden 2025 annotations from ", path_jo_annotations_raw)
  jo_annotations_raw = read_csv(path_jo_annotations_raw, show_col_types = FALSE) %>%
    mutate(labels = str_to_lower(labels), focal_class = str_to_lower(focal_class)) %>%
    # mutate(not_focal = str_detect(labels, "\\bnot_focal\\b")) %>%
    # mutate(unknown = str_detect(labels, "\\bunknown\\b")) %>%
    mutate(file = str_remove(file, "\\.wav$")) %>%
    separate_rows(labels, sep = ", ") %>%
    rename(label = labels)
  
  # Create a tibble to store annotation data with structure:
  # file,   class_1, class_2, class_3, ...
  # file_1, TP,      unknown, TN,      ...
  # file_2, TN,      TN,      TP,      ...
  # ...
  jo_files = unique(jo_predictions$file)
  jo_annotations = tibble(file = jo_files)
  jo_annotations[class_labels$label] = NA
  
  # Aggregate J&O annotations
  progress = progress_bar$new(format = "[:bar] :percent :elapsedfull (ETA :eta)", total = length(jo_annotations$file), clear = FALSE)
  annotation_warnings = list()
  for (prediction_file in jo_annotations$file) {
    
    # Get all annotations for the file
    file_annotations = jo_annotations_raw %>% filter(file == prediction_file)
    
    # Discard any annotations with missing labels
    if (anyNA(file_annotations$label)) {
      annotation_warnings[["NA"]] = prediction_file
      file_annotations = file_annotations %>% filter(!is.na(label)) 
    }
    
    # Get the focal class for which the file was validated
    focal_class = unique(file_annotations$focal_class)
    
    for (class_label in class_labels$label) {
      
      # TODO: just encode TN as 0 and TP as 1
      a = "TN"
      if (nrow(file_annotations) == 0) {
        # True negative if no labels are present (no annotations for the file exist)
        # NOTE: Annotations from Jacuzzi and Olden 2025 are exhaustive and pooled across focal species
        a = "TN"
      } else if (any(file_annotations$label == class_label)) {
        # True positive if the class label is present
        a = "TP"
      } else if (any(file_annotations$label == "unknown")) {
        # Unknown if "unknown" label is present
        a = "unknown"
      } else if (class_label == focal_class & any(file_annotations$label == "not_focal")) {
        # True negative if "not_focal" is present and the file's focal_class == class_label
        a = "TN"
      } else if (any(file_annotations$label == "not_focal")) {
        # Unknown if "not_focal" is present but the file's focal class != class_label
        # NOTE: Annotations from Jacuzzi and Olden 2025 are exhaustive and pooled across focal species
        a = "unknown"
      }
      # If none of the above cases applied we default to TN
      jo_annotations[jo_annotations$file == prediction_file, class_label] = a
      
    }
    progress$tick()
  }
  
  saveRDS(jo_annotations, path_jo_annotations_cache)
  message(crayon::green("Cached", path_jo_annotations_cache))
  
}

# Count number of TP, TN, unknown, and NA per class
jo_counts = jo_annotations %>%
  pivot_longer(
    cols = all_of(class_labels$label),
    names_to = "class_label",
    values_to = "value",
    values_drop_na = FALSE
  ) %>%
  count(class_label, value, .drop = FALSE) %>%
  pivot_wider(
    names_from = value,
    values_from = n,
    values_fill = 0
  )

# Platt scaling -----------------------------------------------------------------------

# DEBUG
stop("DEBUG")
annotations = jo_annotations

# For each class
for (class_label in class_labels$label) {
  
  message("Processing class ", class_label)
  
  # Get all annotations with matching "class" value
  # Rename "class" to "label_predicted"
  class_annotations = annotations %>% select(file, all_of(class_label))
  
  # Discard files (rows) with "unknown"
  class_annotations = class_annotations %>% filter(.data[[class_label]] != "unknown")
  table(class_annotations[[class_label]])
  
  # TODO just cast as numeric
  # Encode as binary 1 (TP) and 0 (TN)
  class_annotations[[class_label]] = recode(
    class_annotations[[class_label]],
    "TP" = 1, "TN" = 0,
    .default = NA_real_
  )
  
  # Merge with prediction scores
  class_predictions = jo_predictions %>% filter(label_predicted == class_label)
  class_calibration_data = left_join(class_annotations, class_predictions, by = c("file")) %>%
    rename(label_truth = !!sym(class_label))
  
  # Sanity check
  summary(class_calibration_data %>% filter(label_truth == 1) %>% pull(confidence_source))
  summary(class_calibration_data %>% filter(label_truth == 0) %>% pull(confidence_source))
  summary(class_calibration_data %>% filter(label_truth == 1) %>% pull(confidence_target))
  summary(class_calibration_data %>% filter(label_truth == 0) %>% pull(confidence_target))
  
  # Platt scaling
  for (model in c("source", "target")) {
    
    # Convert confidence scores to logit scale
    data_model = class_calibration_data
    data_model$score = sapply(data_model[[paste0("confidence_", model)]], conf_to_logit)
    
    # Calculate PR AUC and ROC AUC of classifier (from raw scores)
    auc_pr <- auc_roc <- NA
    tryCatch(
      {
        auc_pr = pr.curve(scores.class0 = subset(data_model, label_truth == 1) %>% pull(score),
                          scores.class1 = subset(data_model, label_truth == 0) %>% pull(score))$auc.integral
        
        auc_roc = roc.curve(scores.class0 = subset(data_model, label_truth == 1) %>% pull(score),
                            scores.class1 = subset(data_model, label_truth == 0) %>% pull(score))$auc
      },
      warning = function(w) {
        model_warning <<- TRUE
        message(crayon::yellow("WARNING:", conditionMessage(w)))
      },
      error = function(e) {
        model_warning <<- TRUE
        message(crayon::yellow("WARNING:", conditionMessage(e)))
      }
    )
    message("  PR-AUC ", round(auc_pr,3), ", ROC-AUC ", round(auc_roc,3))
    
    # Perform Platt scaling to determine threshold for desired probability of true positive
    threshold = Inf # default threshold is infinite (i.e. do not retain detections unless a species is validated)
    precision_threshold <- recall_threshold <- precision_tmin <- recall_tmin <- NA
    model_warning = FALSE
    tryCatch(
      {
        regression = glm(label_truth ~ score, data_model, family = binomial(link = "logit"))
        intercept   = as.numeric(coef(regression)[1])
        coefficient = as.numeric(coef(regression)[2])
        threshold_logit = (log(tp_min_prob / (1 - tp_min_prob)) - intercept) / coefficient # logit scale
        threshold       = logit_to_conf(threshold_logit) # confidence scale
        
        message("  ", round(threshold, 3), " threshold to achieve Pr(TP)>=", tp_min_prob)
        
        # Calculate estimated precision and recall at this threshold
        # Predicted positive/negative based on threshold
        calc_precision_recall = function(d, t_logit) {
          predicted = ifelse(d$score >= t_logit, 1, 0)
          TP = sum(predicted == 1 & d$label_truth == 1)
          FP = sum(predicted == 1 & d$label_truth == 0)
          FN = sum(predicted == 0 & d$label_truth == 1)
          return(data.frame(
            precision = TP / (TP + FP),
            recall    = TP / (TP + FN)
          ))
        }
        perf_t = calc_precision_recall(data_model, threshold_logit)
        precision_threshold  = perf_t$precision
        recall_threshold     = perf_t$recall
        message("  Performance at threshold ", round(threshold,3), ":\n  Precision ", round(precision_threshold,3), "\n  Recall ", round(recall_threshold,3))
        
        # Calculate at minimum threshold
        perf_tmin = calc_precision_recall(data_model, conf_to_logit(threshold_min_classifier_score))
        precision_tmin = perf_tmin$precision
        recall_tmin    = perf_tmin$recall
        message("  Performance at minimum threshold ", threshold_min_classifier_score, ":\n  Precision ", round(precision_tmin,3), "\n  Recall ", round(recall_tmin,3))
        
        if (visualize_calibration_regressions) {
          
          ## Visualize precision and recall performance as a function of score
          data_sorted = data_model %>% mutate(score = logit_to_conf(score)) %>% arrange(desc(score))
          n_pos = sum(data_sorted$label_truth == 1)
          data_sorted = data_sorted %>% mutate(
            tp = cumsum(label_truth == 1),
            fp = cumsum(label_truth == 0),
            recall = tp / n_pos,
            precision = ifelse(tp + fp == 0, 1, tp / (tp + fp))  # handle division by zero
          )
          plt = ggplot(data_sorted, aes(x = score)) +
            geom_line(aes(y = precision, color = "Precision"), linewidth = 1) +
            geom_line(aes(y = recall, color = "Recall"), linewidth = 1) +
            labs(title = paste0(class_label), x = "Score", y = "Performance", color = "Metric"); print(plt)
          
          ## Visualize the logistic regression and data
          x_range = seq(min(data_model$score), max(data_model$score), length.out = 100)
          
          # Predict fitted values and calculate confidence intervals on probability scale
          pred = predict(regression, newdata = data.frame(score = x_range), type = "link", se.fit = TRUE)
          
          # Calculate confidence intervals on probability scale
          z = qnorm(0.975)  # 95% CI
          upper  = logit_to_conf(pred$fit + z * pred$se.fit)
          lower  = logit_to_conf(pred$fit - z * pred$se.fit)
          y_pred = logit_to_conf(pred$fit)
          
          regression_df = data.frame(score = x_range, prob = y_pred, lower = lower, upper = upper, warning = model_warning)
          plt = ggplot() +
            geom_vline(xintercept = 0, color = "gray", linetype = "solid", linewidth = 0.5) +
            geom_point(data = data_model, aes(x = score, y = label_truth), shape = 1, alpha = 0.5) +
            geom_hline(yintercept = tp_min_prob, color = "black", linetype = "dashed", linewidth = 0.5) +
            geom_vline(xintercept = threshold_logit, color = "black", linetype = "dashed", linewidth = 0.5) +
            geom_ribbon(data = regression_df, aes(x = score, ymin = lower, ymax = upper, fill = warning), alpha = 0.2) +
            geom_line(data = regression_df, aes(x = score, y = prob, color = warning), linewidth = 1) +
            scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red")) +
            scale_fill_manual(values = c("FALSE" = "blue", "TRUE" = "red")) +
            labs(
              x = "Score (logit)", y = "True positive probability",
              title = class_label, subtitle = paste0("model '", model, "', threshold p(TP) ≥ ", tp_min_prob, " = ", round(threshold, 3))
            ); print(plt)
        }
      },
      warning = function(w) {
        model_warning <<- TRUE
        message(crayon::yellow("WARNING", conditionMessage(w)))
      },
      error = function(e) {
        stop(crayon::red("ERROR", conditionMessage(e)))
      }
    )
    
  }
}

################################################################################################################################################





# Load validation annotations from Jacuzzi and Olden 2025 ---------------------------------------------------------------------------------
message("Loading Jacuzzi and Olden 2025 annotations from ", path_jo_annotations)
jo_annotations = read_csv(path_jo_annotations, show_col_types = FALSE) %>%
  mutate(labels = str_to_lower(labels), focal_class = str_to_lower(focal_class)) %>%
  mutate(not_focal = str_detect(labels, "\\bnot_focal\\b")) %>%
  mutate(unknown = str_detect(labels, "\\bunknown\\b")) %>%
  separate_rows(labels, sep = ", ") %>%
  rename(label = labels)

jo_files_tp = list()
# Loop over each class
for (lbl in class_labels$label) {
  jo_files_tp[[lbl]] = jo_annotations %>% filter(label == lbl) %>% pull(file)
}
jo_files_tn = list()
for (lbl in class_labels$label) {
  lbl_annotations = jo_annotations %>% filter(!file %in% jo_files_tp[[lbl]])
  lbl_annotations = lbl_annotations %>% filter(unknown == FALSE)
  # TODO: Unannotated files (e.g. Background files) are intepreted as having no relevant signals (i.e. labels) present
  lbl_annotations = lbl_annotations %>% filter(focal_class == lbl)
  lbl_nf = lbl_annotations %>% filter()
  jo_files_tn[[lbl]] = lbl_annotations %>% filter(focal_class == lbl & FALSE) %>% pull(file)
}

stop("FIX ANNOTATIONS ABOVE")

# message("Finding all true positive files for each class")
# jo_files_tp = map(class_labels$label, function(lbl) {
#   jo_annotations %>%
#     # The file contains the label
#     filter(label == lbl) %>%
#     pull(file) %>% unique() %>% str_remove("\\.wav$")
# })
# names(jo_files_tp) = class_labels$label
# 
# jo_n_tp = map_int(jo_files_tp, length)
# jo_n_tp = tibble(class = names(jo_n_tp), n_tp = jo_n_tp)
# 
# # NOTE: Annotations from Jacuzzi and Olden 2025 are exhaustive and pooled across focal species
# message("Finding all true negative files for each class")
# jo_files_tn = map(class_labels$label, function(lbl) {
#   jo_annotations %>%
#     # The file does not contain a true positive
#     group_by(file) %>% filter(!any(label == lbl)) %>% ungroup() %>%
#     # The file does not contain an "unknown"
#     group_by(file) %>% filter(!unknown) %>% ungroup() %>%
#     # TODO: Also filter by not_focal?
#     # "Segments from which the absence of the predicted species could be definitively determined, but the identity of other species could not be, were labeled as “not target” and only included in performance evaluations for the predicted species class."
#     pull(file) %>% unique() %>% str_remove("\\.wav$")
# })
# names(jo_files_tn) = class_labels$label

jo_n_tn = map_int(jo_files_tn, length)
jo_n_tn = tibble(class = names(jo_n_tn), n_tn = jo_n_tn)

message("Number of TP and TN samples per class from Jacuzzi and Olden 2025:")
jo_n_tptn = full_join(jo_n_tp, jo_n_tn, by = "class") %>% separate(
  class,
  into = c("scientific_name", "common_name"),
  sep = "_",
  fill = "right",
  extra = "merge",
  remove = FALSE
)
print(jo_n_tptn, n = Inf)

# Load additional validation annotations ---------------------------------------------------------------------------------

# NOTE: These annotations are NOT exhaustive and should not be pooled across species. Rather, they are species-specific, reflecting the presence or absence of one species only.

message("Loading WADNR annotations from ", path_wadnr_annotations)
annotations_wadnr = list.files(path_wadnr_annotations, pattern = "\\.txt$", recursive = TRUE, full.names = TRUE) %>%
  map_dfr(~ read_tsv(.x, show_col_types = FALSE) %>%
  mutate(source_dir = str_to_lower(basename(dirname(.x))))) %>%
  clean_names() %>% rename(file = begin_file, focal_class = source_dir) %>% select(file, focal_class, label) %>%
  mutate(label = str_to_lower(label)) %>% mutate(label = if_else(label %in% c("non_target", "not_target"), "not_focal", label))

annotations_wadnr = annotations_wadnr %>%
  group_by(file) %>% mutate(unknown = any(label == "unknown")) %>% ungroup()
annotations_wadnr = annotations_wadnr %>%
  group_by(file) %>% mutate(not_focal = any(label == "not_focal")) %>% ungroup()

message("Finding all true positive files for each class")
files_tp = map(class_labels$common_name, function(lbl) {
  annotations_wadnr %>%
    # The file contains the label
    filter(label == lbl) %>%
    pull(file) %>% unique() %>% str_remove("\\.wav$")
})
names(files_tp) = class_labels$label

n_tp = map_int(files_tp, length)
n_tp = tibble(class = names(n_tp), n_tp = n_tp)

message("Finding all true negative files for each class")
files_tn = map(class_labels$common_name, function(lbl) {
  annotations_wadnr %>%
    # The file does not contain a true positive
    group_by(file) %>% filter(!any(label == lbl)) %>% ungroup() %>%
    # The file does not contain an unknown
    group_by(file) %>% filter(!unknown) %>% ungroup() %>%
    # The file contains file contains "not_focal" AND focal_class EQUALS the label
    filter((not_focal & focal_class == lbl) ) %>%
    pull(file) %>% unique() %>% str_remove("\\.wav$")
})
names(files_tn) = class_labels$label

files_tp[["bonasa umbellus_ruffed grouse"]]
files_tn[["bonasa umbellus_ruffed grouse"]]

n_tn = map_int(files_tn, length)
n_tn = tibble(class = names(n_tn), n_tn = n_tn)

message("Number of TP and TN samples per class from additional annotations:")
n_tptn = full_join(n_tp, n_tn, by = "class") %>% separate(
  class,
  into = c("scientific_name", "common_name"),
  sep = "_",
  fill = "right",
  extra = "merge",
  remove = FALSE
)
print(n_tptn, n = Inf)

# Join all validation annotations ---------------------------------------------------------------------------------

stopifnot(all(
  jo_n_tptn$class == n_tptn$class &
  jo_n_tptn$scientific_name == n_tptn$scientific_name &
  jo_n_tptn$common_name == n_tptn$common_name
))

message("Number of TP and TN samples per class from all annotations:")
total_tptn = jo_n_tptn
total_tptn$n_tp = total_tptn$n_tp + n_tptn$n_tp
total_tptn$n_tn = total_tptn$n_tn + n_tptn$n_tn
print(total_tptn, n = Inf)

# message("Classes with < 100 TP:")
# print(total_tptn %>% filter(n_tp < 100), n = Inf)


# DEBUG
# all_files_tp = Map(c, files_tp, jo_files_tp)
# all_files_tn = Map(c, files_tn, jo_files_tn)

# Just use WADNR data for now
all_files_tp = files_tp
all_files_tn = files_tn
# DEBUG

stopifnot(all(names(jo_files_tp) == names(files_tp)))
stopifnot(all(names(jo_files_tn) == names(files_tn)))

all_files_tp_df = lapply(names(all_files_tp), function(l) {
  n = length(all_files_tp[[l]])
  data.frame(file = all_files_tp[[l]], class = rep(l, n), label_truth = rep(l, n), stringsAsFactors = FALSE)
}) %>% bind_rows()

all_files_tn_df = lapply(names(all_files_tn), function(l) {
  n = length(all_files_tn[[l]])
  data.frame(
    file = all_files_tn[[l]], class = rep(l, n), label_truth = rep("TN", n), stringsAsFactors = FALSE)
}) %>% bind_rows()

annotations = bind_rows(all_files_tp_df, all_files_tn_df)

# Calculate species-specific probabilistic score thresholds using Platt scaling (Platt 2000, Wood and Kahl 2024). Takes as input a parquet file containing the validation dataset following the format below:

stop("DEBUG")

# For each class
for (class_label in class_labels$label) {
  
  message("Processing class ", class_label)

  # Get all annotations with matching "class" value
  # Rename "class" to "label_predicted"
  class_annotations = annotations %>% filter(class == class_label) %>% rename(label_predicted = class) %>% arrange(file)
  table(class_annotations$label_truth)

  # Get all predictions with matching "label_predicted" value that we have annotations for
  class_predictions = predictions %>% filter(label_predicted == class_label, file %in% class_annotations$file) %>% arrange(file)

  # Merge into the following format:
  # label_predicted       confidence  file                                  label_truth 
  # <chr>                 <dbl>       <chr>                                 <chr>       
  # townsend's warbler    0.996       101_20210508_090002_1293.0s_1296.0s   townsend's warbler
  # varied thrush         0.0003      101_20210508_090002_1293.0s_1296.0s   0            
  # violet-green swallow  0.0001      101_20210508_090002_1293.0s_1296.0s   0   
  data = class_annotations %>% right_join(class_predictions, by = c("label_predicted", "file"))
  table(data$label_truth)
  
  # Create binary indicator for each unique audio file for true presence (aggregating all validations for each file)
  data = data %>% mutate(label_truth = as.integer(label_truth == class_label))
  stopifnot(all(unique(data$label_truth) %in% c(0, 1)))
  
  # Sanity check
  summary(data %>% filter(label_truth == 1) %>% pull(confidence_source))
  summary(data %>% filter(label_truth == 0) %>% pull(confidence_source))
  
  for (model in c("source", "target")) {
    
    # Convert confidence scores to logit scale
    data_model = data
    data_model$score = sapply(data_model[[paste0("confidence_", model)]], conf_to_logit)
    
    # Calculate PR AUC and ROC AUC of classifier (from raw scores)
    auc_pr = pr.curve(scores.class0 = subset(data_model, label_truth == 1) %>% pull(score),
                      scores.class1 = subset(data_model, label_truth == 0) %>% pull(score))$auc.integral
    
    auc_roc = roc.curve(scores.class0 = subset(data_model, label_truth == 1) %>% pull(score),
                        scores.class1 = subset(data_model, label_truth == 0) %>% pull(score))$auc
    
    message("  PR-AUC ", round(auc_pr,3), ", ROC-AUC ", round(auc_roc,3))
    
    # Perform Platt scaling to determine threshold for desired probability of true positive
    threshold = Inf # default threshold is infinite (i.e. do not retain detections unless a species is validated)
    precision_threshold <- recall_threshold <- precision_tmin <- recall_tmin <- NA
    model_warning = FALSE
    tryCatch(
      {
        regression = glm(label_truth ~ score, data_model, family = binomial(link = "logit"))
        intercept   = as.numeric(coef(regression)[1])
        coefficient = as.numeric(coef(regression)[2])
        threshold_logit = (log(tp_min_prob / (1 - tp_min_prob)) - intercept) / coefficient # logit scale
        threshold       = logit_to_conf(threshold_logit) # confidence scale
        
        message("  ", round(threshold, 3), " threshold to achieve Pr(TP)>=", tp_min_prob)
        
        # Calculate estimated precision and recall at this threshold
        # Predicted positive/negative based on threshold
        calc_precision_recall = function(d, t_logit) {
          predicted = ifelse(d$score >= t_logit, 1, 0)
          TP = sum(predicted == 1 & d$label_truth == 1)
          FP = sum(predicted == 1 & d$label_truth == 0)
          FN = sum(predicted == 0 & d$label_truth == 1)
          return(data.frame(
            precision = TP / (TP + FP),
            recall    = TP / (TP + FN)
          ))
        }
        perf_t = calc_precision_recall(data_model, threshold_logit)
        precision_threshold  = perf_t$precision
        recall_threshold     = perf_t$recall
        message("  Performance at threshold ", round(threshold,3), ":\n  Precision ", round(precision_threshold,3), "\n  Recall ", round(recall_threshold,3))
        
        # Calculate at minimum threshold
        perf_tmin = calc_precision_recall(data_model, conf_to_logit(threshold_min_classifier_score))
        precision_tmin = perf_tmin$precision
        recall_tmin    = perf_tmin$recall
        message("  Performance at minimum threshold ", threshold_min_classifier_score, ":\n  Precision ", round(precision_tmin,3), "\n  Recall ", round(recall_tmin,3))
        
        if (TRUE) {
          
          ## Visualize precision and recall performance as a function of score
          data_sorted = data_model %>% mutate(score = logit_to_conf(score)) %>% arrange(desc(score))
          n_pos = sum(data_sorted$label_truth == 1)
          data_sorted = data_sorted %>% mutate(
            tp = cumsum(label_truth == 1),
            fp = cumsum(label_truth == 0),
            recall = tp / n_pos,
            precision = ifelse(tp + fp == 0, 1, tp / (tp + fp))  # handle division by zero
          )
          plt = ggplot(data_sorted, aes(x = score)) +
            geom_line(aes(y = precision, color = "Precision"), linewidth = 1) +
            geom_line(aes(y = recall, color = "Recall"), linewidth = 1) +
            labs(title = paste0(class_label), x = "Score", y = "Performance", color = "Metric"); print(plt)
          
          ## Visualize the logistic regression and data
          x_range = seq(min(data_model$score), max(data_model$score), length.out = 100)
          
          # Predict fitted values and calculate confidence intervals on probability scale
          pred = predict(regression, newdata = data.frame(score = x_range), type = "link", se.fit = TRUE)
          
          # Calculate confidence intervals on probability scale
          z = qnorm(0.975)  # 95% CI
          upper  = logit_to_conf(pred$fit + z * pred$se.fit)
          lower  = logit_to_conf(pred$fit - z * pred$se.fit)
          y_pred = logit_to_conf(pred$fit)
          
          regression_df = data.frame(score = x_range, prob = y_pred, lower = lower, upper = upper, warning = model_warning)
          plt = ggplot() +
            geom_vline(xintercept = 0, color = "gray", linetype = "solid", linewidth = 0.5) +
            geom_point(data = data_model, aes(x = score, y = label_truth), shape = 1, alpha = 0.5) +
            geom_hline(yintercept = tp_min_prob, color = "black", linetype = "dashed", linewidth = 0.5) +
            geom_vline(xintercept = threshold_logit, color = "black", linetype = "dashed", linewidth = 0.5) +
            geom_ribbon(data = regression_df, aes(x = score, ymin = lower, ymax = upper, fill = warning), alpha = 0.2) +
            geom_line(data = regression_df, aes(x = score, y = prob, color = warning), linewidth = 1) +
            scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red")) +
            scale_fill_manual(values = c("FALSE" = "blue", "TRUE" = "red")) +
            labs(
              x = "Score (logit)", y = "True positive probability",
              title = paste0(class_label), subtitle = paste0("Threshold p(TP) ≥ ", tp_min_prob, " = ", round(threshold, 3))
            ); print(plt)
        }
      },
      warning = function(w) {
        model_warning <<- TRUE
        message(crayon::yellow("WARNING", conditionMessage(w)))
      },
      error = function(e) {
        stop(crayon::red("ERROR", conditionMessage(e)))
      }
    )
    
  }
}

