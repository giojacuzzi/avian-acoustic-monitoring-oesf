## 1c_classifier_calibration.R #########################################################################################
# Calculate species-specific probabilistic score thresholds using Platt scaling (Platt 2000, Wood and Kahl 2024).
# 
## CONFIG:
calibration_class = "all" # class to calibrate (e.g. "all" or "catharus ustulatus_swainson's thrush")
overwrite_annotation_cache = FALSE   # Aggregate annotations from the raw data
overwrite_prediction_cache = FALSE   # Aggregate predictions from the raw data
threshold_min_classifier_score = 0.5 # Naive classifier minimum confidence score threshold to assume binary presence/absence. "For most false-positive models in our study, using a mid-range threshold of 0.50 or above generally yielded stable estimates." (Katsis et al. 2025)
threshold_min_detected_days = 2      # Minimum number of unique days detected to retain species detections at a site
tp_min_prob = 0.95                   # Desired minimum probability of a true positive
#
## OUTPUT: calibration_results.csv
out_cache_dir = "data/cache/1_pam/1_classifier_calibration"
path_annotations_clean_cache         = paste0(out_cache_dir, "/annotations_clean.csv")
path_calibration_results_raw_cache   = paste0(out_cache_dir, "/calibration_results_raw.csv")   # raw optimized calibration results for subsequent modeling
path_calibration_results_clean_cache = paste0(out_cache_dir, "/calibration_results_clean.csv") # clean optimized calibration results for inspection and print
# path_calibration_data_cache    = paste0(out_cache_dir, "/calibration_data.rds")    # raw calibration data and plots
#
## INPUT:
path_site_key = "data/sites/site_key_long.csv"
# Directories containing raw .csv predictions for calibration, e.g.
# confidence_source	confidence_target	label
# 0.566		          0.678	            Loxia curvirostra_Red Crossbill
# 0.843		          0.998	            Melospiza melodia_Song Sparrow
# ...       ...     ...
path_jo_predictions_raw    = "/Users/giojacuzzi/Library/CloudStorage/GoogleDrive-giojacuzzi@gmail.com/My Drive/Research/Projects/C4 - OESF avian communities/data/calibration/predictions/Jacuzzi_Olden_2025"
path_wadnr_predictions_raw = "/Users/giojacuzzi/Library/CloudStorage/GoogleDrive-giojacuzzi@gmail.com/My Drive/Research/Projects/C4 - OESF avian communities/data/calibration/predictions/WADNR"
#
# Directories containing raw annotations for calibration
# e.g. for WADNR .selections.txt Raven Pro files:
# Selection	Begin File                                                  Label
# 1         105_0.594_290_SMA00424_20200527_060059_1941.0s_1944.0s.wav	Brown Creeper
# 2         114_0.513_176_SMA00403_20200605_060010_642.0s_645.0s.wav  	unknown
# ...       ...                                                         ...
#
# e.g. for Jacuzzi and Olden 2025 .csv summary:
# focal_class	                      file	                                                        labels
# Bombycilla cedrorum_Cedar Waxwing	Cedar Waxwing-0.9761422611646536_SMA00346_20200725_145917.wav	"Bombycilla cedrorum_Cedar Waxwing, not_focal, Ixoreus naevius_Varied Thrush"
# Bombycilla cedrorum_Cedar Waxwing	Cedar Waxwing-0.101301789380809_SMA00370_20200716_060341.wav	not_focal, unknown
# ...                               ...                                                           ...
path_jo_annotations_raw    = "/Users/giojacuzzi/Library/CloudStorage/GoogleDrive-giojacuzzi@gmail.com/My Drive/Research/Projects/C4 - OESF avian communities/data/calibration/annotations/Jacuzzi_Olden_2025/test_data_annotations.csv"
path_wadnr_annotations_raw = "/Users/giojacuzzi/Library/CloudStorage/GoogleDrive-giojacuzzi@gmail.com/My Drive/Research/Projects/C4 - OESF avian communities/data/calibration/annotations/WADNR"
#######################################################################################################################

source("src/global.R")

# Helper variables and functions ======================================================================================

if (!dir.exists(out_cache_dir)) dir.create(out_cache_dir, recursive = TRUE)
path_jo_predictions_cache = paste0(out_cache_dir, "/jo_predictions.rds")
path_jo_annotations_cache = paste0(out_cache_dir, "/jo_annotations.rds")
path_wadnr_predictions_cache = paste0(out_cache_dir, "/wadnr_predictions.rds")
path_wadnr_annotations_cache = paste0(out_cache_dir, "/wadnr_annotations.rds")

conf_to_logit = function(c) {
  c = min(max(c, 0.00001), 1.0 - 0.00001) # prevent undefined logit for extreme scores due to rounding error
  return(log(c / (1 - c)))
}

logit_to_conf = function(l) {
  return(1 / (1 + exp(-l)))
}

calibration_class = str_to_lower(calibration_class)
if (calibration_class == "all") calibration_class = class_labels$label

# Aggregate classifier predictions ---------------------------------------------------------------------------------------------

# Aggregate J&O classifier predictions
if (!overwrite_prediction_cache) {
  message("Loading cached Jacuzzi and Olden (2025) predictions from ", path_jo_predictions_cache)
  jo_predictions = readRDS(path_jo_predictions_cache)
  
} else {
  message("Aggregating raw Jacuzzi and Olden (2025) classifier predictions from ", path_jo_predictions_raw) # ETA 2 min
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

# Aggregate WADNR classifier predictions
if (!overwrite_prediction_cache) {
  message("Loading cached WADNR predictions from ", path_wadnr_predictions_cache)
  wadnr_predictions = readRDS(path_wadnr_predictions_cache)
  
} else {
  message("Aggregating raw WADNR classifier predictions from ", path_wadnr_predictions_raw) # ETA 15 min
  wadnr_predictions = list.files(path_wadnr_predictions_raw, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE) %>%
    map_dfr(~ read_csv(.x, show_col_types = FALSE) %>% mutate(file = basename(.x)))
  
  message("Formatting predictions") # ETA 2 min
  wadnr_predictions = wadnr_predictions %>% clean_names() %>%
    mutate(label = str_to_lower(label), common_name = str_to_lower((common_name)), scientific_name = str_to_lower(scientific_name)) %>%
    select(file, label, confidence_source, confidence_target) %>%
    rename(label_predicted = label) %>%
    mutate(file = str_remove(file, "\\.BirdNET\\.results\\.csv$"))
  
  saveRDS(wadnr_predictions, path_wadnr_predictions_cache)
  message(crayon::green("Cached", path_wadnr_predictions_cache))
}

# Aggregate annotations ------------------------------------------------------------------------------------------

# Aggregate J&O classifier annotations
if (!overwrite_annotation_cache) {
  message("Loading cached Jacuzzi and Olden (2025) annotations from ", path_jo_annotations_cache)
  jo_annotations = readRDS(path_jo_annotations_cache)
  
} else {
  message("Loading raw Jacuzzi and Olden (2025) annotations from ", path_jo_annotations_raw)
  jo_annotations_raw = read_csv(path_jo_annotations_raw, show_col_types = FALSE) %>%
    mutate(labels = str_to_lower(labels), focal_class = str_to_lower(focal_class)) %>%
    mutate(file = str_remove(file, "\\.wav$")) %>%
    separate_rows(labels, sep = ", ") %>%
    rename(label = labels)
  
  # Create a tibble to store annotation data with structure:
  # file,   class_1, class_2, class_3, ...
  # file_1, TP,      unknown, TN,      ...
  # file_2, TN,      TN,      TP,      ...
  # ...
  jo_files = unique(jo_predictions$file) # Consider all files, whether annotations exist or not (no annotations == no species presence)
  jo_annotations = tibble(file = jo_files)
  jo_annotations[calibration_class] = NA
  
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
    
    for (class_label in calibration_class) {
      
      a = "0"
      if (nrow(file_annotations) == 0) {
        # True negative if no labels are present (no annotations for the file exist)
        # NOTE: Annotations from Jacuzzi and Olden 2025 are exhaustive and pooled across focal species
        a = "0"
      } else if (any(file_annotations$label == class_label)) {
        # True positive if the class label is present
        a = "1"
      } else if (any(file_annotations$label == "unknown")) {
        # Unknown if "unknown" label is present
        a = "unknown"
      } else if (class_label == focal_class & any(file_annotations$label == "not_focal")) {
        # True negative if "not_focal" is present and the file's focal_class == class_label
        a = "0"
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
    cols = all_of(calibration_class),
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

message("J&0 annotation counts:")
print(jo_counts, n = Inf)

# Aggregate WADNR classifier annotations
if (!overwrite_annotation_cache) {
  message("Loading cached WADNR annotations from ", path_wadnr_annotations_cache)
  wadnr_annotations = readRDS(path_wadnr_annotations_cache)

} else {
  message("Loading raw WADNR annotations from ", path_wadnr_annotations_raw)
  wadnr_annotations_raw = list.files(path_wadnr_annotations_raw, pattern = "\\.txt$", recursive = TRUE, full.names = TRUE) %>%
    map_dfr(~ read_tsv(.x, show_col_types = FALSE) %>%
              mutate(source_dir = str_to_lower(basename(dirname(.x))))) %>%
    clean_names() %>% rename(file = begin_file, focal_class = source_dir) %>% select(file, focal_class, label) %>%
    mutate(file = str_remove(file, "\\.wav$")) %>%
    mutate(label = str_to_lower(label)) %>% mutate(label = if_else(label %in% c("non_target", "not_target"), "not_focal", label))
  
  # Standardize labels for focal_class and label
  wadnr_annotations_raw = wadnr_annotations_raw %>%
    left_join(
      class_labels %>% select(common_name, label),
      by = c("label" = "common_name"),
      suffix = c("", "_label")
    ) %>%
    left_join(
      class_labels %>% select(common_name, label),
      by = c("focal_class" = "common_name"),
      suffix = c("", "_focal")
    ) %>%
    mutate(
      label = coalesce(label_label, label),
      focal_class = coalesce(label_focal, focal_class)
    ) %>%
    select(-label_label, -label_focal)

  wadnr_files = unique(wadnr_annotations_raw$file) # Only consider annotated files
  wadnr_annotations = tibble(file = wadnr_files)
  wadnr_annotations[calibration_class] = NA
  
  progress = progress_bar$new(format = "[:bar] :percent :elapsedfull (ETA :eta)", total = length(wadnr_annotations$file), clear = FALSE)
  annotation_warnings = list()
  for (prediction_file in wadnr_annotations$file) {
    
    # Get all annotations for the file
    file_annotations = wadnr_annotations_raw %>% filter(file == prediction_file)
    
    # Discard any annotations with missing labels
    if (anyNA(file_annotations$label)) {
      annotation_warnings[["NA"]] = prediction_file
      file_annotations = file_annotations %>% filter(!is.na(label))
    }
    
    # Get the focal class(es) for which the file was validated
    focal_class = unique(file_annotations$focal_class)
    
    for (class_label in calibration_class) {
      
      a = "unknown"
      if (nrow(file_annotations) == 0) {
        # Unknown no labels are present (no annotations for the file exist)
        a = "unknown"
      } else if (any(file_annotations$label == class_label)) {
        # True positive if the class label is present
        a = "1"
      } else if (any(file_annotations$label == "unknown")) {
        # Unknown if "unknown" label is present
        a = "unknown"
      } else if (any(file_annotations$label == "not_focal") & (class_label %in% (file_annotations %>% filter(focal_class == class_label) %>% pull(focal_class)))) {
        # True negative if "not_focal" is present and the file's focal_class == class_label
        a = "0"
      } else if (any(file_annotations$label == "not_focal")) {
        # Unknown if "not_focal" is present but the file's focal class != class_label
        a = "unknown"
      }
      # If none of the above cases applied we default to unknown
      wadnr_annotations[wadnr_annotations$file == prediction_file, class_label] = a
      
    }
    progress$tick()
  }
  saveRDS(wadnr_annotations, path_wadnr_annotations_cache)
  message(crayon::green("Cached", path_wadnr_annotations_cache))
}

# Count number of TP, TN, unknown, and NA per class
wadnr_counts = wadnr_annotations %>%
  pivot_longer(
    cols = all_of(calibration_class),
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

message("WADNR annotation counts:")
print(wadnr_counts, n = Inf)

# Join all calibration data ---------------------------------------------------------------------------------

if (!"1" %in% names(jo_counts)) jo_counts$`1` = 0
if (!"0" %in% names(jo_counts)) jo_counts$`0` = 0
if (!"1" %in% names(wadnr_counts)) wadnr_counts$`1` = 0
if (!"0" %in% names(wadnr_counts)) wadnr_counts$`0` = 0

stopifnot(all(jo_counts$class_label == wadnr_counts$class_label))
stopifnot(all(sort(calibration_class) == jo_counts$class_label))

message("Number of TP and TN samples per class from all annotations:")
total_counts = tibble(class_label = sort(calibration_class))
total_counts$`1` = jo_counts$`1` + wadnr_counts$`1`
total_counts$`0` = jo_counts$`0` + wadnr_counts$`0`
total_counts$`unknown` = jo_counts$`unknown` + wadnr_counts$`unknown`
total_counts = total_counts %>% full_join(class_labels %>% rename(class_label = label), by = "class_label")
print(total_counts %>% select(common_name, `1`, `0`, `unknown`) %>% arrange(common_name), n = Inf)

stopifnot(colnames(jo_annotations) == colnames(wadnr_annotations))
predictions = bind_rows(jo_predictions, wadnr_predictions)
annotations = bind_rows(jo_annotations, wadnr_annotations)

if (overwrite_annotation_cache) {
  annotations_clean = annotations %>%
    rename_with(
      ~ ifelse(.x == "file", .x, str_replace(.x, ".*?_", "")) 
    ) %>%
    mutate(
      # Extract serial number with either SMA or S4A
      serialno = str_extract(file, "S(MA|4A)\\d+"),
      # Extract season year (YYYY) if present after the serial
      season = str_extract(file, "S(MA|4A)\\d+_(\\d{4})") %>%  str_remove("^S(MA|4A)\\d+_"),
      # Extract full date (YYYYMMDD) if present after the serial
      yday = str_extract(file, "S(MA|4A)\\d+_(\\d{8})") %>% str_remove("^S(MA|4A)\\d+_") %>%  as.Date(format = "%Y%m%d") %>% yday()
    )
  site_key = read.csv(path_site_key) %>% mutate(site = tolower(site), site_agg = tolower(site_agg), date = as.Date(date, format = "%m/%d/%y"), season = as.character(year(date)))
  site_key$yday = yday(site_key$date)
  annotations_clean = annotations_clean %>%
    left_join(
      site_key %>% select(serialno, yday, site, site_agg, season),
      by = c("serialno", "yday", "season")
    )
  if (anyNA(annotations_clean$site)) {
    message(yellow("WARNING: Missing site ID match for the following annotation files:"))
    annotations_clean %>% filter(is.na(site)) %>% select(file, serialno, season, yday, site) %>% print(n = Inf)
  }
  
  write.csv(annotations_clean, path_annotations_clean_cache)
  message(crayon::green("Cached", path_annotations_clean_cache))
}

# Platt scaling -----------------------------------------------------------------------

# Calculate species-specific probabilistic score thresholds using Platt scaling (Platt 2000, Wood and Kahl 2024). Takes as input a parquet file containing the validation dataset following the format below:
calibrate = function(preds, anno, labels) {
  calibration_data = list()
  plots = list()
  stats = tibble()
  for (class_label in labels) {
    
    message(class_label)
    
    # Get all annotations with matching "class" value
    # Rename "class" to "label_predicted"
    class_annotations = anno %>% select(file, all_of(class_label)) %>%
      rename(label_truth = !!sym(class_label))
    
    # Discard files (rows) with "unknown" and cast to integer
    class_annotations = class_annotations %>% filter(label_truth != "unknown") %>%
      mutate(label_truth = as.integer(label_truth))
    table(class_annotations$label_truth)
    
    # Merge with prediction scores
    class_predictions = preds %>% filter(label_predicted == class_label)
    class_calibration_data = left_join(class_annotations, class_predictions, by = c("file"))
    table(class_calibration_data$label_truth)
    
    # Debugging sanity check
    # summary(class_calibration_data %>% filter(label_truth == 1) %>% pull(confidence_source))
    # summary(class_calibration_data %>% filter(label_truth == 0) %>% pull(confidence_source))
    # summary(class_calibration_data %>% filter(label_truth == 1) %>% pull(confidence_target))
    # summary(class_calibration_data %>% filter(label_truth == 0) %>% pull(confidence_target))
    
    # Platt scaling
    plot_pr = list()
    plot_prauc = list()
    plot_threshold = list()
    for (model in c("source", "target")) {
      
      # Convert confidence scores to logit scale
      class_calibration_data$score = sapply(class_calibration_data[[paste0("confidence_", model)]], conf_to_logit)
      
      # Visualize precision and recall performance as a function of score
      data_sorted = class_calibration_data %>% mutate(score = logit_to_conf(score)) %>% arrange(desc(score))
      n_pos = sum(data_sorted$label_truth == 1)
      data_sorted = data_sorted %>% mutate(
        tp = cumsum(label_truth == 1),
        fp = cumsum(label_truth == 0),
        recall = tp / n_pos,
        precision = ifelse(tp + fp == 0, 1, tp / (tp + fp))  # handle division by zero
      )
      plt_pr = ggplot(data_sorted, aes(x = score)) +
        geom_line(aes(y = precision, color = "Precision"), linewidth = 1) +
        geom_line(aes(y = recall, color = "Recall"), linewidth = 1) +
        lims(x = c(0,1), y = c(0,1)) +
        labs(title = paste0(class_label), subtitle = paste0("model '", model, "'"), x = "Score", y = "Performance", color = "Metric")
      plot_pr[[model]] = plt_pr
      
      # Predicted positive/negative based on threshold
      calc_precision_recall = function(d, t_logit) {
        predicted = ifelse(d$score > t_logit, 1, 0)
        TP = sum(predicted == 1 & d$label_truth == 1)
        FP = sum(predicted == 1 & d$label_truth == 0)
        FN = sum(predicted == 0 & d$label_truth == 1)
        return(data.frame(
          precision = TP / (TP + FP),
          recall    = TP / (TP + FN)
        ))
      }
      
      # Calculate PR AUC and ROC AUC of classifier (from raw scores)
      auc_pr <- auc_roc <- NA
      t_maxp <- precision_tmaxp <- recall_tmaxp <- precision_tmin <- recall_tmin <- NA
      tryCatch(
        {
          auc_pr = pr.curve(scores.class0 = subset(class_calibration_data, label_truth == 1) %>% pull(score) %>% na.omit(),
                            scores.class1 = subset(class_calibration_data, label_truth == 0) %>% pull(score) %>% na.omit())$auc.integral
          
          auc_roc = roc.curve(scores.class0 = subset(class_calibration_data, label_truth == 1) %>% pull(score) %>% na.omit(),
                              scores.class1 = subset(class_calibration_data, label_truth == 0) %>% pull(score) %>% na.omit())$auc
          
          # Visualize precision recall AUC
          pr_curve = data_sorted %>% select(recall, precision) %>% arrange(recall)
          plt_prauc = ggplot(pr_curve, aes(x = recall, y = precision)) +
            geom_line(linewidth = 1) +
            lims(x = c(0,1), y = c(0,1)) +
            labs(
              title = paste0(class_label),
              subtitle = paste0("model '", model, "' PR AUC = ", round(auc_pr, 3)),
              x = "Recall", y = "Precision"
            )
          plot_prauc[[model]] = plt_prauc
          
          # Minimum threhsold to maximize precision (i.e. precision == 1)
          t_maxp = max(class_calibration_data$score[class_calibration_data$label_truth == 0], na.rm = TRUE)
          perf_tmaxp = calc_precision_recall(class_calibration_data %>% filter(!is.na(score)), t_maxp)
          precision_tmaxp = perf_tmaxp$precision
          recall_tmaxp    = perf_tmaxp$recall
          
          # Calculate performance at minimum threshold
          perf_tmin = calc_precision_recall(class_calibration_data %>% filter(!is.na(score)), conf_to_logit(threshold_min_classifier_score))
          precision_tmin = perf_tmin$precision
          recall_tmin    = perf_tmin$recall
          # message("  Performance at minimum threshold ", threshold_min_classifier_score, ":\n  Precision ", round(precision_tmin,3), "\n  Recall ", round(recall_tmin,3))
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
      # message("  PR-AUC ", round(auc_pr,3), ", ROC-AUC ", round(auc_roc,3))
      
      # Perform Platt scaling to determine threshold for desired probability of true positive
      threshold = Inf # default threshold is infinite (i.e. do not retain detections unless a species is validated)
      precision_threshold <- recall_threshold <- NA
      model_warning = FALSE
      tryCatch(
        {
          regression = glm(label_truth ~ score, class_calibration_data %>% filter(!is.na(score)), family = binomial(link = "logit"))
          intercept   = as.numeric(coef(regression)[1])
          coefficient = as.numeric(coef(regression)[2])
          threshold_logit = (log(tp_min_prob / (1 - tp_min_prob)) - intercept) / coefficient # logit scale
          threshold       = logit_to_conf(threshold_logit) # confidence scale
          
          # message("  ", round(threshold, 3), " threshold to achieve Pr(TP)>=", tp_min_prob)
          
          # Calculate estimated precision and recall at this threshold
          perf_t = calc_precision_recall(class_calibration_data %>% filter(!is.na(score)), threshold_logit)
          precision_threshold  = perf_t$precision
          recall_threshold     = perf_t$recall
          # message("  Performance at threshold ", round(threshold,3), ":\n  Precision ", round(precision_threshold,3), "\n  Recall ", round(recall_threshold,3))
          
          ## Visualize the logistic regression and data
          x_range = seq(min(class_calibration_data$score, na.rm = TRUE), max(class_calibration_data$score, na.rm = TRUE), length.out = 100)
          
          # Predict fitted values and calculate confidence intervals on probability scale
          pred = predict(regression, newdata = data.frame(score = x_range), type = "link", se.fit = TRUE)
          
          # Calculate confidence intervals on probability scale
          z = qnorm(0.975)  # 95% CI
          upper  = logit_to_conf(pred$fit + z * pred$se.fit)
          lower  = logit_to_conf(pred$fit - z * pred$se.fit)
          y_pred = logit_to_conf(pred$fit)
          
          regression_df = data.frame(score = x_range, prob = y_pred, lower = lower, upper = upper, warning = model_warning)
          plt_threshold = ggplot() +
            geom_vline(xintercept = 0, color = "gray", linetype = "solid", linewidth = 0.5) +
            geom_point(data = class_calibration_data, aes(x = score, y = label_truth), shape = 1, alpha = 0.5) +
            geom_hline(yintercept = tp_min_prob, color = "black", linetype = "dashed", linewidth = 0.5) +
            geom_vline(xintercept = threshold_logit, color = "black", linetype = "dashed", linewidth = 0.5) +
            geom_ribbon(data = regression_df, aes(x = score, ymin = lower, ymax = upper, fill = warning), alpha = 0.2) +
            geom_line(data = regression_df, aes(x = score, y = prob, color = warning), linewidth = 1) +
            scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red")) +
            scale_fill_manual(values = c("FALSE" = "blue", "TRUE" = "red")) +
            lims(x = c(-11.52,11.52), y = c(0,1)) +
            labs(
              x = "Score (logit)", y = "True positive probability",
              title = class_label, subtitle = paste0("model '", model, "', threshold p(TP) ≥ ", tp_min_prob, " = ", round(threshold, 3))
            ) + theme(legend.position = "bottom")
          plot_threshold[[model]] = plt_threshold
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
      
      stats = bind_rows(stats, tibble(
        class_label          = class_label,                                  # Class
        model                = model,                                        # Classifier model
        n_tp                 = sum(class_calibration_data$label_truth == 1), # Number of true positive examples
        n_tn                 = sum(class_calibration_data$label_truth == 0), # Number of true negative examples
        auc_pr               = auc_pr,                                       # Precision-recall AUC
        auc_roc              = auc_roc,                                      # Receiver operating curve AUC
        warning              = model_warning,                                # Issue fitting logistic regression
        tp_min_p             = tp_min_prob,                                  # Requested minimum probability of true positive
        t_tp                 = threshold,                                    # Calibrated threshold to achieve requested minimum probability of TP
        precision_ttp        = precision_threshold,                          # Precision at the calibrated threshold
        recall_ttp           = recall_threshold,                             # Recall at the calibrated threshold
        t_min                = threshold_min_classifier_score,               # Naive minimum confidence score threshold
        precision_tmin       = precision_tmin,                               # Precision at the naive threshold
        recall_tmin          = recall_tmin,                                  # Recall at the naive threshold
        t_maxp               = logit_to_conf(t_maxp),                        # Threshold to achieve perfect precision (1) while maximizing recall
        precision_tmaxp      = precision_tmaxp,                              # Precision at the perfect precision threshold
        recall_tmaxp         = recall_tmaxp                                  # Recall at the perfect precision threshold
      ))
      
    }
    class_plots = list(plot_pr, plot_prauc, plot_threshold)
    names(class_plots) = c("pr", "prauc", "threshold")
    plots[[class_label]] = class_plots
  }
  calibration_data[["stats"]] = stats
  calibration_data[["plots"]] = plots
  return(calibration_data)
}

# Calibrate combined data -----------------------------------------------------------------------------

message("Calibrating each class:")
calibration_data = calibrate(predictions, annotations, calibration_class)

message("Calibration stats:")
stats = calibration_data[["stats"]] %>% mutate(across(where(is.numeric), ~ round(., 3)))
print(stats, n = Inf)

# saveRDS(calibration_data, path_calibration_data_cache)
# message(crayon::green("Cached", path_calibration_data_cache))

# Determine optimal thresholds per class -----------------------------------------------------------------------------

message("Determining optimal thresholds for each class")

# Drop additional abiotic and biotic classes
class_thresholds = stats
class_thresholds = class_thresholds %>% filter(!str_starts(class_label, "abiotic"))
class_thresholds = class_thresholds %>% filter(!str_starts(class_label, "biotic"))

calibration_results = tibble()
for (l in unique(class_thresholds$class_label)) {
  l_data = class_thresholds %>% filter(class_label == l)
  
  # Determine the model with the highest PR AUC
  l_data = l_data %>% slice_max(order_by = auc_pr, n = 1, with_ties = FALSE)
  
  # The "method" column describes one of four methods to determine the threshold. Most classes use "probabilistic",
  # which leverages logistic regression (Platt scaling) to determine the threshold at which Pr(TP)>=0.99. However,
  # some classes cannot fit such a regression because they:
  #   1) report perfect performance from our calibration data (the random sample contained no true negatives).
  #      These classes default to A) an optimized threshold that maximizes recall while maintaining perfect precision
  #      ("precision"); or B) a naive threshold at the decision boundary 0.5 ("naive"), whichever threshold is higher.
  #   2) report exceptionally poor performance (the random sample contained very few true positives). These classes must
  #      be manually reviewed by a human to determine true presence/absence (model = "manual").
  #   3) are very rare (there are not enough samples available to meaningfully evaluate performance). These classes also
  #      must be manually reviewed by a human to determine true presence/absence (model = "manual").
  
  t = Inf
  m = "manual"
  p <- r <- NA
  model_warning = FALSE
  tryCatch(
    {
      if (l_data$warning) {
        # For classes that did not converge, use max of t_maxp or t_min
        if (l_data$t_maxp > l_data$t_min) {
          t = l_data$t_maxp
          p = l_data$precision_tmaxp
          r = l_data$recall_tmaxp
          m = "precision"
        } else {
          t = l_data$t_min
          p = l_data$precision_tmin
          r = l_data$recall_tmin
          m = "naive"
        }
      } else {
        # Use the probabilistic threshold
        t = l_data$t_tp
        p = l_data$precision_ttp
        r = l_data$recall_ttp
        m = "probabilistic"
      }
    },
    warning = function(w) {
      message(crayon::yellow("WARNING:", l, conditionMessage(w)))
    },
    error = function(e) {
      message(crayon::yellow("WARNING:", l, conditionMessage(e)))
    }
  )
  calibration_results = bind_rows(calibration_results, l_data %>%
                                     select(class_label, model, n_tp, n_tn, auc_pr, auc_roc) %>% mutate(threshold = t, precision = p, recall = r, method = m))
}
calibration_results = calibration_results %>% mutate(across(where(is.numeric), ~ round(., 3)))
calibration_results = left_join(calibration_results, class_labels %>% rename(class_label = label), by = "class_label")

# The following classes are rare or insufficiently performant, and must only be
# detected through manual review (i.e. set their thresholds to Inf, and method to "manual")
manual_classes = c(
  "brown-headed cowbird",
  "cliff swallow",
  "downy woodpecker",
  "eurasian collared-dove",
  "great horned owl",
  "merlin",
  "northern goshawk",
  "northern rough-winged swallow",
  "ruffed grouse",
  "townsend's solitaire",
  "western bluebird",
  "yellow warbler"
)
calibration_results = calibration_results %>%
  mutate(
    threshold = if_else(common_name %in% manual_classes, Inf, threshold),
    precision = if_else(common_name %in% manual_classes, NA, precision),
    recall    = if_else(common_name %in% manual_classes, NA, recall),
    method    = if_else(common_name %in% manual_classes, "manual", method)
  )

# TODO: The following classes require exhaustive manual validation
# cassin's vireo
# great horned owl
# marbled murrelet
# western bluebird

# Format calibration results for output and tables
calibration_results_clean = calibration_results
absent_classes = calibration_results_clean %>% filter(n_tp == 0)
absent_classes$model = "source"
absent_classes$threshold = Inf
absent_classes$method = "manual"
calibration_results_clean = bind_rows(calibration_results_clean %>% filter(n_tp > 0), absent_classes)
calibration_results_clean = calibration_results_clean %>% select(scientific_name, common_name, everything()) %>% select(-class_label)

message("Calibration results:")
print(calibration_results_clean, n = Inf)

write_csv(calibration_results, path_calibration_results_raw_cache, na = "-")
message(crayon::green("Cached", path_calibration_results_raw_cache))

write_csv(calibration_results_clean, path_calibration_results_clean_cache, na = "-")
message(crayon::green("Cached", path_calibration_results_clean_cache))

# Data inspection -----------------------------------------------------------------------------

if (FALSE) {
  
  # Inspect a class
  calibration_data[["stats"]] %>% filter(class_label == calibration_class) %>% mutate(across(where(is.numeric), ~ round(., 2)))
  calibration_data[["plots"]][[l]][["pr"]]
  calibration_data[["plots"]][[l]][["prauc"]]
  calibration_data[["plots"]][[l]][["threshold"]]
  
  # Find specific annotations
  class_annotations = annotations %>% select(file, all_of(l)) %>% rename(label_truth = !!sym(l))
  class_annotations = class_annotations %>% filter(label_truth != "unknown") %>% mutate(label_truth = as.integer(label_truth))
  table(class_annotations$label_truth)
  class_predictions = predictions %>% filter(label_predicted == l)
  class_calibration_data = left_join(class_annotations, class_predictions, by = c("file"))
  View(class_calibration_data %>% arrange(desc(confidence_source)))
  
}
