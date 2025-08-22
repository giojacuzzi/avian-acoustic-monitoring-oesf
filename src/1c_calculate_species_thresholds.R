# Calculate species-specific probabilistic score thresholds using Platt scaling (Platt 2000, Wood and Kahl 2024)

library(dplyr)
library(ggplot2)
library(PRROC)

p_tp = 0.95 # Desired probability of true positive
display_plots = FALSE

conf_to_logit = function(c) {
  # guard against undefined logit for exceptionally low/high scores beyond model rounding limits
  c = min(max(c, 0.00001), 1.0 - 0.00001)
  return(log(c / (1 - c)))
}

logit_to_conf = function(l) {
  return(1 / (1 + exp(-l)))
}

message("Calculating species-specific thresholds and performance metrics (current time ", time_start <- format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ")")

models = c("source", "target")
for (model in models) {
  
  message("Calculating species-specific thresholds and performance metrics for model: ", model)
  
  results = data.frame()
  path_in = paste0("data/debug/predictions_", model, ".parquet")
  predictions = arrow::read_parquet(path_in) %>% arrange(file, label_predicted)
  labels = sort(unique(predictions$label_predicted))
  n_labels = length(labels)
  
  for (l in 1:n_labels) {
    
    # Obtain predictions for this label
    label = labels[l]
    message(label)
    label_predictions = predictions[predictions$label_predicted == label, ]
    
    # If the label is truly present in a given prediction, mark as 1 (otherwise 0)
    label_predictions$label_truth[
      label_predictions$label_truth == label_predictions$label_predicted
    ] <- 1
    label_predictions$label_truth = as.numeric(label_predictions$label_truth)
    
    # Convert confidence scores to logit scale
    label_predictions = label_predictions %>% mutate(score = sapply(confidence, conf_to_logit))
    
    data = data.frame(
      score = label_predictions$score,
      label_truth = label_predictions$label_truth
    )
    stopifnot(all(unique(data$label_truth) %in% c(0, 1)))
    
    # Calculate PR AUC and ROC AUC of classifier (from raw scores)
    auc_pr = pr.curve(scores.class0 = subset(data, label_truth == 1) %>% pull(score),
                      scores.class1 = subset(data, label_truth == 0) %>% pull(score))$auc.integral
    
    auc_roc = roc.curve(scores.class0 = subset(data, label_truth == 1) %>% pull(score),
                        scores.class1 = subset(data, label_truth == 0) %>% pull(score))$auc
    
    message("Classifier performance PR-AUC ", round(auc_pr,3), ", ROC-AUC ", round(auc_roc,3))
    
    # Perform platt scaling to determine threshold for desired probability of true positive
    model_warning = FALSE
    tryCatch(
      {
        regression = glm(label_truth ~ score, data, family = binomial(link = "logit"))
      },
      warning = function(w) {
        model_warning <<- TRUE
        message(crayon::yellow("WARNING", conditionMessage(w)))
      },
      error = function(e) {
        stop(crayon::red("ERROR", conditionMessage(e)))
      }
    )
    
    intercept   = as.numeric(coef(regression)[1])
    coefficient = as.numeric(coef(regression)[2])
    threshold_logit = (log(p_tp / (1 - p_tp)) - intercept) / coefficient # logit scale
    threshold_conf = logit_to_conf(threshold_logit) # confidence scale
    
    message("Threshold to achieve TP probability ", p_tp, ": ", round(threshold_conf, 3))
    
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
    perf_t = calc_precision_recall(data, threshold_logit)
    message("Performance at threshold ", round(threshold_conf,3), ": precision ", round(perf_t$precision,3), ", recall ", round(perf_t$recall,3))
    # "For most false-positive models in our study, using a mid-range threshold of 0.50 or above generally yielded stable estimates." (Katsis et al. 2025)
    perf_t0.5 = calc_precision_recall(data, conf_to_logit(0.5))
    message("Performance at threshold 0.5: precision ", round(perf_t0.5$precision,3), ", recall ", round(perf_t0.5$recall,3))
    
    results = rbind(results, data.frame(
      model         = model,
      species       = label,
      n_pos         = sum(data$label_truth == 1),
      n_neg         = sum(data$label_truth == 0),
      auc_pr        = auc_pr,
      auc_roc       = auc_roc,
      warning       = model_warning,
      p_tp          = p_tp,
      t_conf_tp     = threshold_conf,
      precision_tp  = perf_t$precision,
      recall_tp     = perf_t$recall,
      precision_0.5 = perf_t0.5$precision,
      recall_0.5    = perf_t0.5$recall
    ))
    
    if (display_plots) {
      # Visualize the logistic regression and data
      x_range = seq(min(data$score), max(data$score), length.out = 100)
      
      # Predict fitted values and calculate confidence intervals on probability scale
      pred = predict(model, newdata = data.frame(score = x_range), type = "link", se.fit = TRUE)
      
      # Calculate confidence intervals on probability scale
      z = qnorm(0.975)  # 95% CI
      upper  = logit_to_conf(pred$fit + z * pred$se.fit)
      lower  = logit_to_conf(pred$fit - z * pred$se.fit)
      y_pred = logit_to_conf(pred$fit)
      
      regression_df = data.frame(score = x_range, prob = y_pred, lower = lower, upper = upper, warning = model_warning)
      plt = ggplot() +
        geom_vline(xintercept = 0, color = "gray", linetype = "solid", size = 0.5) +
        geom_point(data = data, aes(x = score, y = label_truth), shape = 1, alpha = 0.5) +
        geom_hline(yintercept = p_tp, color = "black", linetype = "dashed", size = 0.5) +
        geom_vline(xintercept = threshold_logit, color = "black", linetype = "dashed", size = 0.5) +
        geom_ribbon(data = regression_df, aes(x = score, ymin = lower, ymax = upper, fill = warning), alpha = 0.2) +
        geom_line(data = regression_df, aes(x = score, y = prob, color = warning), size = 1) +
        scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red")) +
        scale_fill_manual(values = c("FALSE" = "blue", "TRUE" = "red")) +
        labs(
          x = "Score (logit)", y = "True positive probability",
          title = paste0(label, " (", model, " model)"), subtitle = paste0("Threshold p(TP) ≥ ", p_tp, " = ", round(threshold_conf, 3))
        ); print(plt)
    }
  }
  
  print(results)
  
  path_out = paste0("data/cache/1_calculate_species_thresholds/species_thresholds_", model, ".csv")
  if (!dir.exists(dirname(path_out))) dir.create(dirname(path_out), recursive = TRUE)
  write.csv(results, path_out, row.names = FALSE)
  message(crayon::green("Cached species threshold data to ", path_out))
}

message("Loading cached species threshold data for comparison between models")
species_thresholds_source = read.csv("data/cache/1_calculate_species_thresholds/species_thresholds_source.csv")
species_thresholds_target = read.csv("data/cache/1_calculate_species_thresholds/species_thresholds_target.csv")

# Identify species only supported by the source model
species_only_in_source = species_thresholds_source %>% filter(!species %in% species_thresholds_target$species)

# Compare shared species performance
source_metrics = species_thresholds_source %>% select(species, auc_pr, recall_tp) %>%
  rename(auc_pr_source = auc_pr, recall_tp_source = recall_tp)

target_metrics = species_thresholds_target %>% select(species, auc_pr, recall_tp) %>%
  rename(auc_pr_target = auc_pr, recall_tp_target = recall_tp)
shared_metrics = inner_join(source_metrics, target_metrics, by = "species")

# Compute delta performance values
print(shared_metrics %>% mutate(delta_auc_pr = auc_pr_target - auc_pr_source, delta_recall_tp = recall_tp_target - recall_tp_source) %>% select(species, delta_auc_pr, delta_recall_tp) %>% arrange(desc(delta_auc_pr)))

species_thresholds = rbind(species_thresholds_target, species_thresholds_source)
path_out = paste0("data/cache/1_calculate_species_thresholds/species_thresholds.csv")
if (!dir.exists(dirname(path_out))) dir.create(dirname(path_out), recursive = TRUE)
write.csv(species_thresholds, path_out, row.names = FALSE)
message(crayon::green("Cached species threshold data to ", path_out))

message(crayon::green("Finished calculating species-specific thresholds and performance metrics (", round(as.numeric(difftime(Sys.time(), time_start, units = 'mins')), 2), " min )"))

