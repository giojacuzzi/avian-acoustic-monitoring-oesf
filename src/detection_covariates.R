##############################################################################
# Quantify covariates on detection
#
##############################################################################

path_data_out = "data/cache/detection_covariates/data_detection.rds"

library(dplyr)

files = list.files("data/environment/climate", pattern = "\\.csv$", full.names = TRUE)

data = files %>%
  lapply(function(file) {
    d = read_csv(file, skip = 6, show_col_types = FALSE)
    d$site = str_remove(basename(file), "\\.csv$")
    d
  }) %>% bind_rows() %>% select(site, everything()) %>% janitor::clean_names()

message('Saving detection covariate data cache ', path_data_out)
dir.create(dirname(path_data_out), recursive = TRUE, showWarnings = FALSE)
saveRDS(data, path_data_out)
