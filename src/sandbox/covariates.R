library(tidyverse)

plot_data_path = 'data/environment/PAM_PreHarvest_Habitat_results_DD_WD_TM.xlsx'

# Load and clean plot-level data
plot_data = list(
  readxl::read_xlsx(plot_data_path, sheet = 2, skip = 1) %>% janitor::clean_names(),
  readxl::read_xlsx(plot_data_path, sheet = 4, skip = 1) %>% janitor::clean_names(),
  readxl::read_xlsx(plot_data_path, sheet = 5, skip = 1) %>% janitor::clean_names(),
  readxl::read_xlsx(plot_data_path, sheet = 6, skip = 1) %>% janitor::clean_names()
) %>%
  reduce(full_join, by = c("station", "strata")) %>%
  rename(site = station, stratum = strata) %>%
  mutate(
    tag = str_extract(site, "_.*$") %>% str_remove("^_"),
    site = str_remove(site, "_.*$")
  )

# Coalesce duplicate site entries
plot_data = plot_data %>% group_by(site) %>% summarise(across(everything(), ~ coalesce(.[!is.na(.)][1], NA)))

# NA entries by column
na_values =  t(as.data.frame(plot_data %>% summarise(across(everything(), ~ sum(is.na(.))))))

# TODO: Override NAs with 0 for select columns?
# selected_cols = c(avg_dbh_cm, avg_height_cm)
# plot_data %>% mutate(across(all_of(selected_cols), ~replace(., is.na(.), 0)))

cor_data = plot_data %>% select(-c(site, stratum, tag))

# Drop "unknown" variables
cor_data = cor_data %>% select(-ends_with("_un"))
# Drop variables with sparse observations and low variance
cor_data = cor_data %>% select(-c(avg_dbh_cm_abam, avg_dbh_cm_pisi, avg_height_m_abam, avg_height_m_alru, avg_height_m_pisi, avg_hlc_abam, avg_hlc_alru, avg_hlc_pisi, avg_llc_abam, avg_llc_alru, avg_llc_pisi, avg_lcr_abam, avg_lcr_alru, avg_lcr_pisi))
# i = 75
# cor(cor_data[1:i], use = "pairwise.complete.obs", method = "pearson")
# colnames(cor_data)[i]

cor_matrix = cor(cor_data, use = "pairwise.complete.obs", method = "pearson")

corrplot::corrplot(cor_matrix, method = "color", type = "upper",
                   tl.cex = 0.2, tl.col = "black", tl.srt = 45, # text color and rotation
                   addCoef.col = "black", # add correlation coefficients
                   number.cex = 0.2, # size of the numbers
                   diag = FALSE) # hide diagonal

correlation_threshold = function(cm, threshold) {
  highly_correlated_idx = which(abs(cm) >= threshold & abs(cm) < 1, arr.ind = TRUE)
  highly_correlated_idx = highly_correlated_idx[highly_correlated_idx[,1] < highly_correlated_idx[,2], ]
  data.frame(
    rownames(cm)[highly_correlated_idx[,1]],
    colnames(cm)[highly_correlated_idx[,2]],
    correlation = cm[highly_correlated_idx]
  )
}

correlation_threshold(cor_matrix, 0.8)

# Remove the following variables to reduce multicollinearity
cor_data_reduced = cor_data %>% select(-c(
  # perc_deciduous_shrub
  cv_deciduous_shrub,
  # per_cover_all_shrubs
  per_cover_total_understory, per_cover_medium_shrub, shrub_layer_vol,
  # cv_all_shrubs
  cv_shrub_layer_vol,
  # avg_understory_heights
  max_understory_heights,
  # cv_avg_understory_heights
  cv_max_understory_heights,
  # avg_dbh_cm_psme, avg_dbh_cm_tshe
  avg_dbh_cm_all, avg_height_m_psme,
  # ba_ha_psme
  large_per_hectare_psme,
  # ba_ha_thpl
  avg_dbh_cm_thpl,
  # ba_ha_abam
  large_per_hectare_abam,
  # snags_ha
  decay_1_snags_ha,
  # vol_alldown_m3
  vol_rottendown_m3
))
# Discard all HLC and LCR measurements, but retain LLC
# These measurements are highly correlated and some are interdependent.
# LLC is a directly interpretable measurement of the amount of crown habitat
cor_data_reduced = cor_data_reduced %>% select(-matches("_hlc_|_lcr_"))
# Drop all species-specific large tree average heights and LLC measurements
cor_data_reduced = cor_data_reduced %>% select(-starts_with("avg_height_m_"), avg_height_m_all)
cor_data_reduced = cor_data_reduced %>% select(-matches("_llc_"), avg_llc_all)

cor_matrix_reduced = cor(cor_data_reduced, use = "pairwise.complete.obs", method = "pearson")
corrplot::corrplot(cor_matrix_reduced, method = "color", type = "upper",
                   tl.cex = 0.2, tl.col = "black", tl.srt = 45, # text color and rotation
                   addCoef.col = "black", # add correlation coefficients
                   number.cex = 0.2, # size of the numbers
                   diag = FALSE) # hide diagonal

correlation_threshold(cor_matrix_reduced, 0.7)



