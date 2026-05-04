source("src/global.R")

# Q: Among thinned stands, are higher QMD plots more recently thinned?

path_occurrence_predictor_plot_data = "data/cache/4_msom/1_assemble_msom_data/V3_occurrence_predictor_plot_data.rds"
message("Loading occurrence predictor plot scale data from ", path_occurrence_predictor_plot_data)
occurrence_predictor_plot_data = readRDS(path_occurrence_predictor_plot_data)[[1]]

poly_thinning_treatment = st_read('data/environment/GIS Data/Forest Development Strata/ThinAfter94NoHarvSinceClipByInitBuf3.shp') %>% 
  st_transform(crs_m)

poly_with_aru <- st_filter(poly_thinning_treatment, aru_sites)

aru_in_poly <- st_join(aru_sites, poly_thinning_treatment["FMA_FY"]) |>
  filter(!is.na(FMA_FY))
#
poly_with_aru$FMA_FY

data = left_join(aru_in_poly %>% st_drop_geometry(), occurrence_predictor_plot_data %>% select(site, qmd_6_mean))

data$years_since_thin = max(data$FMA_FY) - data$FMA_FY

cor.test(~ years_since_thin + qmd_6_mean, data = na.omit(data[, c("years_since_thin", "qmd_6_mean")]), method = "spearman")

ggplot(data, aes(x = years_since_thin, y = qmd_6_mean)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE)
