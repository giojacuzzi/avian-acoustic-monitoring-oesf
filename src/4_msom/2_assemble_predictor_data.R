# 2_assemble_predictor_data.R ####################################################################################
# Derive MSOM-ready data
# - TODO
#
# CONFIG:

#
# OUTPUT:
out_cache_dir  = "data/cache/4_msom/2_assemble_predictor_data"
out_path_occurrence_predictor_plot_data      = paste0(out_cache_dir, "/occurrence_predictor_plot_data.rds")
out_path_occurrence_predictor_homerange_data = paste0(out_cache_dir, "/occurrence_predictor_homerange_data.rds")
out_path_detection_predictor_data            = paste0(out_cache_dir, "/detection_predictor_data.rds")
# TODO
#
# INPUT:
path_y                    = "data/cache/4_msom/1_assemble_detection_arrays/y.rds"
path_xyday                = "data/cache/4_msom/1_assemble_detection_arrays/xyday.rds"
path_predictors_detection = "data/cache/detection_covariates/data_detection.rds"
path_plot_scale_data      = "data/cache/habitat_vars/data_plot_scale_sites.rds"
path_homerange_scale_data = "data/cache/habitat_vars/data_homerange_scale_sites.rds"
##################################################################################################################

source("src/global.R")

if (!dir.exists(out_cache_dir)) dir.create(out_cache_dir, recursive = TRUE)

# Load dependencies ----------------------------------------------------------------------------------------------

message("Loading observation y array from ", path_y)
y = readRDS(path_y)

message("Loading survey yday array from ", path_xyday)
x_yday = readRDS(path_xyday)

species = dimnames(y)$species
seasons = dimnames(y)$season
surveys = dimnames(y)$survey
sites   = dimnames(y)$site

message("Loading detection predictor data")
detection_data = readRDS(path_predictors_detection) %>% mutate(year = as.character(year)) %>% rename(season = year)

message("Loading local plot scale data from ", path_plot_scale_data)
occurrence_predictor_plot_data_local = readRDS(path_plot_scale_data) %>% sf::st_drop_geometry() %>% arrange(site) %>% mutate(site = tolower(site)) %>%
  select(site, elev, dist_road_paved, dist_road_all, dist_watercourse_major, dist_watercourse_all)

s = 'median' # TODO: do occurrence_predictor_plot_data_rs at PLOT scale, not median
message("Loading remote sensing plot scale data from '", s, "' layer ", path_homerange_scale_data)
occurrence_predictor_plot_data_rs = readRDS(path_homerange_scale_data)[[s]] %>% arrange(site) %>% mutate(site = tolower(site)) %>%
  select(site, homerange_qmd_mean, homerange_treeden_all_mean)

message("Combining all plot scale occurrence predictor data")
occurrence_predictor_plot_data = occurrence_predictor_plot_data_local %>% full_join(occurrence_predictor_plot_data_rs, by = "site")

message("Loading median homerange scale data")
occurrence_predictor_homerange_data = readRDS(path_homerange_scale_data)
occurrence_predictor_homerange_data = map(occurrence_predictor_homerange_data, ~ .x %>% arrange(site) %>% mutate(site = tolower(site)))
names(occurrence_predictor_homerange_data) = tolower(names(occurrence_predictor_homerange_data))

# Discard any sites not surveyed ---------------------------------------------------------------------------------

# Discard sites with no survey observations
survey_counts = apply(!is.na(x_yday), c(1, 3), sum) # surveys per site per season
site_survey_counts = rowSums(survey_counts)         # surveys per site across all seasons
sites_not_surveyed = names(site_survey_counts)[site_survey_counts == 0]
if (length(sites_not_surveyed) > 0) {
  message("Discarding ", length(sites_not_surveyed), " sites with no surveys")
  x_yday = x_yday[!(dimnames(x_yday)$site %in% sites_not_surveyed), , , drop = FALSE]
  y      = y[!(dimnames(y)$site %in% sites_not_surveyed), , , , drop = FALSE]
}

# Discard sites with no environmental data
sites_missing_environmental_data = setdiff(sites, occurrence_predictor_plot_data$site)
if (length(sites_missing_environmental_data) > 0) {
  message("Discarding ", length(sites_missing_environmental_data), " sites with missing environmental data:")
  print(sites_missing_environmental_data)
  x_yday = x_yday[!(dimnames(x_yday)$site %in% sites_missing_environmental_data), , , drop = FALSE]
  y      = y[!(dimnames(y)$site %in% sites_missing_environmental_data), , , , drop = FALSE]
}

seasons = dimnames(y)$season
surveys = dimnames(y)$survey
sites   = dimnames(y)$site

message("Total number of sites: ", length(sites))
sites_per_season = sapply(seasons, function(season) {
  mat = x_yday[, , season, drop = FALSE]
  sum(rowSums(!is.na(mat)) > 0)
})
message("Sites per season:")
print(sites_per_season)

for (t in seasons) {
  mat = x_yday[, , t, drop = FALSE]
  n_surveys_per_site = apply(!is.na(mat), 1, sum)
  message("Season ", t, ": ", sum(n_surveys_per_site), " total sampling periods (surveys) conducted across ", sites_per_season[t], " sampling units (sites)")
  message("Sampling periods (surveys) conducted per site: median ", median(n_surveys_per_site), ", range ", min(n_surveys_per_site), "–", max(n_surveys_per_site))
  print(table(n_surveys_per_site))
}

# Assemble all occurrence predictor data -------------------------------------------------------------------------

# Conform occurrence data sites with surveyed sites
occurrence_predictor_plot_data = occurrence_predictor_plot_data %>% filter(site %in% sites) # discard data for irrelevant sites
stopifnot(dimnames(y)[["site"]] == occurrence_predictor_plot_data$site)    # check that covariate data are aligned with observation matrix by site
for (r in names(occurrence_predictor_homerange_data)) {
  occurrence_predictor_homerange_data[[r]] = occurrence_predictor_homerange_data[[r]] %>% filter(site %in% sites)
  stopifnot(dimnames(y)[["site"]] == occurrence_predictor_homerange_data[[r]]$site)
}

# Cache
saveRDS(occurrence_predictor_plot_data, out_path_occurrence_predictor_plot_data)
message(crayon::green("Cached occurrence predictor plot data to", out_path_occurrence_predictor_plot_data))
saveRDS(occurrence_predictor_homerange_data, out_path_occurrence_predictor_homerange_data)
message(crayon::green("Cached occurrence predictor homerange data to", out_path_occurrence_predictor_homerange_data))

# Assemble all detection predictor data --------------------------------------------------------------------------

# Reformat x_yday to long
x_yday_long = as.data.frame.table(x_yday, responseName = "yday") %>%
  mutate(
    survey = as.integer(as.character(survey)),
    season = as.character(season)
  ) %>% filter(!is.na(yday))

# Join with detection_data
x_env = x_yday_long %>%
  left_join(detection_data %>% 
              select(site, season, yday, tmax_deg_c, prcp_mm_day),
            by = c("site", "season", "yday"))

# Indices for array assignment
site_idx = match(x_env$site, dimnames(x_yday)[[1]])
survey_idx = x_env$survey
season_idx = match(x_env$season, dimnames(x_yday)[[3]])

# Assign values directly into new arrays
x_tmax = x_yday
x_tmax[cbind(site_idx, survey_idx, season_idx)] = x_env$tmax_deg_c
x_prcp = x_yday
x_prcp[cbind(site_idx, survey_idx, season_idx)] = x_env$prcp_mm_day

# Assemble detection covariate data
detection_predictor_data = list(
  yday        = x_yday,
  prcp_mm_day = x_prcp,
  tmax_deg_c  = x_tmax
)

# Cache
saveRDS(detection_predictor_data, out_path_detection_predictor_data)
message(crayon::green("Cached detection predictor data to", out_path_detection_predictor_data))
