##############################################################################
# Finalize candidate set of variables for occurrence at the home range scale
#
##############################################################################
library(mapview)
library(terra)
library(sf)
library(dplyr)
library(ggrepel)
theme_set(theme_minimal())

crs_m = 32610 # coordinate reference system EPSG:32610, UTM Zone 10N, (meters)

path_data_plot_scale = "data/cache/occurrence_covariates/data_plot_scale.rds"
message('Loading plot scale data from cache ', path_data_plot_scale)
data_plot_scale = readRDS(path_data_plot_scale)
data_plot_scale = data_plot_scale %>% filter(hs == TRUE)

path_data_homerange_scale = "data/cache/occurrence_covariates/data_homerange_scale.rds"
message('Loading homerange scale data from cache ', path_data_homerange_scale)
data_homerange_scale = readRDS(path_data_homerange_scale)

path_rast_cover_clean = "data/cache/occurrence_covariates/rast_cover_clean.tif"
message('Loading raster cover data from cache ', path_rast_cover_clean)
rast_cover_clean = rast(path_rast_cover_clean)

watersheds = st_read('data/environment/GIS Data/watersheds/watersheds.shp') %>% janitor::clean_names() %>% st_transform(crs = crs_m)

mapview(rast_cover_clean,
        alpha.regions = 1.0,
        col.regions = c('#90c6bd', '#3c8273', '#d8c18a', '#9b652b', 'darkgray', '#6495ED')) +
  mapview(data_plot_scale) +
  mapview(layer.name = 'plot', st_buffer(data_plot_scale, 100), alpha.regions = 0.0, lwd = 2) +
  mapview(layer.name = 'median', st_buffer(data_plot_scale, 204), alpha.regions = 0.0, lwd = 2) +
  mapview(layer.name = 'mean', st_buffer(data_plot_scale, 468), alpha.regions = 0.0, lwd = 2) +
  mapview(layer.name = 'nearmax', st_buffer(data_plot_scale, 1000), alpha.regions = 0.0, lwd = 2)

### Principal component analysis ##############################################################



### Collinearity analysis ##############################################################

# TODO: Should I be including all spatial scales together?

# Start with all candidate variables
data_homerange_scale_focus = data_homerange_scale$Mean

data_homerange_scale_plot_candidates = data_homerange_scale_focus %>% st_drop_geometry() %>% select(where(is.numeric))
length(names(data_homerange_scale_plot_candidates))
data_homerange_scale_plot_candidates = as.data.frame(data_homerange_scale_plot_candidates)
data_homerange_scale_plot_candidates = sapply(data_homerange_scale_plot_candidates, as.numeric)

cor_matrix = cor(data_homerange_scale_plot_candidates, use = "pairwise.complete.obs", method = "pearson")
cor_matrix[lower.tri(cor_matrix, diag = TRUE)] = NA
(collinearity_candidates = subset(as.data.frame(as.table(cor_matrix)), !is.na(Freq) & abs(Freq) >= 0.8))
