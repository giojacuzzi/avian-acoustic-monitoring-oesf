#####################################################################################
# Predict habitat use across the landscape
#
# INPUT:
path_msom = "data/cache/models/V4_msom_V4_nofp_nofp_all.rds"

source("src/global.R")

# Load MSOM ---------------------------------------------------------------------------

message("Loading data for multi-species occupancy model ", path_msom)
model_data = readRDS(path_msom)

msom_summary = model_data$msom_summary
msom = model_data$msom
groups = model_data$groups %>% arrange(common_name)
sites = model_data$sites
species = model_data$species
seasons = model_data$seasons
stages = model_data$stages
strata = factor(stages$stratum_4, levels = c("standinit", "compex", "thin", "mature"))

z = msom$sims.list$z
n_iter    = dim(z)[1]
n_sites   = dim(z)[2]
n_seasons = dim(z)[3]
n_species = dim(z)[4]

# Load plot and homerange scale data for all points -----------------------------------

homerange_scale = "pileated woodpecker"

path_homerange = "data/cache/7_landscape/OPT_landscape_data_homerange_scale_2020_clean_strata_4_landscape.rds"
path_plot      = "data/cache/7_landscape/OPT_landscape_data_plot_scale_2020_clean_strata_4_landscape.rds"

data_homerange = read_rds(path_homerange)
names(data_homerange)
str(data_homerange[[homerange_scale]])

points_plot = read_rds(path_plot)
str(points_plot)

points_homerange = st_as_sf(data_homerange[[homerange_scale]], coords = c("x", "y"), crs = crs_m)

points = points_homerange %>% st_join(points_plot, join = st_equals)

# mapview(points)

point_data = cbind(st_coordinates(points), st_drop_geometry(points))
str(point_data)

vars = c(
  "pcnt_standinit", "pcnt_compex", "pcnt_mature", "pcnt_thin",
   "elevation", "dist_road_paved", "dist_watercourse_major")

rasters = lapply(vars, function(v) {
  rast(point_data[, c("X", "Y", v)], type = "xyz", crs = crs_m_rast)
})
names(rasters) = vars

mapview(rasters[["pcnt_thin"]])

rast_stack = rast(rasters)
plot(rast_stack)

# Predict use probability -----------------------------------


