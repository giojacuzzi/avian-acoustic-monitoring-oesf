# Predict occupancy across a landscape

source("src/global.R")

library(dplyr)
library(sf)
library(terra)
library(ggplot2)
library(mapview)
options(mapview.maxpixels = 2117676)

crs_m = 32610 # coordinate reference system EPSG:32610, UTM Zone 10N, (meters)
crs_m_rast = "EPSG:32610"

public_lands = st_read("data/environment/GIS Data/WA_Major_Public_Lands_(non-DNR)/WA_Major_Public_Lands_(non-DNR).shp") %>%
  janitor::clean_names() %>% st_transform(crs = crs_m) %>%
  filter(name %in% c("Olympic National Park", "Olympic National Forest"))

path_rast_cover_clean = "data/cache/occurrence_covariates/rast_cover_clean.tif"
message('Loading raster cover data from cache ', path_rast_cover_clean)
rast_cover_clean = rast(path_rast_cover_clean)
rast_cover_clean = as.factor(rast_cover_clean)

# Manually overwrite understory reinit cells as old forest -- 4 is now LSOG
values(rast_cover_clean)[values(rast_cover_clean) == 3] <- 4

mapview(rast_cover_clean, alpha.regions = 1.0,
        col.regions = c('#90c6bd', '#3c8273', '#9b652b', '#b2675e', 'darkgray', '#6495ed')) +
  mapview(study_landscape_planning_units, alpha.regions = 0.5) +
  # mapview(public_lands, alpha.regions = 0.5) +
  # mapview(watersheds, alpha.regions = 0.5) +
  mapview(aru_sites)
# + mapview(st_buffer(watershed, dist = max_range))


# Initialize results list
results_list = list()
message("Calculating proportional abundance of cover classes for landscape planning units:")
for (i in seq_len(nrow(study_landscape_planning_units))) {
  unit_name = study_landscape_planning_units$unit[i]
  message(unit_name)

  unit_poly = study_landscape_planning_units[i, ] %>% st_make_valid()
  unit_poly = st_collection_extract(unit_poly, "POLYGON") %>% st_cast("MULTIPOLYGON") %>% st_union()

  counts = terra::extract(
    rast_cover_clean, vect(unit_poly),
    fun = table, na.rm = TRUE, ID = TRUE
  )
  counts_long = as.data.frame(counts) %>%
    tidyr::pivot_longer(cols = -ID, names_to = "cover_class", values_to = "count", values_drop_na = TRUE) %>%
    mutate(unit = unit_name, prop_abund = count / sum(count)) %>% select(unit, cover_class, prop_abund)
  
  results_list[[i]] = counts_long
}
cover_results = bind_rows(results_list)

# "A conservation objective of the OESF is to support old-forest ecosystem functions, including that of spotted owl habitat, partly through providing a shifting mosaic of stands that are managed to retain or develop structural complexity... To learn to integrate older forest ecosystem values and their functions with commercial forest activities assuming, as a working hypothesis, that landscapes managed for a fairly even apportionment of forest cover among stands in all stages of development, from stand initiation to old growth, will support desirable levels of both commodities and ecosystem functions." - WADNR HCP
ggplot(cover_results, aes(x = unit, y = prop_abund, fill = cover_class)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  viridis::scale_fill_viridis(discrete = TRUE) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Landscape planning unit", y = "Proportional abundance", fill = "Cover class")

cover_results %>% filter(cover_class == 2)

# "At least 20 percent of DNR-managed lands in the landscape planning unit in the understory-reinitiation to old-growth stages" - WADNR HCP
# "This analysis shows that landscapes with less than 20% coverage by older forest rarely provide suitable habitat for northern spotted owls." - Bart and Forsman 1992
cover_results %>% filter(cover_class == 4)

# "At least 40 percent of DNR-managed lands in the landscape planning unit in the stem-exclusion to old-growth stages"
cover_results %>%
  filter(cover_class %in% c(2,4,5)) %>% group_by(unit) %>%
  summarise(prop_abund = sum(prop_abund, na.rm = TRUE))

stop()

################################################################################

pnts_name = "no_action"
data_plot_scale = readRDS(paste0("data/cache/habitat_vars/data_plot_scale_", pnts_name, ".rds"))
data_homerange_scale = readRDS(paste0("data/cache/habitat_vars/data_homerange_scale_", pnts_name, ".rds"))

# one_to_one = terra::as.points(crop(rast_cover_clean, vect(test_watershed)), na.rm = FALSE)
# nrow(one_to_one)
# mapview(one_to_one)
# 
# 
# # Get bbox and clean up landscape planning units for sampling
# bbox = st_bbox(landscape_planning_units)
# landscape_cleaner <- st_buffer(landscape_planning_units, -30) |> st_buffer(30)
# mapview(bbox) + mapview(landscape_planning_units) + mapview(landscape_cleaner)
# 
# cellsize = 250 # distance between points (meters)
# grid <- st_sf(st_make_grid(
#   x = st_as_sfc(bbox),
#   cellsize = cellsize,     
#   offset = c(bbox["xmin"], bbox["ymin"]),
#   what = "centers"    # generate points (not polygons)
# ))
# grid <- grid[st_within(grid, st_union(landscape_cleaner), sparse = FALSE), ]
# # grid <- st_crop(grid, test_watershed)
# 
# mapview(bbox) + mapview(landscape_cleaner) + mapview(grid) # + mapview(one_to_one)

# TODO: remove riparian buffers?

# convert to raster
grid$value <- runif(nrow(grid), min = 0, max = 1)
r <- rast(vect(grid), resolution = cellsize)
r <- rasterize(vect(grid), r, field = "value")

# Sauer et al. 2012 - "Summary and Prediction Using the Model  We used the model to predict occupancy, by species and for groups of species, under present habitat conditions and under 3 alternative management scenarios. We used point-specific predicted habitat data under each scenario in conjunction with estimated species-specific habitat associations from the model to predict occupancy at a point. We defined pointspecific habitat associated with a management strategy (Hm(j,h)), and used the a(i,h) to predict cmði; jÞ for each species i. From cmði; jÞ, the total occupancy over all the points in the area for management strategy m and species i is  P  j cmði; jÞ, the total summed occupancy for all (or for a  group of i < I) species is P  i  P  j cmði; jÞ. We defined this  total sum as the objective function, which is a metric that reflects the cumulative occupancy of the target species over all habitats. We note that our objective functions can be interpreted as prediction of the cumulative species richness, summed over the points, as summed (among species) predicted occupancy at a point is equivalent to estimated species richness; its maximum value is thus the number of species times the number of points, I J... We considered 4 alternative management scenarios for the Patuxent Wildlife Refuge: 1) maintaining similar amounts of habitats on the refuge; 2) allowing meadows on the refuge to revert to mixed deciduous forest (no meadows scenario); 3) allowing wetlands (primarily impoundments) to revert to mixed deciduous forest (no wetlands scenario); and 4) allowing both meadows and wetlands to revert to mixed deciduous forest (no meadows or wetlands scenario)... To model the habitat consequences of the management scenarios, we modified the habitat information associated with each point, replacing the wetland and meadow habitat at each point with mixed deciduous forest."

# Get environmental covariates for all raster cells of the map

# Generate occupancy predictions for each cell for each species

  # Compute the linear predictor using the posterior draws and back-transform

  # Average across posterior samples to get mean occupancy probability, credible intervals, and posterior sd

# Compute expected community/group richness for each cell as the sum of occupancy probability among species, including mean and credible intervals

# Compute the cumulative occupancy for each posterior draw across the entire landscape as the sum of per cell community/group richness

# Map community richness (and group richness) and uncertainty as a raster

# Map species-specific occupancy probability and uncertainty (sd) as a raster




