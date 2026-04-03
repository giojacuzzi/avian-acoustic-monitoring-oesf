# Calculate distance to nearest edge for each site

cover_classification = "clean_strata_4" # e.g. clean_stage_3, strata_4
t = 2020 # year: 2020, 2021, 2022, 2023
path_rast_cover       = paste0("data/cache/3_gis/2_gen_cover_rasters/rast_cover_", t, "_", cover_classification, ".tif")
path_site_cover_class = paste0("data/cache/3_gis/2_gen_cover_rasters/site_cover_class_", t, "_sf.rds")

message('Loading raster cover data from cache ', path_rast_cover)
rast_cover = rast(path_rast_cover)

sites = read_rds(path_site_cover_class)

mapview(rast_cover) + mapview(sites)

# Find edges in rast_cover
rast_edges = boundaries(rast_cover, inner = TRUE, classes = TRUE)
rast_edges[rast_edges == 0] = NA
# Compute a distance raster (m)
rast_dist_to_edge = distance(rast_edges)
# Get distance at each site
sites = sites %>% mutate(edge_dist_m = terra::extract(rast_dist_to_edge, vect(sites))[, 2])

# ── 4. Quick sanity check ─────────────────────────────────────────────────────
summary(sites$dist_to_edge_m)
mapview(rast_dist_to_edge) + mapview(sites, zcol = "dist_to_edge_m")
