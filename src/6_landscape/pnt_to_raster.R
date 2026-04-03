# Convert pnts to rasters

path_optimize = "data/cache/7_landscape/OPT_landscape_data_homerange_scale_2020_clean_strata_4_landscape.rds"

data = read_rds(path_optimize)

points = st_as_sf(data[["plot"]], coords = c("x", "y"), crs = 32610) # EPSG

mapview(points)

df = cbind(st_coordinates(points), st_drop_geometry(points))

rast_from_points = rast(df[, c("X", "Y", "pcnt_mature")], type = "xyz", crs = "EPSG:32610")

mapview(rast_from_points)
