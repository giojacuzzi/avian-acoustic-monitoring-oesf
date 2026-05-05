source("src/global.R")

# rast_cover = rast("data/cache/3_gis/2_gen_cover_rasters/rast_cover_2020_clean_strata_4.tif")
rast_cover = rast("data/cache/3_gis/2_gen_cover_rasters/rast_cover_2020_clean_strata_5.tif")

{
  # Union the planning units and convert to SpatVector
  lpu_union <- st_union(landscape_planning_units) |> vect()
  
  # Mask raster to the unioned planning units
  rast_masked <- mask(crop(rast_cover, lpu_union), lpu_union)
  
  # Get frequency table
  freq_table <- freq(rast_masked)
  
  # Calculate percentages
  freq_table$percent <- freq_table$count / sum(freq_table$count) * 100
  
  print(freq_table, digits = 2)
}