# Crop RS-FRIS data layers to study area and cache

source("src/global.R")

path_rsfris = "/Volumes/gioj/OESF/gis"

# List all subdirectories
subdirs <- list.dirs(path_rsfris, recursive = FALSE)
# Initialize empty list to store results
data_list <- list()

for (subdir in subdirs) {
  # Extract version from directory name: "RS-FRIS 4.0 rasters"
  version <- sub("RS-FRIS ([0-9\\.]+) rasters", "\\1", basename(subdir))
  
  # List all .img files (exclude auxiliary files)
  img_files <- list.files(subdir, pattern = "\\.img$", full.names = TRUE)
  img_files <- img_files[!grepl("\\.img\\.", basename(img_files))]
  
  # Extract layer name: remove prefix and the final .img
  layers <- basename(img_files)
  layers <- sub("^RS_FRIS_", "", layers)       # remove prefix
  layers <- sub("\\.img$", "", layers)         # remove trailing .img
  
  # Create data.frame for this directory
  df_sub <- data.frame(
    version = version,
    path = img_files,
    layer = layers,
    stringsAsFactors = FALSE
  )
  
  data_list[[length(data_list) + 1]] <- df_sub
}

# Combine all directories
df <- do.call(rbind, data_list)

# View result
print(df)

# Initialize nested list
raster_list <- list()

# Loop over each row in df
for (i in seq_len(nrow(df))) {
  version <- df$version[i]
  layer <- df$layer[i]
  path <- df$path[i]
  
  message(path)
  
  # Initialize sublist for version if it doesn't exist
  if (is.null(raster_list[[version]])) {
    raster_list[[version]] <- list()
  }
  
  # Load raster and store
  raster_list[[version]][[layer]] <- load_raster(path)
}

base_out <- "data/environment/rsfris_study_area"

for (version in names(raster_list)) {
  
  version_dir <- file.path(base_out, version)
  dir.create(version_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (layer in names(raster_list[[version]])) {
    
    r <- raster_list[[version]][[layer]]
    
    out_path <- file.path(version_dir, paste0(layer, ".tif"))
    
    message(out_path)
    
    # Write raster
    terra::writeRaster(r, out_path, overwrite = TRUE)
    # If using raster package instead:
    # raster::writeRaster(r, out_path, overwrite = TRUE)
  }
}
message("Finished")

## Elevation

# Load raster
dem <- rast("/Volumes/gioj/OESF/gis/usgs_elevation/USGS_13.tif")

# Transform study area to raster CRS
sa <- st_transform(study_area, crs(dem))

# Convert sf to SpatVector
sa_vect <- vect(sa)

# Crop raster to bounding box
dem_crop <- crop(dem, sa_vect)

# Mask raster to polygon
dem_mask <- mask(dem_crop, sa_vect)

# Plot
plot(dem_mask)

out_path <- "data/environment/elevation/elevation.tif"
dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
terra::writeRaster(dem_mask, out_path, overwrite = TRUE)
