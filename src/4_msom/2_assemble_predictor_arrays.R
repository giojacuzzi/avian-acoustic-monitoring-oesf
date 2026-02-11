# 2_assemble_predictor_arrays.R ####################################################################################
# Derive MSOM-ready data arrays:
# - TODO
#
# CONFIG:

#
# OUTPUT:
out_cache_dir  = "data/cache/4_msom/2_assemble_predictor_arrays"
# TODO
#
# INPUT:
path_xyday  = "data/cache/4_msom/1_assemble_detection_arrays/xyday.rds"
path_predictors_detection = "data/cache/detection_covariates/data_detection.rds"
##################################################################################################################

source("src/global.R")

if (!dir.exists(out_cache_dir)) dir.create(out_cache_dir, recursive = TRUE)

# Load dependencies ----------------------------------------------------------------------------------------------

message("Loading detection yday arrays from ", path_xyday)
x_yday = readRDS(path_xyday)

# TODO: Load occurrence and detection covariate data

# Discard any sites not surveyed ---------------------------------------------------------------------------------

# TODO: Discard sites with no survey observations

# Assemble all detection predictor data --------------------------------------------------------------------------

# TODO: Create arrays for all detection predictors

# TODO: Ensure detection data is available for all sites

# Assemble all occurrence predictor data -------------------------------------------------------------------------

# TODO: Create arrays for all occurrence predictors

# TODO: Ensure occurrence data is available for all sites

# Save results to cache ------------------------------------------------------------------------------------------

# TODO: Save results to cache
