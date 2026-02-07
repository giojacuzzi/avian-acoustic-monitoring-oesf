# global.R ###################################################################################################
# Global packages, data, and helper functions sourced by multiple scripts 
##############################################################################################################

# Load required packages (automatically install any missing) -------------------------------------------------

if (!exists("pkgs", envir = .GlobalEnv)) {
  message("Loading required packages (automatically installing any missing)")
  pkgs = c(
    # data manipulation
    "arrow",            # cache data compression
    "janitor",          # data cleaning
    "tidyverse",        # general purpose
    "units",            # unit standardization
    # geospatial data
    "geosphere",        # distance metrics
    "leafsync",         # synchronized mapview panels
    "mapview",          # interactive geospatial visualization
    "sf",               # vector data
    "terra",            # raster data
    "tidyterra",        # raster data manipulation
    "tigris",           # political boundaries
    # visualization and plotting
    "ggplot2",
    "ggbeeswarm",       # beeswarm figures
    "ggnewscale",       # multiple scales
    "ggrepel",          # annotations
    "ggeffects",        # marginal effects
    "metR",             # contours
    "patchwork",        # multipanel plots
    "viridis",          # colors
    # statistics
    "car",              # variance inflation factors
    "DHARMa",           # modeling diagnostics
    "MuMIn",            # model selection
    "jagsUI",           # hierarchical bayesian MSOM
    "landscapemetrics", # landscape metrics
    "piecewiseSEM",     # structural equation modeling
    "PRROC",            # classifier performance evaluation
    "vegan",            # community ecology methods
    # utility
    "crayon",           # console warnings
    "progress"          # dynamic progress bar
  )
  print(sapply(pkgs, function(pkg) {
    if (!pkg %in% installed.packages()[, "Package"]) install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
    as.character(packageVersion(pkg)) # print package version
  }))
}

# Class labels ------------------------------------------------------------------------------------

if (!exists("class_labels", envir = .GlobalEnv)) {
  class_labels = readLines("data/models/ensemble/ensemble_class_labels.txt") %>% tolower() %>% tibble(label = .) %>%
    separate(label, into = c("scientific_name", "common_name"), sep = "_", extra = "merge", fill  = "right", remove = FALSE) %>%
    select(label, common_name, scientific_name)
}

# GIS data ----------------------------------------------------------------------------------------

# Coordinate reference system EPSG:32610, UTM Zone 10N, (meters)
crs_m = 32610 # coordinate reference system EPSG:32610, UTM Zone 10N, (meters)
crs_m_rast = "EPSG:32610"

# Conversion factors
conv_in_to_cm =                  2.54 # inches to centimeters
conv_ft_to_m =                 0.3048 # feet to meters
conv_ft2peracre_to_m2perha = 0.229568 # square feet/acre to square meters/ha
conv_ft3_to_m3 =            0.0283168 # cubic feet to cubic meters
conv_ft3peracre_to_m3perha = 0.069968 # cubic feet/acre to cubic meters/hectare
conv_peracre_to_perha =       2.47105 # units per acre to units per hectare
conv_perha_to_peracre =      0.404686 # units per hectare to units per acre
conv_m2_to_ha =                0.0001 # square meters to ha

# ARU locations (sites)
if (!exists("aru_sites", envir = .GlobalEnv)) {
  aru_sites = st_read('data/environment/GIS Data/AcousticStations.shp', quiet = TRUE) %>%
    st_drop_geometry() %>% janitor::clean_names() %>% as.data.frame() %>%
    st_as_sf(coords = c("utm_e", "utm_n"), crs = crs_m) %>%
    select(name, ces, treatment, geometry) %>% rename(site = name) %>%
    mutate(site = tolower(site))
}

# Table linking unique sampling unit IDs with season/serialno/deploy combinations ("unit_key.csv")
path_unit_key = "data/unit_key.csv"

# Helper functions -----------------------------------------------------------------------------

pairwise_collinearity = function(vars, threshold = 0.0) {
  cor_matrix = cor(vars, use = "pairwise.complete.obs", method = "pearson")
  cor_matrix[lower.tri(cor_matrix, diag = TRUE)] = NA
  return(collinearity_candidates = subset(as.data.frame(as.table(cor_matrix)), !is.na(Freq) & abs(Freq) >= threshold))
}

# Set ggplot theme -----------------------------------------------------------------------------------------------

# Adapted from James Robinson
theme_sleek = function(base_size = 11, base_family = "") {
  half_line = base_size/2
  theme_light(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.length = unit(half_line / 2.2, "pt"),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(colour = "grey30"),
      strip.text.y = element_text(colour = "grey30"),
      axis.text = element_text(colour = "grey30"),
      axis.title = element_text(colour = "grey30"),
      legend.title = element_text(colour = "grey30", size = rel(0.9)),
      panel.border = element_rect(fill = NA, colour = "grey70", linewidth = 1),
      legend.key.size = unit(0.9, "lines"),
      legend.text = element_text(size = rel(0.7), colour = "grey30"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA),
      plot.title = element_text(colour = "grey30", size = rel(1)),
      plot.subtitle = element_text(colour = "grey30", size = rel(.85))
    )
}
theme_set(theme_sleek())
