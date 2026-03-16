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
    "FD",               # functional diversity
    "ade4",             # fourth-corner statistic
    "jagsUI",           # hierarchical bayesian MSOM
    "landscapemetrics", # landscape metrics
    "piecewiseSEM",     # structural equation modeling
    "PRROC",            # classifier performance evaluation
    "betapart",         # beta diversity
    "vegan",            # community ecology methods
    "indicspecies",     # indicator species
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

class_labels = readLines("data/models/ensemble/ensemble_class_labels.txt") %>% tolower() %>% tibble(label = .) %>%
  separate(label, into = c("scientific_name", "common_name"), sep = "_", extra = "merge", fill  = "right", remove = FALSE) %>%
  select(label, common_name, scientific_name)

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

# Define study area as bounding rectangle containing max species home range buffer (6500 m) around sites
study_area = st_as_sfc(st_bbox(st_buffer(aru_sites, 6500)))

# Table linking unique sampling unit IDs with season/serialno/deploy combinations ("unit_key.csv")
path_site_key = "data/site_key.csv"
site_key = read_csv(path_site_key, show_col_types = FALSE) %>% mutate(
  site = str_to_lower(site),
  site_agg = str_to_lower(site_agg)
)

# RS-FRIS versions per year
# 4.0 - "as flown" 2019-2020
# 5.0 - "as flown" 2021-2022
# 5.1 - depletion updates for DNR lands through 9/1/2024
rsfris_version_years = data.frame(
  year =    c(2020,  2021,  2022,  2023),
  version = c("4.0", "5.0", "5.0", "5.0")
)

# Developmental stage classifications (e.g. O'Hara et al. 1996, Oliver and Larson 1996, and Spies 1997)
stages_3 = tibble(
  class   = c("standinit", "compex", "mature"),
  age_min = c(0, 25, 80),
  age_max = c(25, 80, Inf),
  idx = 1:3
); print(stages_3)
stages_4 = tibble(
  class   = c("standinit", "compex", "underdev", "old"),
  age_min = c(0, 25, 80, 200),
  age_max = c(25, 80, 200, Inf),
  idx = 1:4
); print(stages_4)

## Management strata (e.g. Minkova)
strata_4 = tibble(
  class   = c("thin", "standinit", "compex", "mature"),
  age_min = c(NA, 0, 25, 80),
  age_max = c(NA, 25, 80, Inf),
  idx = 1:4
); print(strata_4)
strata_5 = tibble(
  class   = c("thin", "standinit", "compex", "underdev", "old"),
  age_min = c(NA, 0, 25,  80, 200),
  age_max = c(NA, 25, 80, 200, Inf),
  idx = 1:5
); print(strata_5)

# Helper functions -----------------------------------------------------------------------------

pairwise_collinearity = function(vars, threshold = 0.8) {
  cor_matrix = cor(vars %>% select(where(is.numeric)), use = "pairwise.complete.obs", method = "pearson")
  cor_matrix[lower.tri(cor_matrix, diag = TRUE)] = NA
  return(collinearity_candidates = subset(as.data.frame(as.table(cor_matrix)), !is.na(Freq) & abs(Freq) >= threshold))
}

pairwise_collinearity_by_group = function(data_vars, group_var, threshold = 0.8) {
  candidates_by_group = split(data_vars %>% select(where(is.numeric)), data_vars[[group_var]])
  collinearity_results = do.call(rbind, lapply(names(candidates_by_group), function(group) {
    v = candidates_by_group[[group]]
    v = v[, apply(v, 2, sd, na.rm = TRUE) > 0, drop = FALSE] # Drop zero-variance variables in group
    if (ncol(v) < 2) return()
    cor_matrix = cor(v, use = "pairwise.complete.obs", method = "pearson")
    cor_matrix[lower.tri(cor_matrix, diag = TRUE)] = NA # Keep only upper triangle
    df = subset(as.data.frame(as.table(cor_matrix)), !is.na(Freq) & abs(Freq) >= threshold)
    if (nrow(df) > 0) {
      df[[group_var]] = group
    }
    return(df)
  }))
  return(collinearity_results)
}

# Load a raster, cropped, projected, and masked to the study area
load_raster = function(path_rast) {
  r = rast(path_rast)
  sa_r = vect(st_transform(study_area, crs(r)))
  r = crop(r, sa_r)
  r = mask(r, sa_r)
  r = project(r, crs_m_rast)
  return(r)
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
