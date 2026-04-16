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
    "MCMCvis",          # MCMC visual inspection
    "landscapemetrics", # landscape metrics
    "piecewiseSEM",     # structural equation modeling
    "PRROC",            # classifier performance evaluation
    "betapart",         # beta diversity
    "vegan",            # community ecology methods
    "indicspecies",     # indicator species
    # utility
    "benchmarkme",      # runtime benchmarking
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

# Conservation priority species -------------------------------------------------------------------

# ACAD Global - Species of Continental Importance
acad_continental = readxl::read_xlsx("data/traits/ACAD Global 2024.05.23.xlsx", sheet = 1) %>% clean_names() %>%
  mutate(common_name = str_to_lower(common_name)) %>%
  select(
    common_name,
    pt_c,
    ccs_b, # Continental Combined Score for breeding season
    continental_importance,
    iucn_red_list_2023
  )
# This includes species listed as Near Threatened, Vulnerable, and Endangered on the IUCN red list (2023) as well as common birds in steep decline whose populations have declined continentally by an estimated 50% or more since 1970, or are currently experiencing accelerating short-term decline (ACAD 2024)
conpri_con_species = acad_continental %>% filter(!is.na(iucn_red_list_2023) | !is.na(continental_importance)) %>% pull(common_name) %>% sort()

# ACAD Regional - Species of Regional Importance
acad_regional = read.csv("data/traits/ACAD Regional 2024.06.03-filtered.csv") %>% clean_names() %>%
  mutate(common_name = str_to_lower(common_name)) %>%
  select(
    common_name,
    pt_r,
    rcs_b, # Regional Combined Score for breeding season
    regional_importance,
    action_code
  )
# Critical recovery, immediate management, and management attention
conpri_reg_species = acad_regional %>% filter(action_code %in% c("CR", "IM", "MA")) %>% pull(common_name) %>% sort()

conpri_con_species %in% conpri_reg_species

# Listed at the state or federal level endangered, threatened, candidate, and species of concern
conpri_local_species = c("northern spotted owl", "marbled murrelet", "pileated woodpecker", "golden eagle", "northern goshawk", "peregrine falcon", "vaux's swift")

conpri_species = union(union(conpri_con_species, conpri_reg_species), conpri_local_species) %>% sort()

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
)
stages_4 = tibble(
  class   = c("standinit", "compex", "underdev", "old"),
  age_min = c(0, 25, 80, 200),
  age_max = c(25, 80, 200, Inf),
  idx = 1:4
)

## Management strata (e.g. Minkova)
strata_4 = tibble(
  class   = c("thin", "standinit", "compex", "mature"),
  age_min = c(NA, 0, 25, 80),
  age_max = c(NA, 25, 80, Inf),
  idx = 1:4
)
strata_5 = tibble(
  class   = c("thin", "standinit", "compex", "underdev", "old"),
  age_min = c(NA, 0, 25,  80, 200),
  age_max = c(NA, 25, 80, 200, Inf),
  idx = 1:5
)

# WADNR landscape planning units
wadnr_units = st_read("data/environment/GIS Data/WA_DNR_Units/WA_DNR_Units.shp", quiet = TRUE) %>%
  janitor::clean_names() %>% st_transform(crs = crs_m) %>%
  filter(jurisdic_2 %in% c("Willy - Huel", "Kalaloch", "Copper Mine", "Upper Clearwate", "Queets")) %>% # select units that were sampled
  mutate(jurisdic_2 = ifelse(jurisdic_2 == "Upper Clearwate", "Upper Clearwater", jurisdic_2))

wadnr_parcels = st_read("data/environment/GIS Data/WA_DNR_Managed_Land_Parcels/WA_DNR_Managed_Land_Parcels.shp", quiet = TRUE) %>%
  janitor::clean_names() %>% st_transform(crs = crs_m)
wadnr_parcels = st_filter(wadnr_parcels, st_union(wadnr_units)) # select parcels within units

landscape_planning_units = st_intersection(wadnr_units, st_union(wadnr_parcels)) %>%
  st_make_valid() %>% select(jurisdic_2) %>% rename(unit = jurisdic_2)
landscape_planning_units_clean = st_buffer(landscape_planning_units, -30) |> st_buffer(30)
# mapview(landscape_planning_units) + mapview(landscape_planning_units_clean)

# Helper functions -----------------------------------------------------------------------------

progress_bar_format = "[:bar]:percent :elapsedfull (ETA :eta)"

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

stage_colors = c(
  "Stand Initiation"      = "orange",
  "Competitive Exclusion" = "forestgreen",
  "Thinned"               = "purple",
  "Mature"                = "tan4"
)
