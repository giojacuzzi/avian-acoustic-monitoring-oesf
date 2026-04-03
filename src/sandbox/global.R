message("Loading global data")

library(dplyr)
library(sf)

# Coordinate reference system EPSG:32610, UTM Zone 10N, (meters)
crs_m = 32610
crs_m_rast = "EPSG:32610"

# Conversion factors
conv_in_to_cm =                  2.54 # inches to centimeters
conv_ft_to_m =                 0.3048 # feet to meters
conv_ft2peracre_to_m2perha = 0.229568 # square feet/acre to square meters/ha
conv_ft3_to_m3 = 0.0283168            # cubic feet to cubic meters
conv_ft3peracre_to_m3perha = 0.069968 # cubic feet/acre to cubic meters/hectare
conv_peracre_to_perha = 2.47105       # units per acreto units per hectare
conv_perha_to_peracre = 0.404686      # units per hectare to units per acre
conv_m2_to_ha         = 0.0001        # square meters to ha # TODO: CHECK THIS

# ARU locations (sites)
aru_sites = st_read('data/environment/GIS Data/AcousticStations.shp', quiet = TRUE) %>%
  st_drop_geometry() %>% janitor::clean_names() %>% as.data.frame() %>%
  st_as_sf(coords = c("utm_e", "utm_n"), crs = crs_m) %>%
  select(name, ces, treatment, geometry) %>% rename(site = name) %>%
  mutate(site = tolower(site))

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

# Species trait data
species_trait_data = read.csv("data/cache/species_traits/species_traits.csv")
max_range = max(species_trait_data$home_range_radius_m)

# Study area boundary buffer
study_area = st_buffer(st_as_sfc(st_bbox(landscape_planning_units)), dist = max_range)

# DEBUG
watersheds = st_read("/Users/giojacuzzi/repos/avian-acoustic-monitoring-oesf/data/environment/GIS Data/watersheds/watersheds.shp", quiet = TRUE) %>% janitor::clean_names() %>% st_transform(crs = crs_m)
test_watershed = watersheds %>% filter(exu == "Ap")
