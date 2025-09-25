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
aru_sites = st_read('data/environment/GIS Data/AcousticStations.shp') %>%
  st_drop_geometry() %>% janitor::clean_names() %>% as.data.frame() %>%
  st_as_sf(coords = c("utm_e", "utm_n"), crs = crs_m) %>%
  select(name, ces, treatment, geometry) %>% rename(site = name) %>%
  mutate(site = tolower(site))

# Study area boundary buffer
study_area = st_buffer(st_as_sfc(st_bbox(aru_sites)), dist = 10000)
