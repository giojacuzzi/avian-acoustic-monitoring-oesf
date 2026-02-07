## 1_preprocess_gis_data.R #########################################################################################
# Load and preprocess raw GIS data
#
# OUTPUT: (various sf objects in local memory)
#
## INPUT:
path_watercourses = "data/environment/DNR_Hydrography/DNR_Hydrography_-_Watercourses_-_Forest_Practices_Regulation/DNR_Hydrography_-_Watercourses_-_Forest_Practices_Regulation.shp"
path_waterbodies  = "data/environment/DNR_Hydrography/DNR_Hydrography_-_Water_Bodies_-_Forest_Practices_Regulation/DNR_Hydrography_-_Water_Bodies_-_Forest_Practices_Regulation.shp"
path_roads_wadnr  = "data/environment/GIS Data/roads/T3roads.shp"
path_roads_wsdot  = "data/environment/WSDOT_-_Local_Agency_Public_Road_Routes/WSDOT_-_Local_Agency_Public_Road_Routes.shp"
###########################################################################################################

source("src/global.R")

message("Loading watercourses from ", path_watercourses)
# Watercourses/waterbodies (rivers and streams) from Type 1-3. These support distinct riparian vegetation which constitutes different habitat. Type 4 (non-fish-bearing streams) and type 5 are not boundaries because they are small features that do not change the vegetation, support less aquatic biota, and often are not permanent.
watercourses = st_read(path_watercourses, quiet = TRUE)
watercourses = watercourses %>% st_crop(st_transform(study_area, st_crs(watercourses))) %>%
  st_transform(crs_m) %>% janitor::clean_names() %>% select(geometry, sl_wtrty_c)
# mapview(watercourses, zcol = 'sl_wtrty_c')

boundary_watercourses = watercourses %>% filter(sl_wtrty_c %in% c(1, 2, 3))
boundary_watercourses$sl_wtrty_c = 1

message("Loading waterbodies from ", path_waterbodies)
waterbodies = st_read(path_waterbodies, quiet = TRUE)
waterbodies = waterbodies %>% st_crop(st_transform(study_area, st_crs(waterbodies))) %>%
  st_transform(crs_m) %>% janitor::clean_names()
# mapview(waterbodies, zcol = 'sl_wtrty_c')

boundary_waterbodies = waterbodies %>% filter(sl_wtrty_c %in% c(1))

# Paved roads. Unpaved roads do not constitute a boundary of a patch because they are narrow, rarely driven and likely experienced by the birds as gaps in the forest.
message("Loading WADNR roads ", path_roads_wadnr)
roads_wadnr = st_read(path_roads_wadnr, quiet = TRUE) %>%
  st_zm(drop = TRUE) %>% janitor::clean_names() %>% st_transform(crs = crs_m) %>%
  select(road_usgs1, road_surfa, geometry) %>%
  mutate(geometry = st_cast(geometry, "MULTILINESTRING"))

message("Loading WSDOT roads from ", path_roads_wsdot)
roads_wsdot = st_read(path_roads_wsdot, quiet = TRUE) %>%
  filter(RouteIdent %in% c("400000220i", "031265969i")) %>%
  st_transform(crs_m) %>% select(geometry) # get Hoh Mainline Road / Clearwater Road from WSDOT
roads_wsdot$road_usgs1 = 'Primary Highway'
roads_wsdot$road_surfa = NA

roads = rbind(roads_wsdot, roads_wadnr)
# mapview(roads, zcol = 'road_usgs1') + mapview(aru_sites, label = aru_sites$site)

paved_roads = st_make_valid(roads %>% filter(road_usgs1 %in% c("Primary Highway", "Light-Duty Road")))

