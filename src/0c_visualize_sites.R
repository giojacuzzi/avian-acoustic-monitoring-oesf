####################################################################################
# Visualize site locations on a map and derive GIS environmental data
#
# INPUT:
# Raw gpx data containing the coordinates for each site
path_sites_gpx = "data/sites/PAMlocations20230510.gpx"
####################################################################################

library(mapview)
library(tidyverse)
library(elevatr)
library(sf)

gpx_data = gpx::read_gpx(path_sites_gpx)$waypoints %>% janitor::clean_names() %>% rename(site = name) %>% select(site, latitude, longitude)

site_data = st_as_sf(gpx_data, coords = c("longitude", "latitude"), crs = 4269)
site_data = get_elev_point(locations = site_data, units = "meters", src = "aws")

mapview(site_data, zcol="elevation", grid=F, legend=T)
