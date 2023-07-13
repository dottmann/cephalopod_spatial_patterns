
## Script name: Combine ices survey data with cobalt & hexagons
##
## Authors: Daniel Ottmann 
## Email: daniel.ottmann.riera@gmail.com
##
## Date created: January 2023
## Last update:  July 2023
##
## ---------------------------
##
## Readme:
##
## This script links cobalt and survey data with a hexagon.
## Through the hexa_ID, it also links survey data with a cobalt point
##
## ---------------------------


#####################################################################
# Load libraries:
library(tidyverse)
library(ggthemes)
library(sf)
library(data.table)
library(maptools)
library(ggmap)
library(maps)
library(raster)
library(mapproj)
library(geosphere)


# Clear environment:
rm(list = ls())

# Load data:
load(file = "data/hexa_centroids.Rdata")
load(file = "data/data_cobalt.Rdata")
load(file = "data/ICESsurveys30Nov_withq_ceph.RData")


#################################################################

# Create equal size hexagons:
sfc <- st_sfc(st_polygon(list(rbind(c(-180,-50), c(54,-50), c(54,90), c(-180,90), c(-180,-50)))))
st_crs(sfc) = 4326
sfc <- st_transform(sfc, crs = st_crs(5070))
cellarea <- 6200 * (1e+6)
cellsize <- 2 * sqrt(cellarea / ((3 * sqrt(3)/2))) * sqrt(3)/2
hexa <- st_make_grid(x = sfc, cellsize = cellsize, square = FALSE)
hexa <- st_transform(hexa, crs = st_crs(4326))
hexa <- st_make_valid(hexa)

# Get centroids of each hexagon:
hexa_centroids <- st_centroid(hexa)
hexa_centroids <- lapply(hexa_centroids, as.matrix)
hexa_centroids <- lapply(hexa_centroids, as.data.frame)
hexa_centroids <- rbindlist((hexa_centroids))
hexa_centroids <- hexa_centroids %>%
  rename(lon = V1,
         lat = V2) %>%
  mutate(hexa_id = 1:nrow(hexa_centroids))


#############################################
# Fit survey coordenates within exagons:

# Get coordenates from each survey sample:
survey3_latlon <- survey3 %>%
  rename(lat = ShootLat,
         lon = ShootLong)

# Convert coordenates into geometric points:
coords <- survey3_latlon %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# Overlay grid with trawl coordinates:
newgrid <- st_intersects(hexa, coords)

# Convert to data frame and edit:
overlap <- as.data.frame(newgrid)

overlap <- overlap %>%
  rename(hexa_id = row.id,
         survey3_latlon_id = col.id) %>%
  arrange(survey3_latlon_id)

# Add hexa_id to each survey entry:
hexa_id <- overlap$hexa_id
survey3 <- cbind(survey3, hexa_id)
  
# Plot to see overlap between centroids and hexagons that have survey data:
grid_master <- hexa[unique(overlap$hexa_id)]

p <- ggplot() +
  geom_polygon(data = map_data("world"), aes(x = long, y = lat, group = group), fill="lightgray", colour = "white") + 
  geom_sf(data = grid_master) +
  geom_point(data = hexa_centroids, aes(x = lon, y = lat), size = .5, alpha = .5) +
  ylim(c(30, 70)) +
  xlim(c(-80, 60)) +
  theme_bw() 

p


# Some hexagons have no cobalt point. Others have more than one.
# Get the closest cobalt point to each hexa centroid:

# First generate a table with coordenates and id of cobalt:
cobalt_id <- data_cobalt %>%
  dplyr::select(lat = lat,
                long = long) %>%
  mutate(long = case_when(long > 180 ~ (long - 360),
                          T ~ long),
         cobalt_id = row_number())

# Remove id's because geosphere::distGeo crushes with more than 2 columns:
cobalt_xy <- cobalt_id %>%
  dplyr::select(long, lat) %>%
  as.data.frame()

# Same for hexa centroids:
hexa_centroids_xy <- hexa_centroids %>%
  dplyr::select(lon, lat) %>% 
  na.omit() %>%
  as.data.frame()


# Run a loop to get closest values.
# It takes some time, so we have run and stored it beforehand:

load(file = "data/hexa_centroids_xy.Rdata")

# # Run a loop to get closest values (takes some time):
# for(i in 1:nrow(hexa_centroids_xy)){
#   distances <- geosphere::distGeo(hexa_centroids_xy[i,], cobalt_xy)/1000
#   ranking <- rank(distances, ties.method = "first")   # rank the calculated distances
#   hexa_centroids_xy$cobalt_id[i] <- which.min(ranking) # Find shortest distance
# 
#   if(i %in% seq(from = 0, to = nrow(hexa_centroids_xy), by = 500)){  # Show % progress
#     print(paste0(round(100 * i/nrow(hexa_centroids_xy), 1), " %"))
#   }
# }


# Add cobalt_id to hexa_centroids df:
hexa_centroids <- hexa_centroids %>%
  left_join(hexa_centroids_xy, by = c("lat", "lon")) 

hexa_centroids1 <- hexa_centroids %>%
  dplyr::select(-lon, -lat)

# Link id's in survey data frame:
survey3 <- survey3 %>% left_join(hexa_centroids1)

# Filter variables of interest in data_cobalt:
data_cobalt <- data_cobalt %>%
  mutate(cobalt_id = row_number()) %>%
  dplyr::select(cobalt_id, lz_prod, mz_prod, ben_prod, depth, Temp_C)


# Create new table with cobalt data from each hexa centroid:

hexa_centroids <- hexa_centroids %>%
  left_join(data_cobalt, by = "cobalt_id")


data_cobalt <- hexa_centroids

# Rename lat lon labels:
survey3 <- survey3 %>%
  rename(lat = ShootLat,
         lon = ShootLong)

# Outfile:
# save(survey3, file = "data/survey3_hexa.Rdata")
# save(data_cobalt, file = "data/data_cobalt_hexa_europe.Rdata")


#                          END OF SCRIPT
####################################################################