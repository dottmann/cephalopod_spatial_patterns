
## Script name: Cobalt to ecoregions
##
## Authors: Daniel Ottmann 
## Email: daniel.ottmann.riera@gmail.com
##
## Date created: December 2022
## Last update:  July 2023
##
## ---------------------------
##
## Readme:
##
## This script allocates each cobalt data point to a given MEOW ecoregion
## Then it takes the average values of temperature, zooplankton productivity, and depth
##
## ---------------------------


#####################################################################
# Load libraries:
library(tidyverse)
library(ggthemes)
library(viridis)
library(sf)
library(data.table)
library(maptools)
library(ggmap)
library(maps)
library(raster)
library(mapproj)
library(rgdal)
library(patchwork)
library(viridis)


# Clear environment:
rm(list = ls())

# Load cobalt data:
load(file = "data/data_cobalt.Rdata")
load(file = "data/sau_predictions.Rdata")

# Load shapes of ecoregions:
shape <- readOGR(dsn = "data/meow shapefiles/MEOW shapefiles" ,layer="meow_ecos")
shape_lme <- readOGR(dsn = "data/LME shapefiles" ,layer="LME66")

# Create "exclude" function:
"%ni%" <- Negate("%in%")

#####################################################################################################

# Convert positive lon values to negative:
data_cobalt <- data_cobalt %>%
  mutate(long = case_when(long > 180 ~ (long - 360),
                          T ~ long))

# Group sampling stations by ecoregions in cobalt data:
coord <- data.frame(Longitude = data_cobalt$long, Latitude = data_cobalt$lat)
coordinates(coord)<- ~ Longitude + Latitude
crs(coord) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
shape <- spTransform(shape, CRS(proj4string(coord))) # make it similar to bargrid
shape@proj4string # check coordinates reference system again
tr <- over(coord, shape)
data_cobalt$eco_code <- tr[, 1]
data_cobalt$eco_reg  <- tr[, 2]
data_cobalt$province  <- tr[, 4]
data_cobalt$realm  <- tr[, 6]
data_cobalt$lat_zone  <- tr[, 9]

# Clean-up the environment:
rm(list = c("tr", "shape", "coord"))

# Fix the name of one ecoregion:
data_cobalt$eco_reg[data_cobalt$eco_reg == "North American Pacific Fijordland"] <- "North American Pacific Fjordland"


# Repeat for LME:
coord <- data.frame(Longitude = data_cobalt$long, Latitude = data_cobalt$lat)
coordinates(coord)<- ~ Longitude + Latitude
crs(coord) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
shape_lme <- spTransform(shape_lme, CRS(proj4string(coord))) # make it similar to bargrid
shape_lme@proj4string # check coordinates reference system again
tr <- over(coord, shape_lme)
data_cobalt$lme <- tr[, 3]

# Clean-up the environment:
rm(list = c("tr", "shape_lme", "coord"))


# Filter and select data:
df <-  data_cobalt %>%
  dplyr::select(eco_reg, eco_code, province, realm, lat_zone, lme, long, lat, depth, Temp_C, lz_prod, mz_prod, ben_prod) %>%
  mutate(z_prod = lz_prod + mz_prod) %>%
  filter(!is.na(eco_reg)) 


# LME:
lme_remove <- c("Beaufort Sea", "Canadian High Arctic - North Greenland", "Central Arctic", "Greenland Sea",
                "Kara Sea", "Laptev Sea", "Northern Bering - Chukchi Seas", "East Siberian Sea", "Antarctica")

eco_reg_lme <- df %>%
  dplyr::select(eco_reg, lme) %>%
  filter(!is.na(lme) & lme %ni% lme_remove) %>%
  group_by(eco_reg) %>%
  slice(1)

# Now group the parameters by ecoregion:
df2 <- df %>%
  filter(eco_reg %in% eco_reg_lme$eco_reg) %>%
  group_by(eco_reg, realm, province, lat_zone) %>%
  summarise(Temp_C = mean(Temp_C),
            z_prod = mean(z_prod),
            depth = mean(depth)) %>%
  left_join(eco_reg_lme, by = "eco_reg")

df <- df %>%
  dplyr::select(-depth, -z_prod, -Temp_C) %>%
  left_join(df2, by = "eco_reg")


#--------------------------------------------------------
# Plot spatial predictions:
world_map = map_data("world") %>% 
  filter(! long > 180)

# Add model predictions:
stats_df2 <- stats_df2 %>%
  dplyr::select(eco_reg, cephalopod_fish_proportion, pred_proportion, trunc_vals, model_residuals)

df <- df %>%
  left_join(stats_df2, by = "eco_reg") %>%
  filter(!is.na(pred_proportion))

df <- df %>%
  mutate(cephalopod_fish_proportion = 100 * cephalopod_fish_proportion,
         pred_proportion = 100 * pred_proportion,
         
    ceph_observ_cat = case_when(cephalopod_fish_proportion > 0 & cephalopod_fish_proportion <= 1 ~ "1",
                                     cephalopod_fish_proportion > 1 & cephalopod_fish_proportion <= 2 ~ "2",
                                     cephalopod_fish_proportion > 2 & cephalopod_fish_proportion <= 5 ~ "3",
                                     cephalopod_fish_proportion > 5 & cephalopod_fish_proportion <= 10 ~ "4",
                                     cephalopod_fish_proportion > 10 & cephalopod_fish_proportion <= 20 ~ "5",
                                     cephalopod_fish_proportion > 20 ~ "6",
                                   T ~ "0"),
    
         pred_prop_cat = case_when(pred_proportion > 0 & pred_proportion <= 1 ~ "1",
                                   pred_proportion > 1 & pred_proportion <= 2 ~ "2",
                                   pred_proportion > 2 & pred_proportion <= 5 ~ "3",
                                   pred_proportion > 5 & pred_proportion <= 10 ~ "4",
                                   pred_proportion > 10 & pred_proportion <= 20 ~ "5",
                                   pred_proportion > 20 ~ "6",
                                   T ~ "0"))


p1 <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "gray60", colour = "white", size = .2)  +
  geom_point(data = df, aes(x = long, y = lat, color = ceph_observ_cat), size = .2) +
  scale_color_viridis_d(name = "Observations %",
                        labels = c("0", "0-1", "1-2", "2-5", "5-10", "10-20", ">20")) +
  coord_map("moll") +
  theme_map() +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(size = 4)))

p1

p2 <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "gray60", colour = "white", size = .2)  +
  geom_point(data = df, aes(x = long, y = lat, color = pred_prop_cat), size = .2) +
  scale_color_viridis_d(begin = 0,
                        end = 4/7,
                        name = "Predictions %",
                        labels = c("0", "0-1", "1-2", "2-5", "5-10", "10-20", ">20")) +
  coord_map("moll") +
  theme_map()  +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(size = 4)))

p2

df <- df %>%
  mutate(resids_trunc = case_when(model_residuals > .2 ~ .2,
                                  T ~ model_residuals))

p3 <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "gray60", colour = "white", size = .2)  +
  geom_point(data = df, aes(x = long, y = lat, color = 100 * model_residuals), size = .2) +
  scale_color_gradient2(low = "blue", mid = 'lightyellow', high = 'red',
                        name = "Difference %") +
  coord_map("moll") +
  theme_map() +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.position = "right")

p3

p4 <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "gray60", colour = "white", size = .2)  +
  geom_point(data = df, aes(x = long, y = lat, color = 100 * resids_trunc), size = .2) +
  scale_color_gradient2(low = "blue", mid = 'lightyellow', high = 'red',
                        limits = 100 * c(-.07, .2),
                        name = "Difference %") +
  coord_map("moll") +
  theme_map() +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.position = "right")

p4

p <- (p1 /p2 / p4) +
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 10)) & 
  theme(plot.background = element_blank())

ggsave("plots/ecoregions_world_combined.png", p, height = 95 , width = 70, units = "mm", scale = 3)


p <- (p1 / p4) +
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 10)) & 
  theme(plot.background = element_blank())


ggsave("plots/ecoregions_world_combined2.png", p, height = 55 , width = 70, units = "mm", scale = 2)
