
## Script name: North America survey data at 400m depth cutoff
##
## Authors: Daniel Ottmann 
## Email: daniel.ottmann.riera@gmail.com
##
## Date created: October 2022
## Last update:  February 2023
##
## ---------------------------
##
## Readme:
##
## This script merges survey data with COBALT data.
## It cuts the depth off at 400 m.
## It also enables to select different catchability coefficients.
##
## ---------------------------


#####################################################################
# Clear environment:
rm(list = ls())

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

# Load data:
load(file = "data/data_cobalt_hexa_id_namerica.Rdata")
load(file = "data/data_cobalt.Rdata")
load(file = "data/namerica_survey_data.Rdata")  # Already include the traits

# Get shapes of ecoregions:
shape <- readOGR(dsn = "data/MEOW shapefiles" ,layer="meow_ecos")

# Create "exclude" function:
"%ni%" <- Negate("%in%")


#########################################################################
# Make lists of east, west coasts:
west_ecoreg <- c("Aleutian Islands", "Eastern Bering Sea", "Gulf of Alaska", "North American Pacific Fijordland", "Northern California",  
                 "Oregon, Washington, Vancouver Coast and Shelf",  "Southern California Bight")

east_ecoreg <- c("Carolinian", "Floridian", "Gulf of Maine/Bay of Fundy", "Gulf of St. Lawrence - Eastern Scotian Shelf", 
                 "Northern Grand Banks - Southern Labrador", "Southern Grand Banks - South Newfoundland", "Northern Gulf of Mexico", "Southern Gulf of Mexico", "Puget Trough/Georgia Basin", 
                 "Scotian Shelf", "Virginian")

# Group sampling stations by ecoregions in survey data:
coord <- data.frame(Longitude = survey_data$lon, Latitude = survey_data$lat)
coordinates(coord)<- ~ Longitude + Latitude 
crs(coord) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
shape <- spTransform(shape, CRS(proj4string(coord))) # make it similar to bargrid
shape@proj4string # check coordinates reference system again
tr <- over(coord, shape)
survey_data$eco_code <- tr[,1]
survey_data$eco_reg  <- tr[,2]

# All the missing ecoregions are Alutian Islands. Fix it:
survey_data <- survey_data %>%
  mutate(eco_reg = case_when(is.na(eco_reg) ~  "Aleutian Islands",
                             T ~ eco_reg),
         eco_code = case_when(is.na(eco_code) ~  20053,
                              T ~ eco_code))

range(table(survey_data$eco_reg))

rm(list = c("tr", "coord"))


# Group sampling stations by ecoregions in cobalt data:
coord <-data.frame(Longitude = data_cobalt_hexa_id_n_america$lon, Latitude = data_cobalt_hexa_id_n_america$lat)
coordinates(coord)<- ~ Longitude + Latitude
crs(coord) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
shape <- spTransform(shape, CRS(proj4string(coord))) # make it similar to bargrid
shape@proj4string # check coordinates reference system again
tr <- over(coord, shape)
data_cobalt_hexa_id_n_america$eco_code <- tr[,1]
data_cobalt_hexa_id_n_america$eco_reg  <- tr[,2]

# All the missing ecoregions are Alutian Islands. Fix it:
data_cobalt_hexa_id_n_america <- data_cobalt_hexa_id_n_america %>%
  mutate(eco_reg = case_when(is.na(eco_reg) ~  "Aleutian Islands",
                             T ~ eco_reg),
         eco_code = case_when(is.na(eco_code) ~  20053,
                              T ~ eco_code))

range(table(data_cobalt_hexa_id_n_america$eco_reg))

rm(list = c("tr", "shape", "coord"))

#########################################################################################
# Define "pelagic" fish families:
pel_family <- c("Clupeidae" , "Osmeridae",  "Exocoetidae" , "Atherinidae" , "Engraulidae",
                "Hemiramphidae", "Inermiidae","Belonidae","Scomberesocidae", "Echeneidae",
                "Carangidae","Bramidae","Scombridae","Centrolophidae","Istiophoridae","Ammodytidae")

# Areas where cephalopods were not ampled:
eco_reg_unsampled_cephalopods <- c("Eastern Bering Sea",  "Northern Grand Banks - Southern Labrador", 
                                   "Southern Grand Banks - South Newfoundland", "Gulf of St. Lawrence - Eastern Scotian Shelf")

wrong_data <-  11536

#########################################################################
# Correct for catchability:
survey_data$wtcpue_q <- survey_data$wtcpue / survey_data$Efficiency 

#--------------------------------------------------------------------------------------
# For sensitivity analysis only.
# Calculate remove catchability effect:
# survey_data$wtcpue_q <- survey_data$wtcpue 
#--------------------------------------------------------------------------------------


# Restrict to 400 m:
survey_data <- survey_data %>%
  filter(depth <= 400) %>%
  filter(eco_reg %ni% eco_reg_unsampled_cephalopods, hexa_id != wrong_data)

data_cephalopoda_hexa_id400 <- survey_data %>%
  filter(class == "Cephalopoda") 

data_fish_hexa_id400 <- survey_data %>%
  filter(class != "Cephalopoda") 

data_demersal_fish_hexa_id400 <- survey_data %>%
  filter(class != "Cephalopoda", family %ni% pel_family) 

data_cobalt_hexa_id_n_america400 <- data_cobalt_hexa_id_n_america  %>%
  filter(depth <= 400)

data_cobalt_hexa_id_n_america <- data_cobalt_hexa_id_n_america %>%
  dplyr::select(-depth)

data_cobalt400 <- data_cobalt %>%
  filter(depth <= 400)


# Take a look at the data status:
p <- ggplot() +
  geom_polygon(data = map_data("world"), aes(x = long, y = lat, group = group), fill="lightgray", colour = "white") +
  geom_point(data = data_cobalt, aes(x = long, y = lat), size = .5) +
  coord_map(xlim = c(-170, -50), ylim = c(25, 72)) +
  theme_base()  + 
  theme(plot.background = element_blank()) +
  geom_point(data = survey_data, aes(x = lon, y = lat), size = .2, alpha = .01, color = "yellow") +
  theme(legend.position = "none")

p

# Take average biomass by hexagon:
# Fish:
df_fish <- data_fish_hexa_id400 %>%
  
  # First take average of each haul, hexagon per year
  group_by(haulid, hexa_id, year) %>%
  summarise(depth = mean(depth),
    haul_wtcpue_q = sum(wtcpue_q)) %>%
  
  # Average of each haul, hexagon per year
  group_by(hexa_id, year) %>%
  summarise(depth = mean(depth),
            mean_wt_cpue_fish = mean(haul_wtcpue_q)) %>%
  
  # Average of each hexagon across years
  group_by(hexa_id) %>%
  summarise(depth = mean(depth),
            mean_wt_cpue_fish = mean(mean_wt_cpue_fish))

# All hauls have fish, but not necessarily squid. Keep HaulIDs to generate 0 squids:
haul_ids <- data_fish_hexa_id400 %>%
  group_by(haulid, hexa_id, year) %>%
  dplyr::select(haulid, hexa_id, year) %>%
  slice(1)


# Cephaopods:
df_cephalopoda <- data_cephalopoda_hexa_id400 %>%
  right_join(haul_ids, by = c("haulid", "hexa_id", "year")) %>%
  mutate(wtcpue_q = replace_na(wtcpue_q, 0))%>%
  # First take average of each haul, hexagon per year
  group_by(haulid, hexa_id, year) %>%
  summarise(haul_wtcpue_q = sum(wtcpue_q),
            asymptotic_weight = exp(sum(wtcpue_q * log(asymptotic_weight)) / sum(wtcpue_q)),
            lifespan = exp(sum(wtcpue_q * log(lifespan)) / sum(wtcpue_q)),
            A_growth = exp(sum(wtcpue_q * log(A_growth)) / sum(wtcpue_q))) %>%
  
  # Average of each haul, hexagon per year
  group_by(hexa_id, year) %>%
  summarise(mean_wtcpue_q = mean(haul_wtcpue_q),
            mean_a_size_weighted_wt_cpue_cephalopoda = mean(asymptotic_weight, na.rm = T),
            mean_lifespan_weighted_wt_cpue_cephalopoda = mean(lifespan, na.rm = T),
            mean_A_growth_weighted_wt_cpue_cephalopoda = mean(A_growth, na.rm = T)) %>%
  
  # Average of each hexagon across years
  group_by(hexa_id) %>%
  summarise(mean_wt_cpue_cephalopoda = mean(mean_wtcpue_q),
            mean_a_size_weighted_wt_cpue_cephalopoda = mean(mean_a_size_weighted_wt_cpue_cephalopoda, na.rm = T),
            mean_lifespan_weighted_wt_cpue_cephalopoda = mean(mean_lifespan_weighted_wt_cpue_cephalopoda, na.rm = T),
            mean_A_growth_weighted_wt_cpue_cephalopoda = mean(mean_A_growth_weighted_wt_cpue_cephalopoda, na.rm = T))

# Squid:
df_squid <- data_cephalopoda_hexa_id400 %>%
  filter(superorder == "decapodiformes") %>%
  right_join(haul_ids, by = c("haulid", "hexa_id", "year")) %>%
  mutate(wtcpue_q = replace_na(wtcpue_q, 0)) %>%
  
  # First take average of each haul, hexagon per year
  group_by(haulid, hexa_id, year) %>%  
  summarise(haul_wtcpue_q = sum(wtcpue_q)) %>%
  
  # Average of each haul, hexagon per year
  group_by(hexa_id, year) %>%
  summarise(#depth = mean(depth),
    mean_wt_cpue_squid = mean(haul_wtcpue_q)) %>%
  
  # Average of each hexagon across years
  group_by(hexa_id) %>%
  summarise(#depth = mean(depth),
            mean_wt_cpue_squid = mean(mean_wt_cpue_squid))


# Octopus:
df_octopus <- data_cephalopoda_hexa_id400 %>%
  filter(superorder == "octopodiformes") %>%
  right_join(haul_ids, by = c("haulid", "hexa_id", "year")) %>%
  mutate(wtcpue_q = replace_na(wtcpue_q, 0)) %>%
  
  # First take average of each haul, hexagon per year
  group_by(haulid, hexa_id, year) %>%  
  summarise(haul_wtcpue_q = sum(wtcpue_q)) %>%
  
  # Average of each haul, hexagon per year
  group_by(hexa_id, year) %>%
  summarise(#depth = mean(depth),
    mean_wt_cpue_octopus = mean(haul_wtcpue_q)) %>%
  
  # Average of each hexagon across years
  group_by(hexa_id) %>%
  summarise(#depth = mean(depth),
    mean_wt_cpue_octopus = mean(mean_wt_cpue_octopus))


# Demersals:
df_demersal_fish <- data_demersal_fish_hexa_id400 %>%
  right_join(haul_ids, by = c("haulid", "hexa_id", "year")) %>%
  mutate(wtcpue_q = replace_na(wtcpue_q, 0)) %>%
  
  # First take average of each haul, hexagon per year
  group_by(haulid, hexa_id, year, depth) %>%  
  summarise(haul_wtcpue_q = sum(wtcpue_q)) %>%
  
  # Average of each haul, hexagon per year
  group_by(hexa_id, year) %>%
  summarise(depth = mean(depth),
            mean_wt_cpue_demersal_fish = median(haul_wtcpue_q)) %>%
  
  # Average of each hexagon across years
group_by(hexa_id) %>%
  summarise(depth = mean(depth),
            mean_wt_cpue_demersal_fish = median(mean_wt_cpue_demersal_fish))

# Put them all together and get the cephalopod to fish ratio:
df_cobalt <- df_fish %>%
  left_join(data_cobalt_hexa_id_n_america) %>%
  left_join(df_demersal_fish) %>%
  left_join(df_cephalopoda) %>%
  left_join(df_squid) %>%
  left_join(df_octopus) %>%
  mutate(cephalopoda_fish_ratio = mean_wt_cpue_cephalopoda / mean_wt_cpue_fish,
         cephalopoda_fish_proportion = mean_wt_cpue_cephalopoda / (mean_wt_cpue_cephalopoda + mean_wt_cpue_fish),
         cephalopoda_demersal_fish_proportion = mean_wt_cpue_cephalopoda / (mean_wt_cpue_cephalopoda + mean_wt_cpue_demersal_fish),
         z_prod = lz_prod + mz_prod,
         coast = case_when(lon > -110 & lon < 100 ~ "East coast",
                           T ~ "West coast"))


###################
# Plot ecoregions:
p <- ggplot() +
  coord_map("gilbert", xlim = c(-190, -50), ylim = c(25, 72)) +
  geom_point(data = df_cobalt, aes(x = lon, y = lat, color = eco_reg)) +
  theme_base() + 
  theme(plot.background = element_blank()) +
  theme(legend.position = "none")
p

# Put NAs to 0:
df_cobalt[is.na(df_cobalt$mean_wt_cpue_cephalopoda), ]$mean_wt_cpue_cephalopoda <- 0
df_cobalt[is.na(df_cobalt$mean_wt_cpue_fish), ]$mean_wt_cpue_fish <- 0
df_cobalt[is.na(df_cobalt$mean_wt_cpue_squid), ]$mean_wt_cpue_squid <- 0
df_cobalt[is.na(df_cobalt$mean_wt_cpue_octopus), ]$mean_wt_cpue_octopus <- 0
df_cobalt[is.na(df_cobalt$cephalopoda_fish_ratio), ]$cephalopoda_fish_ratio <- 0
df_cobalt[is.na(df_cobalt$cephalopoda_fish_proportion), ]$cephalopoda_fish_proportion <- 0
df_cobalt[is.na(df_cobalt$cephalopoda_demersal_fish_proportion), ]$cephalopoda_demersal_fish_proportion <- 0

#########################################################
# Biomass distribution:
# In some regions, the biomass of cephalopod was not registered. 
# Put them in a vector:
eco_reg_unsampled_cephalopods <- c("Eastern Bering Sea",  "Northern Grand Banks - Southern Labrador", 
                                   "Southern Grand Banks - South Newfoundland", "Gulf of St. Lawrence - Eastern Scotian Shelf")

# Remove them:
df_cobalt <- df_cobalt %>%
  filter(eco_reg %ni% eco_reg_unsampled_cephalopods)


###########################################################################
# Prepare data set to run some stats:

stats_df <- df_cobalt %>%
  dplyr::select(hexa_id, depth, z_prod, Temp_C, lat, lon, coast, mean_wt_cpue_cephalopoda, mean_wt_cpue_fish, mean_wt_cpue_demersal_fish,
  eco_reg, mean_a_size_weighted_wt_cpue_cephalopoda, mean_lifespan_weighted_wt_cpue_cephalopoda, mean_A_growth_weighted_wt_cpue_cephalopoda,
  cephalopoda_fish_proportion, cephalopoda_demersal_fish_proportion) %>%
  mutate(coast = as.factor(coast),
         eco_reg = as.factor(eco_reg))

# save(stats_df, file = "data/us_cephalopods_stats_data_sensitivity_03q.Rdata")
# save(stats_df, file = "data/us_cephalopods_stats_data_sensitivity_noq.Rdata")


#                                  END OF SCRIPT
###############################################################################