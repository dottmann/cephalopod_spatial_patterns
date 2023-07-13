
## Script name: ICES survey data at 400m depth cutoff
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
## This script merges survey data with COBALT data and with traits data.
## It cuts the depth off at 400 m.
## It also enables to select different catchability coefficients.
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

# Clear environment:
rm(list = ls())

# Load biomass and traits data:
load(file = "data/data_cobalt_hexa_ices.Rdata")
load(file = "data/ices_dat_hexa.Rdata")
cephalopod_traits <- read.delim('data/ices_cephalopoda_traits2.txt', sep = '\t', header = T, stringsAsFactors = F)


# Get shapes of ecoregions:
shape <- readOGR(dsn = "data/MEOW shapefiles" ,layer="meow_ecos")

# Create "exclude" function:
"%ni%" <- Negate("%in%")


#########################################################################
# Edit survey data:
survey_data <- survey3 %>%
  rename(depth = Depth,
         spp = Taxon) %>%
  filter(lon > -40)

# Edit cephalopod traits data:
cephalopod_traits <- cephalopod_traits %>%
  rename(spp = scientificname) %>%
  dplyr::select(spp, asymptotic_weight, asympototic_length, lifespan)  %>%
  mutate(A_growth = (3/(lifespan/12)) * asymptotic_weight^(1/3)) # Calculate growth parameter A

# Group sampling stations by ecoregions in survey data:
coord <- data.frame(Longitude = survey_data$lon, Latitude = survey_data$lat)
coordinates(coord) <- ~ Longitude + Latitude 
crs(coord) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
shape <- spTransform(shape, CRS(proj4string(coord))) # make it similar to bargrid
shape@proj4string # check coordinates reference system again
tr <- over(coord, shape)
survey_data$eco_code <- tr[,1]
survey_data$eco_reg  <- tr[,2]
table(survey_data$eco_reg)
survey_data <- subset(survey_data, !is.na(eco_reg))

rm(list = c("tr", "coord"))

# Group sampling stations by ecoregions in cobalt data:
# First select relevant hexa_ids:
data_cobalt_europe <- data_cobalt %>%
  filter(hexa_id %in% unique(survey_data$hexa_id)) %>%
  # rename(lon = long) %>%
  mutate(lon = case_when(lon >180 ~ lon - 360,
                         T ~ lon))

coord <-data.frame(Longitude = data_cobalt_europe$lon, Latitude = data_cobalt_europe$lat)
coordinates(coord)<- ~ Longitude + Latitude
crs(coord) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
shape <- spTransform(shape, CRS(proj4string(coord))) # make it similar to bargrid
shape@proj4string # check coordinates reference system again
tr <- over(coord, shape)
data_cobalt_europe$eco_code <- tr[,1]
data_cobalt_europe$eco_reg  <- tr[,2]
data_cobalt_europe <- subset(data_cobalt_europe, !is.na(eco_reg))
table(data_cobalt_europe$eco_reg)

rm(list = c("tr", "shape", "coord"))

data_cobalt <- data_cobalt %>%
  # rename(lon = long) %>%
  mutate(lon = case_when(lon >180 ~ lon - 360,
                         T ~ lon))

#--------------------------------------------------------------------------------------

# Catchability correction:

# For catchability coefficient = 0.3:
# q_efficiency <- 0.3
# survey_data <- survey_data %>%
#   mutate(wgtlencpue_q = case_when(class == "Cephalopoda" ~ wgtlencpue/q_efficiency,
#                                   T ~ wgtlencpue_q))

# For no catchability coefficient:
# survey_data <- survey_data %>%
# mutate(wgtlencpue_q = wgtlencpue)


# For size-speciffic catchability coefficient:
# Do nothing

#--------------------------------------------------------------------------------------

#########################################################################################
# Define "pelagic" fish families:
pel_family <- c("Clupeidae" , "Osmeridae",  "Exocoetidae" , "Atherinidae" , "Engraulidae",
                "Hemiramphidae", "Inermiidae","Belonidae","Scomberesocidae", "Echeneidae",
                "Carangidae","Bramidae","Scombridae","Centrolophidae","Istiophoridae","Ammodytidae")


#########################################################################
# Restrict to 400 m:
survey_data <- survey_data %>%
  filter(depth <= 400)


# Create df for each taxa:
data_cephalopoda_hexa_id400 <- survey_data %>%
  filter(class == "Cephalopoda") 

data_fish_hexa_id400 <- survey_data %>%
  filter(class != "Cephalopoda")  %>%
  filter(!is.na(wgtlencpue_q))

data_demersal_fish_hexa_id400 <- survey_data %>%
  filter(class != "Cephalopoda", Family %ni% pel_family) 

data_cobalt_europe400 <- data_cobalt_europe  %>%
  filter(depth <= 400)

data_cobalt_europe <- data_cobalt_europe %>%
  dplyr::select(-depth)

data_cobalt400 <- data_cobalt %>%
  filter(depth <= 400)

# Take a look at the data status:
p <- ggplot() +
  geom_polygon(data = map_data("world"), aes(x = long, y = lat, group = group), fill="lightgray", colour = "white") + 
  geom_point(data = data_cobalt_europe, aes(x = lon, y = lat), color = "grey") +
  geom_point(data = data_cobalt_europe400, aes(x = lon, y = lat), color = "red") +
  geom_point(data = data_cobalt, aes(x = lon, y = lat), size = .2) +
  geom_point(data = data_cobalt400, aes(x = lon, y = lat), size = .2, color = "blue") +
  geom_point(data = survey_data, aes(x = lon, y = lat), size = .2, alpha = .01, color = "yellow") +  
  coord_map("gilbert", xlim = c(-40, 40), ylim = c(25, 72)) +
  theme_bw() 

p

# Take average biomass by hexagon:
# Fish:
df_fish <- data_fish_hexa_id400 %>%
  # Add cpue by haul
  
  group_by(HaulID, hexa_id, Year, depth) %>%
  summarise(haul_wgtlencpue_q = sum(wgtlencpue_q)) %>%
  
  # Average of each haul, hexagon per year
  group_by(hexa_id, Year) %>%
  summarise(depth = mean(depth),
            wgtlencpue_q = mean(haul_wgtlencpue_q)) %>%
  
  # Average of each hexagon across years
  group_by(hexa_id) %>%
  summarise(depth = mean(depth),
            mean_wt_cpue_fish = mean(wgtlencpue_q))

# All hauls have fish, but not necessarily squid. Keep HaulIDs to generate 0 squids:
haul_ids <- data_fish_hexa_id400 %>%
  group_by(HaulID, hexa_id, Year) %>%
  dplyr::select(HaulID, hexa_id, Year) %>%
  slice(1)


# Cephalopoda:
df_cephalopoda <- data_cephalopoda_hexa_id400 %>%
  left_join(cephalopod_traits, by = "spp") %>%
  right_join(haul_ids, by = c("HaulID", "hexa_id", "Year")) %>%
  mutate(wgtlencpue_q = replace_na(wgtlencpue_q, 0)) %>%
  
  # Add cpue by haul
  group_by(HaulID, hexa_id, Year) %>%  # First take average of each haul, hexagon per year
  summarise(haul_wgtlencpue_q = sum(wgtlencpue_q),
            asymptotic_weight = exp(sum(wgtlencpue_q * log(asymptotic_weight)) / sum(wgtlencpue_q)),
            lifespan = exp(sum(wgtlencpue_q * log(lifespan)) / sum(wgtlencpue_q)),
            A_growth = exp(sum(wgtlencpue_q * log(A_growth)) / sum(wgtlencpue_q))) %>%
  
  #  # Average of each haul, hexagon per year
  group_by(hexa_id, Year) %>%
  summarise(mean_wgtlencpue_q = mean(haul_wgtlencpue_q),
            mean_a_size_weighted_wt_cpue_cephalopoda = mean(asymptotic_weight, na.rm = T),
            mean_lifespan_weighted_wt_cpue_cephalopoda = mean(lifespan, na.rm = T),
            mean_A_growth_weighted_wt_cpue_cephalopoda = mean(A_growth, na.rm = T)) %>%
  
  # Average of each hexagon across years
  group_by(hexa_id) %>%  
  summarise(mean_wt_cpue_cephalopoda = mean(mean_wgtlencpue_q),
            mean_a_size_weighted_wt_cpue_cephalopoda = mean(mean_a_size_weighted_wt_cpue_cephalopoda, na.rm = T),
            mean_lifespan_weighted_wt_cpue_cephalopoda = mean(mean_lifespan_weighted_wt_cpue_cephalopoda, na.rm = T),
            mean_A_growth_weighted_wt_cpue_cephalopoda = mean(mean_A_growth_weighted_wt_cpue_cephalopoda, na.rm = T))


# Demersals:
df_demersal_fish <- data_demersal_fish_hexa_id400 %>%
  
  # Add cpue by haul
  group_by(HaulID, hexa_id, Year, depth) %>%
  summarise(haul_wgtlencpue_q = sum(wgtlencpue_q, na.rm = T)) %>%
  
  # Average of each haul, hexagon per year
  group_by(hexa_id, Year) %>%
  summarise(depth = mean(depth),
            wgtlencpue_q = mean(haul_wgtlencpue_q, na.rm = T)) %>%
  
  # Average of each hexagon across years
  group_by(hexa_id) %>%
  summarise(depth = mean(depth),
            mean_wt_cpue_demersal_fish = mean(wgtlencpue_q))

# Put them all together and get the cephalopod to fish ratio:
df_cobalt <- df_fish %>%
  left_join(data_cobalt_europe) %>%
  left_join(df_demersal_fish) %>%
  left_join(df_cephalopoda) %>%
  mutate(cephalopoda_fish_ratio = mean_wt_cpue_cephalopoda / mean_wt_cpue_fish,
         cephalopoda_fish_proportion = mean_wt_cpue_cephalopoda / (mean_wt_cpue_cephalopoda + mean_wt_cpue_fish),
         cephalopoda_demersal_fish_proportion = mean_wt_cpue_cephalopoda / (mean_wt_cpue_cephalopoda + mean_wt_cpue_demersal_fish),
         z_prod = lz_prod + mz_prod) %>%
  filter(!is.na(eco_reg), eco_reg %ni% c("Faroe Plateau", "Baltic Sea"))

# Put NAs to 0:
df_cobalt[is.na(df_cobalt$mean_wt_cpue_cephalopoda), ]$mean_wt_cpue_cephalopoda <- 0
df_cobalt[is.na(df_cobalt$mean_wt_cpue_fish), ]$mean_wt_cpue_fish <- 0
df_cobalt[is.na(df_cobalt$cephalopoda_fish_ratio), ]$cephalopoda_fish_ratio <- 0
df_cobalt[is.na(df_cobalt$cephalopoda_fish_proportion), ]$cephalopoda_fish_proportion <- 0
df_cobalt[is.na(df_cobalt$cephalopoda_demersal_fish_proportion), ]$cephalopoda_demersal_fish_proportion <- 0

###########################################################################
# Prepare data set to run some stats:

stats_df <- df_cobalt %>%
  dplyr::select(hexa_id, depth, z_prod, Temp_C, lat, lon, mean_wt_cpue_cephalopoda, mean_wt_cpue_fish, mean_wt_cpue_demersal_fish,
  eco_reg, mean_a_size_weighted_wt_cpue_cephalopoda, mean_lifespan_weighted_wt_cpue_cephalopoda, mean_A_growth_weighted_wt_cpue_cephalopoda,
  cephalopoda_fish_proportion, cephalopoda_demersal_fish_proportion) %>%
  mutate(eco_reg = as.factor(eco_reg))

# Save according to which catchability coefficient is used:

# save(stats_df, file = "data/ices_cephalopods_stats_data_sensitivity_03q.Rdata")
# save(stats_df, file = "data/ices_cephalopods_stats_data_sensitivity_noq.Rdata")
# save(stats_df, file = "data/ices_cephalopods_stats_data_sensitivity_size_q.Rdata")


#                                 END OF SCRIPT
##############################################################################