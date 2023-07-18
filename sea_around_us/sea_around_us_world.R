## Script name: Sea around us world
##
## Authors: Daniel Ottmann 
## Email: daniel.ottmann.riera@gmail.com
##
## Date created: November 2022
## Last update:  April 2023
##
## ---------------------------
##
## Readme:
##
## This script makes an overview of the data and prepares a data files to run stats.
##
## ---------------------------



#####################################################################
# Clear environment:
rm(list = ls())

# Load libraries:
library(tidyverse)
library(ggthemes)

# Load data:
load(file = "data/data_sau_world.Rdata")
load(file = "data/ecoregion_cobalt_variables.Rdata")

# Create "exclude" function:
"%ni%" <- Negate("%in%") 

#######################################################################
# Create data selection criteria:
cephalopod_and_fish <- c("Cephalopods", "Large bathydemersals (>=90 cm)", "Large benthopelagics (>=90 cm)", "Large demersals (>=90 cm)",
                         "Large flatfishes (>=90 cm)", "Large pelagics (>=90 cm)", "Large reef assoc. fish (>=90 cm)", 
                         "Medium bathydemersals (30 - 89 cm)",  "Medium benthopelagics (30 - 89 cm)", "Medium demersals (30 - 89 cm)", 
                         "Medium pelagics (30 - 89 cm)", "Medium reef assoc. fish (30 - 89 cm)", "Small pelagics (<30 cm)", 
                         "Small to medium flatfishes (<90 cm)", "Large bathypelagics (>=90 cm)", "Medium bathypelagics (30 - 89 cm)", 
                         "Small bathypelagics (<30 cm)", "Small benthopelagics (<30 cm)", "Small demersals (<30 cm)", 
                         "Small reef assoc. fish (<30 cm)", "Large sharks (>=90 cm)", "Large rays (>=90 cm)", "Small to medium rays (<90 cm)",
                         "Small to medium sharks (<90 cm)") 

cephalopod_and_demersals <- c("Cephalopods", "Large bathydemersals (>=90 cm)", "Large benthopelagics (>=90 cm)", "Large demersals (>=90 cm)",
                              "Large flatfishes (>=90 cm)", "Large reef assoc. fish (>=90 cm)", "Medium bathydemersals (30 - 89 cm)",  
                              "Medium benthopelagics (30 - 89 cm)", "Medium demersals (30 - 89 cm)", "Medium reef assoc. fish (30 - 89 cm)", 
                              "Small to medium flatfishes (<90 cm)", "Small benthopelagics (<30 cm)", "Small demersals (<30 cm)", "Small reef assoc. fish (<30 cm)",  
                              "Large rays (>=90 cm)", "Small to medium rays (<90 cm)") 

time_window <- (max(unique(data_sau_world$year)) - 20):max(unique(data_sau_world$year))

unique(data_sau_world$area_name)


# Modify data:
data_sau_1 <- data_sau_world %>% 
  filter(functional_group %in% cephalopod_and_fish & year %in% time_window) %>%
  mutate(class = case_when(functional_group == "Cephalopods" ~ "Cephalopod",
                           T ~ "Fish")) %>%
  group_by(area_name, class, year) %>%
  summarise(me_biomass = sum(tonnes)) %>%
  group_by(area_name, class) %>%
  summarise(me_biomass = mean(me_biomass))

# Split in 2 data frames for fish and cephalopod:
data_sau_fish <- data_sau_1 %>%
  filter(class == "Fish") %>%
  rename(fish_biomass = me_biomass) %>%
  dplyr::select(-class)

data_sau_cephalopod <- data_sau_1 %>%
  filter(class == "Cephalopod") %>%
  rename(cephalopod_biomass = me_biomass) %>%
  dplyr::select(-class)

# Join both data frames back by marine ecosystem and add marine ecosystem info:
me_biomass <- data_sau_fish %>% 
  left_join(data_sau_cephalopod, by = "area_name") %>% 
  mutate(cephalopod_fish_proportion = cephalopod_biomass / (fish_biomass + cephalopod_biomass)) 

# Plot it out:
p <- ggplot() +
  geom_histogram(data = me_biomass, aes(x = cephalopod_fish_proportion), binwidth = 0.01) +
  theme_base()

p

# Pick top ecosystems with highest proportion of cephalopods:
top <- me_biomass  %>%
  ungroup() %>%
  top_n(20, cephalopod_fish_proportion)


p <- ggplot(data = top) +
  geom_bar(aes(x = reorder(area_name, cephalopod_fish_proportion), y = cephalopod_fish_proportion), stat = "identity") + 
  coord_flip() +
  theme_base() +
  theme(axis.title.y = element_blank())

p


names(me_biomass)
names(ecoregion_cobalt_df)

ecoregion_cobalt_df <- ecoregion_cobalt_df %>%
  mutate(area_name = eco_reg) %>%
  dplyr::select(area_name, province, realm, lat_zone, lme, Temp_C, z_prod, depth)

me_biomass <- me_biomass %>%
  left_join(ecoregion_cobalt_df, by = "area_name")


# Remove ecosystems of the polar regions:
lme_remove <- c("Beaufort Sea", "Canadian High Arctic - North Greenland", "Central Arctic", "Greenland Sea",
                "Kara Sea", "Laptev Sea", "Northern Bering - Chukchi Seas", "East Siberian Sea", "Antarctica")

me_biomass2 <- me_biomass %>%
  filter(!is.na(lme), lme %ni% lme_remove)


sau_data_stats <- me_biomass2
save(sau_data_stats, file = "data/sau_data_stats.Rdata")

#                  END OF SCRIPT
#####################################################