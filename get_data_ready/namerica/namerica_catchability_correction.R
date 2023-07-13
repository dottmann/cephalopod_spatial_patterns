
## Script name: Add catchability to data
##
## Authors: Daniel Ottmann 
## Email: daniel.ottmann.riera@gmail.com
##
## Date created: November 2022
## Last update:  July 2023
##
## ---------------------------
##
## Readme:
##
## This script corrects the CPUE by the species catchability
## It also converts CPUE from Kg/HA to Kg/Km2
##
## ---------------------------


#####################################################################
# Clear environment:
rm(list = ls())

# Load libraries:
library(tidyverse)

# Load data:
adapt <- readRDS("data/all-regions-full-oceanadapt.rds")
load(file = "data/data_cephalopoda_hexa_id.Rdata")
load(file = "data/data_fish_hexa_id.Rdata")
load("data/oceanadapt_species_qvalues.RData")
squid_spp <- read.delim("data/cephalopoda_a_size4.txt", sep = "\t", stringsAsFactors = T)

################################################################################

# Calculate growth parameter A:
squid_spp <- squid_spp %>%
  mutate(A_growth = (3/(lifespan/12)) * asymptotic_weight^(1/3))


# Pick the species info from each observation:
data_cephalopoda_hexa_id2 <- data_cephalopoda_hexa_id %>%
  dplyr::select(spp, genus, superorder) %>%
  group_by(spp) %>%
  slice(1) %>%
  dplyr::select(-genus, -superorder)

# Merge it with the biological info from each species:
squid_spp <- squid_spp %>%
  left_join(data_cephalopoda_hexa_id2, by = "spp")

# Edit catchability table:
gd_new <- gd_new %>% 
  rename(spp = spec)

# Add squid:
gd_new <- bind_rows(gd_new, squid_spp)

# Get unique haul ID for ceophalopods and fish:
unique_cephalopoda_haulid <- data_cephalopoda_hexa_id %>%
  dplyr::select(haulid, hexa_id, coast) %>%
  group_by(haulid) %>%
  slice(1)

unique_fish_haulid <- data_fish_hexa_id %>%
  dplyr::select(haulid, hexa_id, coast) %>%
  group_by(haulid) %>%
  slice(1)

# Merge tables:
adapt <- cbind(adapt, gd_new[match(adapt$spp, gd_new$spp), c("class", "superorder", "family", "Efficiency", "asymptotic_weight", "lifespan", "A_growth")])
adapt <- subset(adapt, !(is.na(adapt$Efficiency))) # all others are invertebrates

# Convert all to kg/km2
adapt$wtcpue   <- adapt$wtcpue * 100 # kg/HA to kg/km2

# Calculate cpue corrected by catchability:
adapt$wtcpue_q <- adapt$wtcpue / adapt$Efficiency 

# Add haul id's for fish and cephalopods:
adapt_cephalopoda <- adapt %>%
  filter(class == "Cephalopoda")  %>%
  left_join(unique_cephalopoda_haulid, by = "haulid")

adapt_fish <- adapt %>%
  filter(class != "Cephalopoda")

adapt_fish <- adapt_fish %>%
  left_join(unique_fish_haulid, by = "haulid")

# Merge them again:
survey_data <- rbind(adapt_cephalopoda, adapt_fish)

save(survey_data, file = "data/namerica_survey_data.Rdata")

#                                  END OF SCRIPT
###############################################################################