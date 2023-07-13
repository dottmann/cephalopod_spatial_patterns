
## Script name: Spatial plotting of model predictions
##
## Authors: Daniel Ottmann 
## Email: daniel.ottmann.riera@gmail.com
##
## Date created: March 2023
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
library(patchwork)
library(sf)
library(data.table)
library(maptools)
library(ggmap)
library(maps)
library(raster)
library(mapproj)
library(geosphere)
library(rnaturalearth)
library(viridis)

# Clear environment:
rm(list = ls())

# Load data:
load(file = "data/predictions_biomass_q03.Rdata")
load(file = "data/predictions_traits_q03.Rdata")

biomass_df <- stats_df3
traits_df <- stats_df4

rm(list = c("stats_df3", "stats_df4"))

ctrys <- rnaturalearth::ne_countries(scale = 50, returnclass = "sf")

# Create "exclude" function:
"%ni%" <- Negate("%in%")

biomass_df %>%
  group_by(coast) %>%
  summarise(mean_ceph = mean(mean_wt_cpue_cephalopoda),
            mean_prop = 100 * mean(cephalopoda_fish_proportion),
            mean_a_size = mean(mean_a_size_weighted_wt_cpue_cephalopoda, na.rm = T),
            mean_a_growth = mean(mean_A_growth_weighted_wt_cpue_cephalopoda, na.rm = T),
            
            min_ceph = min(mean_wt_cpue_cephalopoda),
            min_prop = 100 * min(cephalopoda_fish_proportion),
            min_a_size = min(mean_a_size_weighted_wt_cpue_cephalopoda, na.rm = T),
            min_a_growth = min(mean_A_growth_weighted_wt_cpue_cephalopoda, na.rm = T),
            
            max_ceph = max(mean_wt_cpue_cephalopoda),
            max_prop = 100 * max(cephalopoda_fish_proportion),
            max_a_size = max(mean_a_size_weighted_wt_cpue_cephalopoda, na.rm = T),
            max_a_growth = max(mean_A_growth_weighted_wt_cpue_cephalopoda, na.rm = T))

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


#############################################
# Edit data:
biomass_df <- biomass_df %>%
  mutate(diffs_proportion = cephalopoda_fish_proportion - pred_proportion,
         diffs_ceph_biomass = mean_wt_cpue_cephalopoda - pred_ceph_biomass,
         diffs_fish_biomass = mean_wt_cpue_fish - pred_fish_biomass) %>%
  dplyr::select(-mean_a_size_weighted_wt_cpue_cephalopoda, -mean_lifespan_weighted_wt_cpue_cephalopoda, -mean_A_growth_weighted_wt_cpue_cephalopoda)

traits_df <- traits_df %>%
  mutate(pred_a_size = 10^pred_a_size,
         diffs_a_size = mean_a_size_weighted_wt_cpue_cephalopoda - pred_a_size,
         diffs_lifespan = mean_lifespan_weighted_wt_cpue_cephalopoda - pred_lifespan,
         diffs_a_growth = mean_A_growth_weighted_wt_cpue_cephalopoda - pred_a_growth) %>%
  dplyr::select(hexa_id, mean_a_size_weighted_wt_cpue_cephalopoda, mean_lifespan_weighted_wt_cpue_cephalopoda, 
                mean_A_growth_weighted_wt_cpue_cephalopoda, 
                pred_a_size, pred_lifespan, pred_a_growth, diffs_a_size, diffs_lifespan, diffs_a_growth)

biomass_df <- biomass_df %>%
  left_join(traits_df, by = "hexa_id")

biomass_df_eu <- biomass_df %>%
  filter(coast == "Europe") %>%
  mutate(
    cephalopoda_fish_proportion = 100 * cephalopoda_fish_proportion,
    diffs_proportion = 100 * diffs_proportion,
    pred_proportion = 100 * pred_proportion)


biomass_df_usa <- biomass_df %>%
  filter(coast %in% c("East coast", "West coast")) %>%
  mutate(
    cephalopoda_fish_proportion = 100 * cephalopoda_fish_proportion,
    diffs_proportion = 100 * diffs_proportion,
    pred_proportion = 100 * pred_proportion)



# Fit survey coordenates within exagons:
grid_eu <- hexa[unique(biomass_df_eu$hexa_id)]
grid_usa <- hexa[unique(biomass_df_usa$hexa_id)]



#---------------------------------------------------------------------------------------------

# Combine multiple plots in one:

# Set the same limits for all regions:
vars <- c(6, 7, 10, 13:27)
legend_title <- names(as.data.frame(biomass_df)[vars])

upper_limit <- c(max(biomass_df$mean_wt_cpue_cephalopoda, na.rm = T),
                 max(biomass_df$mean_wt_cpue_fish, na.rm = T),
                 max(biomass_df$cephalopoda_fish_proportion, na.rm = T),
                 
                 max(biomass_df$cephalopoda_fish_proportion, na.rm = T),
                 max(biomass_df$mean_wt_cpue_cephalopoda, na.rm = T),
                 max(biomass_df$mean_wt_cpue_fish, na.rm = T),
                 
                 max(biomass_df$cephalopoda_fish_proportion, na.rm = T),
                 max(biomass_df$diffs_ceph_biomass, na.rm = T),
                 max(biomass_df$diffs_fish_biomass, na.rm = T),
                 
                 max(biomass_df$mean_a_size_weighted_wt_cpue_cephalopoda, na.rm = T),
                 max(biomass_df$mean_lifespan_weighted_wt_cpue_cephalopoda, na.rm = T),
                 max(biomass_df$mean_A_growth_weighted_wt_cpue_cephalopoda, na.rm = T),
                 
                 max(biomass_df$mean_a_size_weighted_wt_cpue_cephalopoda, na.rm = T),
                 max(biomass_df$mean_lifespan_weighted_wt_cpue_cephalopoda, na.rm = T),
                 max(biomass_df$mean_A_growth_weighted_wt_cpue_cephalopoda, na.rm = T),
                 
                 max(biomass_df$mean_a_size_weighted_wt_cpue_cephalopoda, na.rm = T),
                 max(biomass_df$mean_lifespan_weighted_wt_cpue_cephalopoda, na.rm = T),
                 max(biomass_df$mean_A_growth_weighted_wt_cpue_cephalopoda, na.rm = T))

lower_limit <- c(min(biomass_df$mean_wt_cpue_cephalopoda, na.rm = T),
                 min(subset(biomass_df, mean_wt_cpue_fish > 0)$mean_wt_cpue_fish, na.rm = T),
                 min(biomass_df$cephalopoda_fish_proportion, na.rm = T),
                 
                 min(biomass_df$cephalopoda_fish_proportion, na.rm = T),
                 min(biomass_df$mean_wt_cpue_cephalopoda, na.rm = T),
                 min(subset(biomass_df, mean_wt_cpue_fish > 0)$mean_wt_cpue_fish, na.rm = T),
                 
                 min(biomass_df$diffs_proportion, na.rm = T),
                 min(biomass_df$diffs_ceph_biomass, na.rm = T),
                 min(biomass_df$diffs_fish_biomass, na.rm = T),
                 
                 min(biomass_df$mean_a_size_weighted_wt_cpue_cephalopoda, na.rm = T),
                 min(biomass_df$mean_lifespan_weighted_wt_cpue_cephalopoda, na.rm = T),
                 min(biomass_df$mean_A_growth_weighted_wt_cpue_cephalopoda, na.rm = T),
                 
                 min(biomass_df$mean_a_size_weighted_wt_cpue_cephalopoda, na.rm = T),
                 min(biomass_df$mean_lifespan_weighted_wt_cpue_cephalopoda, na.rm = T),
                 min(biomass_df$mean_A_growth_weighted_wt_cpue_cephalopoda, na.rm = T),
                 
                 min(biomass_df$diffs_a_size, na.rm = T),
                 min(biomass_df$diffs_lifespan, na.rm = T),
                 min(biomass_df$diffs_a_growth, na.rm = T))


# Europe:
df <- st_sf(biomass_df_eu, grid_eu)

# Categorize cephalopod bimoass and %:
df <- df %>%
  mutate(ceph_biom_cat = case_when(mean_wt_cpue_cephalopoda > 0 & mean_wt_cpue_cephalopoda <= 50 ~ "1",
                                   mean_wt_cpue_cephalopoda > 50 & mean_wt_cpue_cephalopoda <= 100 ~ "2",
                                   mean_wt_cpue_cephalopoda > 100 & mean_wt_cpue_cephalopoda <= 250 ~ "3",
                                   mean_wt_cpue_cephalopoda > 250 & mean_wt_cpue_cephalopoda <= 500 ~ "4",
                                   mean_wt_cpue_cephalopoda > 500 & mean_wt_cpue_cephalopoda <= 1000 ~ "5",
                                   mean_wt_cpue_cephalopoda > 1000 ~ "6",
                                   T ~ "0"),
         ceph_prop_cat = case_when(cephalopoda_fish_proportion > 0 & cephalopoda_fish_proportion <= .25 ~ "1",
                                   cephalopoda_fish_proportion > .25 & cephalopoda_fish_proportion <= .5 ~ "2",
                                   cephalopoda_fish_proportion > .5 & cephalopoda_fish_proportion <= 1 ~ "3",
                                   cephalopoda_fish_proportion > 1 & cephalopoda_fish_proportion <= 2.5 ~ "4",
                                   cephalopoda_fish_proportion > 2.5 & cephalopoda_fish_proportion <= 5 ~ "5",
                                   cephalopoda_fish_proportion > 5 ~ "6",
                                   T ~ "0"))


# Ceph biomass:
i <- 1
my_title <- legend_title[i]
l_limit <- lower_limit[i]
u_limit <- upper_limit[i]


pa <- ggplot() +
  geom_sf(data = df, aes(fill = ceph_biom_cat)) +
  scale_fill_viridis_d(begin = 0,
                       end = 6/7,
                       labels = c("0", "0-50", "50-100","100-250", "250-500","500-1000", ">1000")) +
  # scale_fill_gradient2(low = 'blue', mid = 'lightyellow', high = 'red',
  #                      limits = c(0, 6),
  #                      labels = c("0", "0-100", "100-500","500-1000", "1000-2000","2000-5000", ">5000")) +
  labs(fill = element_blank()) +
  geom_sf(data = ctrys, fill = "grey", colour = NA) +
  # geom_polygon(data = map_data("world"), aes(x = long, y = lat, group = group), fill = "lightgray", colour = "white") + 
  coord_sf(crs =  "+proj=aea +lat_0=23 +lon_0=10 +lat_1=29.5 +lat_2=45.5 +x_0=0 
           +y_0=0 +datum=NAD83 +units=m +no_defs") +
  ylim(c(1500000, 5900000)) +
  xlim(c(-1700000, 2100000)) +
  theme_base() +
  ggtitle("Cephalopod biomass") +
  theme(plot.background = element_blank(),
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),  
        axis.ticks.y=element_blank(),
        plot.title = element_text(size = 15),
        legend.position = "none") 


# Ceph proportion:
i <- 3
my_title <- legend_title[i]
l_limit <- 100 * min(biomass_df$cephalopoda_fish_proportion, na.rm = T)
u_limit <- 100 * max(biomass_df$cephalopoda_fish_proportion, na.rm = T)


pb <- ggplot() +
  geom_sf(data = df, aes(fill = ceph_prop_cat)) +
  scale_fill_viridis_d(begin = 0,
                       end = 6/7,
                       labels = c("0", "0-0.25", "0.25-0.5", "0.5-1", "1-2.5", "2.5-5", ">5")) +
  # scale_fill_gradient2(low = 'blue', mid = 'lightyellow', high = 'red',
  #                      limits = c(0, 6),
                       # labels = c("0", "0-1", "1-2", "2-5", "5-10", "10-15", ">15")) +
  labs(fill = element_blank()) +
  geom_sf(data = ctrys, fill = "grey", colour = NA) +
  # geom_polygon(data = map_data("world"), aes(x = long, y = lat, group = group), fill = "lightgray", colour = "white") + 
  coord_sf(crs =  "+proj=aea +lat_0=23 +lon_0=10 +lat_1=29.5 +lat_2=45.5 +x_0=0 
           +y_0=0 +datum=NAD83 +units=m +no_defs") +
  ylim(c(1500000, 5900000)) +
  xlim(c(-1700000, 2100000)) +
  theme_base() +
  ggtitle("% cephalopod biomass")  +
  theme(plot.background = element_blank(),
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),  
        axis.ticks.y=element_blank(),
        plot.title = element_text(size = 15),
        legend.position = "none")


# A size:
i <- 10
my_title <- legend_title[i]
l_limit <- min(traits_df$mean_a_size_weighted_wt_cpue_cephalopoda, na.rm = T)
u_limit <- max(traits_df$mean_a_size_weighted_wt_cpue_cephalopoda, na.rm = T)
my_breaks <- c(100, 1000, 10000)

pc <- ggplot() +
  geom_sf(data = df, aes(fill = mean_a_size_weighted_wt_cpue_cephalopoda)) +
  scale_fill_viridis(trans = "log",
                     breaks = my_breaks,
                     limits = c(l_limit, u_limit)) +

  labs(fill = element_blank()) +
  geom_sf(data = ctrys, fill = "grey", colour = NA) +
  # geom_polygon(data = map_data("world"), aes(x = long, y = lat, group = group), fill = "lightgray", colour = "white") + 
  coord_sf(crs =  "+proj=aea +lat_0=23 +lon_0=10 +lat_1=29.5 +lat_2=45.5 +x_0=0 
           +y_0=0 +datum=NAD83 +units=m +no_defs") +
  ylim(c(1500000, 5900000)) +
  xlim(c(-1700000, 2100000)) +
  theme_base()  +
  ggtitle("Asymptotic size") +
  theme(plot.background = element_blank(),
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),  
        axis.ticks.y=element_blank(),
        plot.title = element_text(size = 15),
        legend.position = "none") 


# A growth:
i <- 12
my_title <- legend_title[i]
l_limit <- min(traits_df$mean_A_growth_weighted_wt_cpue_cephalopoda, na.rm = T)
u_limit <- max(biomass_df$mean_A_growth_weighted_wt_cpue_cephalopoda, na.rm = T)

pd <- ggplot() +
  geom_sf(data = df, aes(fill = mean_A_growth_weighted_wt_cpue_cephalopoda)) +
  # scale_fill_gradient(low = 'lightyellow', high = 'red', 
  #                     limits = c(l_limit, u_limit)) +
  scale_fill_viridis(limits = c(l_limit, u_limit)) +
  labs(fill = element_blank()) +
  geom_sf(data = ctrys, fill = "grey", colour = NA) +
  # geom_polygon(data = map_data("world"), aes(x = long, y = lat, group = group), fill = "lightgray", colour = "white") + 
  coord_sf(crs =  "+proj=aea +lat_0=23 +lon_0=10 +lat_1=29.5 +lat_2=45.5 +x_0=0 
           +y_0=0 +datum=NAD83 +units=m +no_defs") +
  ylim(c(1500000, 5900000)) +
  xlim(c(-1700000, 2100000)) +
  theme_base() +
  ggtitle("Growth parameter A") +
  theme(plot.background = element_blank(),
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),  
        axis.ticks.y=element_blank(),
        plot.title = element_text(size = 15),
        legend.position = "none")



# USA:
df <- st_sf(biomass_df_usa, grid_usa)

# Categorize cephalopod bimoass and %:
df <- df %>%
  mutate(ceph_biom_cat = case_when(mean_wt_cpue_cephalopoda > 0 & mean_wt_cpue_cephalopoda <= 50 ~ "1",
                                   mean_wt_cpue_cephalopoda > 50 & mean_wt_cpue_cephalopoda <= 100 ~ "2",
                                   mean_wt_cpue_cephalopoda > 100 & mean_wt_cpue_cephalopoda <= 250 ~ "3",
                                   mean_wt_cpue_cephalopoda > 250 & mean_wt_cpue_cephalopoda <= 500 ~ "4",
                                   mean_wt_cpue_cephalopoda > 500 & mean_wt_cpue_cephalopoda <= 1000 ~ "5",
                                   mean_wt_cpue_cephalopoda > 1000 ~ "6",
                                   T ~ "0"),
         ceph_prop_cat = case_when(cephalopoda_fish_proportion > 0 & cephalopoda_fish_proportion <= .25 ~ "1",
                                   cephalopoda_fish_proportion > .25 & cephalopoda_fish_proportion <= .5 ~ "2",
                                   cephalopoda_fish_proportion > .5 & cephalopoda_fish_proportion <= 1 ~ "3",
                                   cephalopoda_fish_proportion > 1 & cephalopoda_fish_proportion <= 2.5 ~ "4",
                                   cephalopoda_fish_proportion > 2.5 & cephalopoda_fish_proportion <= 5 ~ "5",
                                   cephalopoda_fish_proportion > 5 ~ "6",
                                   T ~ "0"))
# Ceph biomass:
i <- 1
my_title <- legend_title[i]
l_limit <- lower_limit[i]
u_limit <- upper_limit[i]


pe <- ggplot() +
  geom_sf(data = df, aes(fill = ceph_biom_cat)) +
  # scale_fill_gradient2(low = 'blue', mid = 'lightyellow', high = 'red', 
  #                      limits = c(0, 6),
  #                      labels = c("0", "0-100", "100-500","500-1000", "1000-2000","2000-5000", ">5000")) +
  scale_fill_viridis_d(labels = c("0", "0-50", "50-100","100-250", "250-500","500-1000", ">1000")) +
  labs(fill = element_blank()) +
  geom_sf(data = ctrys, fill = "grey", colour = NA) +
  # geom_polygon(data = map_data("world"), aes(x = long, y = lat, group = group), fill = "lightgray", colour = "white") + 
  coord_sf(crs = 5070, xlim = c(-5200000, 2600000), ylim = c(307076, 6200000)) +
  # ylim(c(4000000, 7200000)) +
  # xlim(c(-10000000, -5000000)) +
  theme_base() +
  theme(plot.background = element_blank(),
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),  
        axis.ticks.y=element_blank()) +
  labs(fill = expression(paste("kg ", km^-2, " ", yr^-1)))



# Ceph proportion:
i <- 3
my_title <- legend_title[i]
l_limit <- 100 * min(biomass_df$cephalopoda_fish_proportion, na.rm = T)
u_limit <- 100 * max(biomass_df$cephalopoda_fish_proportion, na.rm = T)


pf <- ggplot() +
  geom_sf(data = df, aes(fill = ceph_prop_cat)) +
  # scale_fill_gradient2(low = 'blue', mid = 'lightyellow', high = 'red',
  #                      limits = c(0, 6),
  #                      labels = c("0", "0-1", "1-2", "2-5", "5-10", "10-15", ">15")) +
  scale_fill_viridis_d(labels = c("0", "0-0.25", "0.25-0.5", "0.5-1", "1-2.5", "2.5-5", ">5")) +
  labs(fill = element_blank()) +
  geom_sf(data = ctrys, fill = "grey", colour = NA) +
  # geom_polygon(data = map_data("world"), aes(x = long, y = lat, group = group), fill = "lightgray", colour = "white") + 
  coord_sf(crs = 5070, xlim = c(-5200000, 2600000), ylim = c(307076, 6200000)) +
  # ylim(c(4000000, 7200000)) +
  # xlim(c(-10000000, -5000000)) +
  theme_base() +
  theme(plot.background = element_blank(),
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),  
        axis.ticks.y=element_blank()) +
  labs(fill = "%")



# A size:
i <- 10
my_title <- legend_title[i]
l_limit <- min(traits_df$mean_a_size_weighted_wt_cpue_cephalopoda, na.rm = T)
u_limit <- max(traits_df$mean_a_size_weighted_wt_cpue_cephalopoda, na.rm = T)
my_breaks <- c(100, 1000, 10000)

pg <- ggplot() +
  geom_sf(data = df, aes(fill = mean_a_size_weighted_wt_cpue_cephalopoda)) +
  scale_fill_viridis(trans = "log",
                     breaks = my_breaks,
                     limits = c(l_limit, u_limit)) +
  labs(fill = element_blank()) +
  geom_sf(data = ctrys, fill = "grey", colour = NA) +
  # geom_polygon(data = map_data("world"), aes(x = long, y = lat, group = group), fill = "lightgray", colour = "white") + 
  coord_sf(crs = 5070, xlim = c(-5200000, 2600000), ylim = c(307076, 6200000)) +
  # ylim(c(4000000, 7200000)) +
  # xlim(c(-10000000, -5000000)) +
  theme_base() +
  theme(plot.background = element_blank(),
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),  
        axis.ticks.y=element_blank()) +
  labs(fill = "g")


# A growth:
i <- 12
my_title <- legend_title[i]
l_limit <- min(traits_df$mean_A_growth_weighted_wt_cpue_cephalopoda, na.rm = T)
u_limit <- max(traits_df$mean_A_growth_weighted_wt_cpue_cephalopoda, na.rm = T)

ph <- ggplot() +
  geom_sf(data = df, aes(fill = mean_A_growth_weighted_wt_cpue_cephalopoda)) +
  # scale_fill_gradient(low = 'lightyellow', high = 'red', 
  #                     limits = c(l_limit, u_limit)) +
  scale_fill_viridis(limits = c(l_limit, u_limit)) +
  labs(fill = element_blank()) +
  geom_sf(data = ctrys, fill = "grey", colour = NA) +
  # geom_polygon(data = map_data("world"), aes(x = long, y = lat, group = group), fill = "lightgray", colour = "white") + 
  coord_sf(crs = 5070, xlim = c(-5200000, 2600000), ylim = c(307076, 6200000)) +
  # ylim(c(4000000, 7200000)) +
  # xlim(c(-10000000, -5000000)) +
  theme_base() +
  theme(plot.background = element_blank(),
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),  
        axis.ticks.y=element_blank()) +
  labs(fill = expression(paste(g^-(1/3), yr^-1)))


# Combine all panels in a plot:
p <- (pa | pe) /
  (pb | pf) /
  (pc | pg) /
  (pd | ph) +
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 15)) & 
  theme(plot.background = element_blank())

ggsave("plots/prediction_maps/Figure_2.png", p, height = 90 , width = 53, units = "mm", scale = 5)

rm(list = c("pa", "pb", "pc", "pd", "pe", "pf", "pg", "ph"))


#---------------------------------------------------------------------------------------------

# Plot model differences:

# Europe:
df <- st_sf(biomass_df_eu, grid_eu)

df <- df %>%
  mutate(diffs_ceph_biomass = case_when(diffs_ceph_biomass > 2000 ~ 2000,
                                        diffs_ceph_biomass < -2000 ~ -2000,
                                        T ~ diffs_ceph_biomass),
         diffs_proportion = case_when(diffs_proportion > 8 ~ 8,
                                      T ~ diffs_proportion),
         diffs_a_size = case_when(diffs_a_size > 10000 ~ 10000,
                                  T ~ diffs_a_size))



# Ceph biomass:
my_title <- "diffs_ceph_biomass"
l_limit <- -2000
u_limit <- 2000

pa <- ggplot() +
  geom_sf(data = df, aes(fill = diffs_ceph_biomass)) +
  scale_fill_gradient2(low = 'blue', mid = 'lightyellow', high = 'red',
                       limits = c(l_limit, u_limit)) +
  labs(fill = element_blank()) +
  geom_sf(data = ctrys, fill = "grey", colour = NA) +
  # geom_polygon(data = map_data("world"), aes(x = long, y = lat, group = group), fill = "lightgray", colour = "white") +
  coord_sf(crs =  "+proj=aea +lat_0=23 +lon_0=10 +lat_1=29.5 +lat_2=45.5 +x_0=0
           +y_0=0 +datum=NAD83 +units=m +no_defs") +
  ylim(c(1500000, 5900000)) +
  xlim(c(-1700000, 2100000)) +
  theme_base() +
  ggtitle("Cephalopod biomass") +
  theme(plot.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(size = 15),
        legend.position = "none")



# Ceph proportion:
my_title <- "diffs_proportion"
l_limit <- 100 * min(biomass_df$diffs_proportion, na.rm = T)
u_limit <- 8

pb <- ggplot() +
  geom_sf(data = df, aes(fill = diffs_proportion)) +
  scale_fill_gradient2(low = 'blue', mid = 'lightyellow', high = 'red',
                       limits = c(l_limit, u_limit)) +
  labs(fill = element_blank()) +
  geom_sf(data = ctrys, fill = "grey", colour = NA) +
  # geom_polygon(data = map_data("world"), aes(x = long, y = lat, group = group), fill = "lightgray", colour = "white") +
  coord_sf(crs =  "+proj=aea +lat_0=23 +lon_0=10 +lat_1=29.5 +lat_2=45.5 +x_0=0
           +y_0=0 +datum=NAD83 +units=m +no_defs") +
  ylim(c(1500000, 5900000)) +
  xlim(c(-1700000, 2100000)) +
  theme_base() +
  ggtitle("% cephalopod biomass")  +
  theme(plot.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(size = 15),
        legend.position = "none")



# A size:
my_title <- "diffs_a_size"
l_limit <- min(traits_df$diffs_a_size, na.rm = T)
u_limit <- 10000
my_breaks <- c(1, 10, 100, 1000, 10000)

pc <- ggplot() +
  geom_sf(data = df, aes(fill = diffs_a_size)) +
  scale_fill_gradient2(low = 'blue', mid = 'lightyellow', high = 'red', 
                      # trans = "log",
                      # breaks = my_breaks,
                      limits = c(l_limit, u_limit)) +
  labs(fill = element_blank()) +
  geom_sf(data = ctrys, fill = "grey", colour = NA) +
  # geom_polygon(data = map_data("world"), aes(x = long, y = lat, group = group), fill = "lightgray", colour = "white") + 
  coord_sf(crs =  "+proj=aea +lat_0=23 +lon_0=10 +lat_1=29.5 +lat_2=45.5 +x_0=0 
           +y_0=0 +datum=NAD83 +units=m +no_defs") +
  ylim(c(1500000, 5900000)) +
  xlim(c(-1700000, 2100000)) +
  theme_base()  +
  ggtitle("Asymptotic size") +
  theme(plot.background = element_blank(),
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),  
        axis.ticks.y=element_blank(),
        plot.title = element_text(size = 15),
        legend.position = "none") 


# A growth:
my_title <- "diffs_a_growth"
l_limit <- min(traits_df$diffs_a_growth, na.rm = T)
u_limit <- max(traits_df$diffs_a_growth, na.rm = T)

pd <- ggplot() +
  geom_sf(data = df, aes(fill = diffs_a_growth)) +
  scale_fill_gradient2(low = 'blue', mid = 'lightyellow', high = 'red', 
                      limits = c(l_limit, u_limit)) +
  labs(fill = element_blank()) +
  geom_sf(data = ctrys, fill = "grey", colour = NA) +
  # geom_polygon(data = map_data("world"), aes(x = long, y = lat, group = group), fill = "lightgray", colour = "white") + 
  coord_sf(crs =  "+proj=aea +lat_0=23 +lon_0=10 +lat_1=29.5 +lat_2=45.5 +x_0=0 
           +y_0=0 +datum=NAD83 +units=m +no_defs") +
  ylim(c(1500000, 5900000)) +
  xlim(c(-1700000, 2100000)) +
  theme_base() +
  ggtitle("Growth parameter A") +
  theme(plot.background = element_blank(),
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),  
        axis.ticks.y=element_blank(),
        plot.title = element_text(size = 15),
        legend.position = "none")



# USA:
df <- st_sf(biomass_df_usa, grid_usa)

df <- df %>%
  mutate(diffs_ceph_biomass = case_when(diffs_ceph_biomass > 2000 ~ 2000,
                                        diffs_ceph_biomass < -2000 ~ -2000,
                                        T ~diffs_ceph_biomass),
         diffs_proportion = case_when(diffs_proportion > 8 ~ 8,
                                      T ~ diffs_proportion),
         diffs_a_size = case_when(diffs_a_size > 10000 ~ 10000,
                                  T ~ diffs_a_size))

# Ceph biomass:
my_title <- "diffs_ceph_biomass"
l_limit <- -2000
u_limit <- 2000

pe <- ggplot() +
  geom_sf(data = df, aes(fill = diffs_ceph_biomass)) +
  scale_fill_gradient2(low = 'blue', mid = 'lightyellow', high = 'red',
                       limits = c(l_limit, u_limit)) +
  labs(fill = element_blank()) +
  geom_sf(data = ctrys, fill = "grey", colour = NA) +
  coord_sf(crs = 5070, xlim = c(-5200000, 2600000), ylim = c(307076, 6200000)) +
  ylim(c(1500000, 5900000)) +
  xlim(c(-1700000, 2100000)) +
  theme_base() +
  ggtitle("Cephalopod biomass") +
  theme(plot.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(size = 15))



# Ceph proportion:
my_title <- "diffs_proportion"
l_limit <- 100 * min(biomass_df$diffs_proportion, na.rm = T)
u_limit <- 8

pf <- ggplot() +
  geom_sf(data = df, aes(fill = diffs_proportion)) +
  scale_fill_gradient2(low = 'blue', mid = 'lightyellow', high = 'red',
                       limits = c(l_limit, u_limit)) +
  labs(fill = element_blank()) +
  geom_sf(data = ctrys, fill = "grey", colour = NA) +
  coord_sf(crs = 5070, xlim = c(-5200000, 2600000), ylim = c(307076, 6200000)) +
  ylim(c(1500000, 5900000)) +
  xlim(c(-1700000, 2100000)) +
  theme_base() +
  ggtitle("% cephalopod biomass")  +
  theme(plot.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(size = 15))



# A size:
my_title <- "diffs_a_size"
l_limit <- min(df$diffs_a_size, na.rm = T)
u_limit <- 10000
my_breaks <- c(1, 10, 100, 1000, 10000)

pg <- ggplot() +
  geom_sf(data = df, aes(fill = diffs_a_size)) +
  scale_fill_gradient2(low = 'blue', mid = 'lightyellow', high = 'red', 
                       # trans = "log",
                       # breaks = my_breaks,
                       limits = c(l_limit, u_limit)) +
  labs(fill = element_blank()) +
  geom_sf(data = ctrys, fill = "grey", colour = NA) +
  # geom_polygon(data = map_data("world"), aes(x = long, y = lat, group = group), fill = "lightgray", colour = "white") + 
  coord_sf(crs = 5070, xlim = c(-5200000, 2600000), ylim = c(307076, 6200000)) +
  ylim(c(1500000, 5900000)) +
  xlim(c(-1700000, 2100000)) +
  theme_base()  +
  ggtitle("Asymptotic size") +
  theme(plot.background = element_blank(),
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),  
        axis.ticks.y=element_blank(),
        plot.title = element_text(size = 15)) 


# A growth:
my_title <- "diffs_a_growth"
l_limit <- min(traits_df$diffs_a_growth, na.rm = T)
u_limit <- max(traits_df$diffs_a_growth, na.rm = T)

ph <- ggplot() +
  geom_sf(data = df, aes(fill = diffs_a_growth)) +
  scale_fill_gradient2(low = 'blue', mid = 'lightyellow', high = 'red', 
                       limits = c(l_limit, u_limit)) +
  labs(fill = element_blank()) +
  geom_sf(data = ctrys, fill = "grey", colour = NA) +
  coord_sf(crs = 5070, xlim = c(-5200000, 2600000), ylim = c(307076, 6200000)) +
  ylim(c(1500000, 5900000)) +
  xlim(c(-1700000, 2100000)) +
  theme_base() +
  ggtitle("Growth parameter A") +
  theme(plot.background = element_blank(),
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),  
        axis.ticks.y=element_blank(),
        plot.title = element_text(size = 15))



# Combine all panels in a plot:
p <- (pa | pe) /
  (pb | pf) /
  (pc | pg) /
  (pd | ph) +
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 15)) & 
  theme(plot.background = element_blank())

ggsave("plots/prediction_maps/diffs.png", p, height = 90 , width = 53, units = "mm", scale = 5)

rm(list = c("pa", "pb", "pc", "pd", "pe", "pf", "pg", "ph"))

#                          END OF SCRIPT
####################################################################