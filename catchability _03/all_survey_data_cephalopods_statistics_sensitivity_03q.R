
## All survey data statistics with catchability coefficient = 0.3
##
## Authors: Daniel Ottmann 
## Email: daniel.ottmann.riera@gmail.com
##
## Date created: February 2023
## Last update:  July 2023
##
## ---------------------------
##
## Readme:
##
## This script evalautes with statistical means (GAMMs and Ranfom Forests) the effect that environmental variables have
## on cephalopod biomass and traits distribution.
## It uses biomass estimates of cephalopods calculated with a cathcability coefficient = 0.3 and no spatial dimension. 
##
## ---------------------------


#####################################################################
# Load libraries:
library(tidyverse)
library(ggthemes)
library(patchwork)
library(viridis)
library(mgcv)
library(lme4)
library(lattice)
library(sjstats)
library(gratia)
library(car)
library(randomForest)
# library(MixRF)
library(pdp)
library(scales)

# Clear environment:
rm(list = ls())

# Create "exclude" function:
"%ni%" <- Negate("%in%")

# Load data:
load(file = "data/ices_cephalopods_stats_data_sensitivity_03q.Rdata")
stats_df1 <- stats_df %>%
  mutate(coast = "Europe")

load(file = "data/norway_cephalopods_stats_data_sensitivity_03q.Rdata") 
stats_df2 <- stats_df %>%
  mutate(coast = "Europe")

load(file = "data/us_cephalopods_stats_data_sensitivity_03q.Rdata")
stats_df3 <- stats_df %>%
  relocate(coast, .after = last_col())

temp <- stats_df3 %>% filter(mean_wt_cpue_cephalopoda == max(stats_df3$mean_wt_cpue_cephalopoda))


#########################################################
# Merge tables:
stats_df <- rbind(stats_df1, stats_df2, stats_df3)
stats_df$coast <- as.factor(stats_df$coast)

rm(list = c("stats_df1", "stats_df2", "stats_df3"))

# In some regions, the biomass of cephalopod was not registered. 
# Put them in a vector:
eco_reg_unsampled_cephalopods <- c("Eastern Bering Sea",  "Northern Grand Banks - Southern Labrador", 
                                   "Southern Grand Banks - South Newfoundland", "Gulf of St. Lawrence - Eastern Scotian Shelf")

# We only have 1 data point from the southern GOM. Let's remove it:
eco_reg_undersampled <- "Southern Gulf of Mexico"

stats_df <- stats_df %>%
  filter(eco_reg %ni% eco_reg_unsampled_cephalopods, eco_reg %ni% eco_reg_undersampled)

stats_df2 <- stats_df %>%
  filter(eco_reg %ni% eco_reg_unsampled_cephalopods) %>%
  group_by(eco_reg, coast) %>%
  summarise(n = n(),
            
            mean_wt_cpue_cephalopoda_me = mean(mean_wt_cpue_cephalopoda, na.rm = T),
            sd_wt_cpue_cephalopoda_me = sd(mean_wt_cpue_cephalopoda, na.rm = T),
            se_wt_cpue_cephalopoda_me = sd_wt_cpue_cephalopoda_me / sqrt(n),
            
            mean_wt_cpue_fish_me = mean(mean_wt_cpue_fish, na.rm = T),
            sd_wt_cpue_fish_me = sd(mean_wt_cpue_fish, na.rm = T),
            se_wt_cpue_fish_me = sd_wt_cpue_fish_me / sqrt(n),
            
            mean_cephalopoda_fish_proportion_me = mean(100 * cephalopoda_fish_proportion, na.rm = T),
            sd_cephalopoda_fish_proportion_me = sd(100 * cephalopoda_fish_proportion, na.rm = T),
            se_cephalopoda_fish_proportion_me = sd_cephalopoda_fish_proportion_me / sqrt(n),
            
            mean_a_size_weighted_wt_cpue_cephalopoda_me = mean(mean_a_size_weighted_wt_cpue_cephalopoda, na.rm = T),
            sd_a_size_weighted_wt_cpue_cephalopoda_me = sd(mean_a_size_weighted_wt_cpue_cephalopoda, na.rm = T),
            se_a_size_weighted_wt_cpue_cephalopoda_me = sd_a_size_weighted_wt_cpue_cephalopoda_me / sqrt(n),
            
            mean_lifespan_weighted_wt_cpue_cephalopoda_me = mean(mean_lifespan_weighted_wt_cpue_cephalopoda, na.rm = T),
            sd_lifespan_weighted_wt_cpue_cephalopoda_me = sd(mean_lifespan_weighted_wt_cpue_cephalopoda, na.rm = T),
            se_lifespan_weighted_wt_cpue_cephalopoda_me = sd_lifespan_weighted_wt_cpue_cephalopoda_me / sqrt(n),
            
            mean_A_growth_weighted_wt_cpue_cephalopoda_me = mean(mean_A_growth_weighted_wt_cpue_cephalopoda, na.rm = T),
            sd_A_growth_weighted_wt_cpue_cephalopoda_me = sd(mean_A_growth_weighted_wt_cpue_cephalopoda, na.rm = T),
            se_A_growth_weighted_wt_cpue_cephalopoda_me = sd_A_growth_weighted_wt_cpue_cephalopoda_me / sqrt(n),
            
            sd_Temp_C = sd(Temp_C, na.rm = T),
            se_Temp_C = sd_Temp_C / sqrt(n),
            Temp_C = mean(Temp_C),
            sd_depth = sd(depth, na.rm = T),
            se_depth = sd_depth / sqrt(n),
            depth = mean(depth),
            sd_z_prod = sd(z_prod, na.rm = T),
            se_z_prod = sd_z_prod / sqrt(n),
            z_prod = mean(z_prod)) 

stats_df4 <- stats_df %>%
  filter(mean_wt_cpue_cephalopoda > 0)


#-------------------------------------------------------------------------------------------

# Check variance inflation factors:
vif <- lm(mean_wt_cpue_cephalopoda ~ depth + z_prod + Temp_C + mean_wt_cpue_fish,
  data = stats_df)

vif(vif)

# All VIF values are lower than 2

# Check correlation between variables:
cor_table <- stats_df %>%
  dplyr::select(depth, z_prod, Temp_C, mean_wt_cpue_fish)

cor(cor_table)

# Correlation =< 0.5

# Calculate mean ceph biomass per region:
stats_df %>%
  group_by(coast) %>%
  summarise(mean_ceph = mean(mean_wt_cpue_cephalopoda),
            mean_fish = mean(mean_wt_cpue_fish),
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

stats_df %>%
  summarise(total_fish = sum(mean_wt_cpue_fish),
            total_ceph = sum(mean_wt_cpue_cephalopoda),
            proportion = 100 * total_ceph / (total_fish + total_ceph))


#-------------------------------------------------------------------------------------------


############
# Run stats!

#
# GAM models -----------------------------------------------------------
#

###########################
# Cephalopod proportion:
###########################
gam0 <- gam(cephalopoda_fish_proportion ~ 
              s(Temp_C, k = 3) +
              s(depth, k = 3) +
              s(z_prod, k = 3) + 
              # s(log10(mean_wt_cpue_fish), k = 3) +
              s(coast, bs = "re"),
            correlation = corExp(form = ~ lat + lon),
            data = stats_df,
            method = "REML",
            family = betar())

summary(gam0)
gam.check(gam0)

# Check smooths:
p <- draw(gam0) &
  theme_classic()

p
# ggsave("plots/statistics/ceph_proportion/GAMM_ceph_proportion_smooths.png", p, height = 40 , width = 40, units = "mm", scale = 4)

# Save panels for joint plot:
c1 <-  draw(gam0,          
            ci_alpha = 0.2,       
            ci_col = "orange",         
            smooth_col = "orange")[[3]] +
  ylim(c(-1.5, 1.5)) +
  ylab(expression("P"["Ceph"])) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank()) +
  geom_hline(yintercept = 0, color = "gray")
  
c2 <-  draw(gam0,      
            ci_alpha = 0.2,         
            ci_col = "orange",          
            smooth_col = "orange")[[2]] +
  ylim(c(-1.5, 1.5)) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")

c3 <-  draw(gam0,       
            ci_alpha = 0.2,    
            ci_col = "orange",    
            smooth_col = "orange")[[4]] +
  ylim(c(-1.5, 1.5)) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")


# Check residuals:
# Residuals vs fitted vals:
stats_df %>% 
  mutate(pred_dist = fitted(gam0)) %>% 
  mutate(resid = residuals(gam0, type = "response")) %>% 
  ggplot(aes(x = pred_dist, y = resid)) +
  theme_classic() +
  geom_point(size = 3, shape = 1) + 
  geom_hline(yintercept = 0, color = "gray")


# Residuals vs temperature:
stats_df %>% 
  mutate(resid = residuals(gam0, type = "response")) %>% 
  ggplot(aes(x = Temp_C, y = resid, color = coast)) +
  theme_classic() +
  geom_point(size = 3, shape = 1) + 
  geom_hline(yintercept = 0, color = "gray")

# Residuals vs z_prod:
stats_df %>% 
  mutate(resid = residuals(gam0, type = "response")) %>% 
  ggplot(aes(x = z_prod, y = resid, color = coast)) +
  theme_classic() +
  geom_point(size = 3, shape = 1) + 
  geom_hline(yintercept = 0, color = "gray")

# Residuals vs depth:
stats_df %>% 
  mutate(resid = residuals(gam0, type = "response")) %>% 
  ggplot(aes(x = depth, y = resid, color = coast)) +
  theme_classic() +
  geom_point(size = 3, shape = 1) + 
  geom_hline(yintercept = 0, color = "gray")


#######################
# Cephalopod biomass:
#######################
gam1 <- gam(mean_wt_cpue_cephalopoda ~ 
              s(Temp_C, k = 3) +
              s(depth, k = 3) +
              s(z_prod, k = 3) + 
              # s(log10(mean_wt_cpue_fish), k = 3) +
              s(coast, bs = "re"),
            data = stats_df,
            method = "REML",
            family = nb())

summary(gam1)
gam.check(gam1)

# Check smooths:
p <- draw(gam1) &
  theme_classic()

p
# ggsave("plots/statistics/ceph_biomass/GAMM_ceph_biomass_smooths.png", p, height = 40 , width = 40, units = "mm", scale = 4)

# Save panels for joint plot:
a1 <-  draw(gam1,        
            ci_alpha = 0.2,     
            ci_col = "orange",      
            smooth_col = "orange")[[3]] +
  ylim(c(-6.5, 3))  +
  ylab(expression(paste("B"["Ceph"], " (ln(kg ", km^-2,"))"))) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")

a2 <-  draw(gam1,      
            ci_alpha = 0.2, 
            ci_col = "orange",      
            smooth_col = "orange")[[2]] +
  ylim(c(-6.5, 3))  +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")

a3 <-  draw(gam1,     
            ci_alpha = 0.2,   
            ci_col = "orange",   
            smooth_col = "orange")[[4]] +
  ylim(c(-6.5, 3)) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")


# Check residuals:
# Residuals vs fitted vals:
stats_df %>% 
  mutate(pred_dist = fitted(gam1)) %>% 
  mutate(resid = residuals(gam1, type = "response")) %>% 
  ggplot(aes(x = pred_dist, y = resid)) +
  theme_classic() +
  geom_point(size = 3, shape = 1) + 
  geom_hline(yintercept = 0, color = "gray")


# Residuals vs temperature:
stats_df %>% 
  mutate(resid = residuals(gam1, type = "response")) %>% 
  ggplot(aes(x = Temp_C, y = resid, color = coast)) +
  theme_classic() +
  geom_point(size = 3, shape = 1) + 
  geom_hline(yintercept = 0, color = "gray")

# Residuals vs z_prod:
stats_df %>% 
  mutate(resid = residuals(gam1, type = "response")) %>% 
  ggplot(aes(x = z_prod, y = resid, color = coast)) +
  theme_classic() +
  geom_point(size = 3, shape = 1) + 
  geom_hline(yintercept = 0, color = "gray")

# Residuals vs depth:
stats_df %>% 
  mutate(resid = residuals(gam1, type = "response")) %>% 
  ggplot(aes(x = depth, y = resid, color = coast)) +
  theme_classic() +
  geom_point(size = 3, shape = 1) + 
  geom_hline(yintercept = 0, color = "gray")



#######################
# Fish biomass:
#######################
gam2 <- gam(mean_wt_cpue_fish ~ 
              s(Temp_C, k = 3) +
              s(depth, k = 3) +
              s(z_prod, k = 3) + 
              # s(mean_wt_cpue_cephalopoda, k = 3) +
              s(coast, bs = "re"),
            # correlation = corExp(form = ~ lat + lon),
            data = stats_df,
            method = "REML",
            family = nb())

summary(gam2)
gam.check(gam2)

# Check smooths:
p <- draw(gam2) &
  theme_classic()

p
# ggsave("plots/statistics/fish_biomass/GAMM_fish_biomass_smooths.png", p, height = 40 , width = 40, units = "mm", scale = 4)

# Save panels for joint plot:
b1 <-  draw(gam2, 
            ci_alpha = 0.2,  
            ci_col = "orange",         
            smooth_col = "orange")[[3]] +
  ylim(c(-1.5, 1.6)) +
  ylab(expression(paste("B"["Fish"], " (ln(kg ", km^-2,"))"))) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")

b2 <-  draw(gam2,
            ci_alpha = 0.2,
            ci_col = "orange", 
            smooth_col = "orange")[[2]]  +
  ylim(c(-1.5, 1.6)) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")

b3 <-  draw(gam2,
            ci_alpha = 0.2,
            ci_col = "orange",
            smooth_col = "orange")[[4]] +
  ylim(c(-1.5, 1.6)) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")



# Check residuals:
# Residuals vs fitted vals:
stats_df %>% 
  mutate(pred_dist = fitted(gam2)) %>% 
  mutate(resid = residuals(gam2, type = "response")) %>% 
  ggplot(aes(x = pred_dist, y = resid)) +
  theme_classic() +
  geom_point(size = 3, shape = 1) + 
  geom_hline(yintercept = 0, color = "gray")


# Residuals vs temperature:
stats_df %>% 
  mutate(resid = residuals(gam2, type = "response")) %>% 
  ggplot(aes(x = Temp_C, y = resid, color = coast)) +
  theme_classic() +
  geom_point(size = 3, shape = 1) + 
  geom_hline(yintercept = 0, color = "gray")

# Residuals vs z_prod:
stats_df %>% 
  mutate(resid = residuals(gam2, type = "response")) %>% 
  ggplot(aes(x = z_prod, y = resid, color = coast)) +
  theme_classic() +
  geom_point(size = 3, shape = 1) + 
  geom_hline(yintercept = 0, color = "gray")

# Residuals vs depth:
stats_df %>% 
  mutate(resid = residuals(gam2, type = "response")) %>% 
  ggplot(aes(x = depth, y = resid, color = coast)) +
  theme_classic() +
  geom_point(size = 3, shape = 1) + 
  geom_hline(yintercept = 0, color = "gray")


#######################
# Asymptotic size:
#######################
gam3 <- gam(log10(mean_a_size_weighted_wt_cpue_cephalopoda) ~ 
              s(Temp_C, k = 3) +
              s(depth, k = 3) +
              s(z_prod, k = 3) + 
              # s(log10(mean_wt_cpue_fish), k = 3) +
              s(coast, bs = "re"),
            data = stats_df4,
            method = "REML",
            family = gaussian())

summary(gam3)
gam.check(gam3)

# Check smooths:
p <- draw(gam3) &
  theme_classic()

p
# ggsave("plots/statistics/a_size/GAMM_ceph_a_size_smooths.png", p, height = 40 , width = 40, units = "mm", scale = 4)

# Save panels for joint plot:
d1 <-  draw(gam3,
            ci_alpha = 0.2,
            ci_col = "orange",
            smooth_col = "orange")[[3]] +
  ylim(c(-1.2, .5))  +
  ylab(expression(paste("W"[infinity], " log(g)"))) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")

d2 <-  draw(gam3,
            ci_alpha = 0.2,
            ci_col = "orange",
            smooth_col = "orange")[[2]] +
  ylim(c(-1.2, .5))  +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")

d3 <-  draw(gam3,
            ci_alpha = 0.2,
            ci_col = "orange",
            smooth_col = "orange")[[4]] +
  ylim(c(-1.2, .5)) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")


# Check residuals:
# Residuals vs fitted vals:
stats_df4 %>% 
  mutate(pred_dist = fitted(gam3)) %>% 
  mutate(resid = residuals(gam3, type = "response")) %>% 
  ggplot(aes(x = pred_dist, y = resid)) +
  theme_classic() +
  geom_point(size = 3, shape = 1) + 
  geom_hline(yintercept = 0, color = "gray")


# Residuals vs temperature:
stats_df4 %>% 
  mutate(resid = residuals(gam3, type = "response")) %>% 
  ggplot(aes(x = Temp_C, y = resid, color = coast)) +
  theme_classic() +
  geom_point(size = 3, shape = 1) + 
  geom_hline(yintercept = 0, color = "gray")

# Residuals vs z_prod:
stats_df4 %>% 
  mutate(resid = residuals(gam3, type = "response")) %>% 
  ggplot(aes(x = z_prod, y = resid, color = coast)) +
  theme_classic() +
  geom_point(size = 3, shape = 1) + 
  geom_hline(yintercept = 0, color = "gray")

# Residuals vs depth:
stats_df4 %>% 
  mutate(resid = residuals(gam3, type = "response")) %>% 
  ggplot(aes(x = depth, y = resid, color = coast)) +
  theme_classic() +
  geom_point(size = 3, shape = 1) + 
  geom_hline(yintercept = 0, color = "gray")



#######################
# Lifespan:
#######################
gam4 <- gam(mean_lifespan_weighted_wt_cpue_cephalopoda ~ 
              s(Temp_C, k = 3) +
              s(depth, k = 3) +
              s(z_prod, k = 3) + 
              # s(mean_wt_cpue_fish, k = 3) +
              s(coast, bs = "re"),
            data = stats_df4,
            method = "REML",
            family = gaussian())

summary(gam4)
gam.check(gam4)

# Check smooths:
p <- draw(gam4) &
  theme_classic()

p
# ggsave("plots/statistics/lifespan/GAMM_ceph_lifespan_smooths.png", p, height = 40 , width = 40, units = "mm", scale = 4)

# Save panels for joint plot:
e1 <-  draw(gam4,
            ci_alpha = 0.2,
            ci_col = "orange",
            smooth_col = "orange")[[3]] +
  ylim(c(-10, 8.5)) +
  ylab(expression(paste("R"[0]," (months)"))) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")

e2 <-  draw(gam4,
            ci_alpha = 0.2,
            ci_col = "orange",
            smooth_col = "orange")[[2]] +
  ylim(c(-10, 8.5)) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")

e3 <-  draw(gam4,
            ci_alpha = 0.2,
            ci_col = "orange",
            smooth_col = "orange")[[4]] +
  ylim(c(-10, 8.5)) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")


# Check residuals:
# Residuals vs fitted vals:
stats_df4 %>%
  mutate(pred_dist = fitted(gam4)) %>% 
  mutate(resid = residuals(gam4, type = "response")) %>% 
  ggplot(aes(x = pred_dist, y = resid)) +
  theme_classic() +
  geom_point(size = 3, shape = 1) + 
  geom_hline(yintercept = 0, color = "gray")


# Residuals vs temperature:
stats_df4 %>% 
  mutate(resid = residuals(gam4, type = "response")) %>% 
  ggplot(aes(x = Temp_C, y = resid, color = coast)) +
  theme_classic() +
  geom_point(size = 3, shape = 1) + 
  geom_hline(yintercept = 0, color = "gray")

# Residuals vs z_prod:
stats_df4 %>% 
  mutate(resid = residuals(gam4, type = "response")) %>% 
  ggplot(aes(x = z_prod, y = resid, color = coast)) +
  theme_classic() +
  geom_point(size = 3, shape = 1) + 
  geom_hline(yintercept = 0, color = "gray")

# Residuals vs depth:
stats_df4 %>% 
  mutate(resid = residuals(gam4, type = "response")) %>% 
  ggplot(aes(x = depth, y = resid, color = coast)) +
  theme_classic() +
  geom_point(size = 3, shape = 1) + 
  geom_hline(yintercept = 0, color = "gray")


#######################
# A growth:
#######################
gam5 <- gam(mean_A_growth_weighted_wt_cpue_cephalopoda ~ 
              s(Temp_C, k = 3) +
              s(depth, k = 3) +
              s(z_prod, k = 3) + 
              # s(log10(mean_wt_cpue_fish), k = 3) +
              s(coast, bs = "re"),
            data = stats_df4,
            method = "REML",
            family = gaussian())

summary(gam5)
gam.check(gam5)

# Check smooths:
p <- draw(gam5) &
  theme_classic()

p
# ggsave("plots/statistics/a_growth/GAMM_a_growth_gam0_smooths.png", p, height = 40 , width = 40, units = "mm", scale = 4)

# Save panels for joint plot:
f1 <-  draw(gam5,
            ci_alpha = 0.2,
            ci_col = "orange",
            smooth_col = "orange")[[3]] +
  ylim(c(-7, 10.5)) +
  ylab(expression(paste("A (", g^(1/3), yr^-1,")"))) +
  xlab("Temperature (ºC)") +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")

f2 <-  draw(gam5,
            ci_alpha = 0.2,
            ci_col = "orange",
            smooth_col = "orange")[[2]] +
  ylim(c(-7, 10.5)) +
  ylab("") +
  xlab("Depth (m)") +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")

f3 <-  draw(gam5,
            ci_alpha = 0.2,
            ci_col = "orange",
            smooth_col = "orange")[[4]] +
  ylim(c(-7, 10.5)) +
  ylab("") +
  xlab(expression(paste("Z. Product (g ", m^-2, " ", yr^-1, ")"))) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")


# Combine all smooths in a plot:
p <- (a1 | a2 | a3) /
  (b1 | b2 | b3) / 
  (c1 | c2 | c3) /
  (d1 | d2 | d3) /
  (e1 | e2 | e3) /
  (f1 | f2 | f3) +
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 10)) & 
  theme(plot.background = element_blank())

p
ggsave("plots/statistics/all_smooths.png", p, height = 80 , width = 60, units = "mm", scale = 3)

rm(list = c("a1", "a2", "a3", "b1", "b2", "b3", "c1", "c2", "c3", "d1", "d2", "d3", "e1", "e2", "e3", "f1", "f2", "f3"))

# Check residuals:
# Residuals vs fitted vals:  
stats_df4 %>% 
  mutate(pred_dist = fitted(gam5)) %>% 
  mutate(resid = residuals(gam5, type = "response")) %>% 
  ggplot(aes(x = pred_dist, y = resid)) +
  theme_classic() +
  geom_point(size = 3, shape = 1) + 
  geom_hline(yintercept = 0, color = "gray")


# Residuals vs temperature:
stats_df4 %>% 
  mutate(resid = residuals(gam5, type = "response")) %>% 
  ggplot(aes(x = Temp_C, y = resid, color = coast)) +
  theme_classic() +
  geom_point(size = 3, shape = 1) + 
  geom_hline(yintercept = 0, color = "gray")

# Residuals vs z_prod:
stats_df4 %>% 
  mutate(resid = residuals(gam5, type = "response")) %>% 
  ggplot(aes(x = z_prod, y = resid, color = coast)) +
  theme_classic() +
  geom_point(size = 3, shape = 1) + 
  geom_hline(yintercept = 0, color = "gray")

# Residuals vs depth:
stats_df4 %>% 
  mutate(resid = residuals(gam5, type = "response")) %>% 
  ggplot(aes(x = depth, y = resid, color = coast)) +
  theme_classic() +
  geom_point(size = 3, shape = 1) + 
  geom_hline(yintercept = 0, color = "gray")


#-------------------------------------------------------------------------

# Add predictions to the data frame:
stats_df3 <- stats_df %>% 
  mutate(pred_proportion = predict(gam0, type = "response"), # Same output as fitted()
         pred_ceph_biomass = fitted(gam1),
         pred_fish_biomass = fitted(gam2))

save(stats_df3, file = "data/predictions_biomass_q03.Rdata")

stats_df4 <- stats_df4 %>% 
  mutate(pred_a_size = fitted(gam3),
         pred_lifespan = fitted(gam4),
         pred_a_growth = fitted(gam5))

save(stats_df4, file = "data/predictions_traits_q03.Rdata")

#-------------------------------------------------------------------------

#-----------------------------------------------------------
  
# Make cleaner plots grouping observations by Marine Ecosystems:

# Define color-blind colors:
my_colors <- c("#F8766D", "gray50", "#619CFF")

####################
# Cephalopod proportion:

# Make predictions of range of a variable keeping all other variables constant:
# Loop to obtain partial effect predictions for all the variables:
variables <- c("Temp_C", "depth", "z_prod")
partial_plot_proportion <- data.frame()

test <- stats_df3 %>%
  dplyr::select(Temp_C, depth, z_prod, coast) %>%
  as.data.frame()

for(i in 1:length(variables)){
  var1 <- seq(min(test[,which(colnames(test) == variables[i])]),
              max(test[,which(colnames(test) == variables[i])]), length = nrow(test))
  
  vars <- variables[-i]
  
  var2 <- rep(mean(test[,which(colnames(test) == vars[1])]),length(var1)) 
  var3 <- rep(mean(test[,which(colnames(test) == vars[2])]),length(var1))
  
  coast <- rep("West coast", length(var1)) # set the random effect to the same coast
  
  newdat <- data.frame(var1, var2, var3, coast) 
  
  colnames(newdat)[1:3] <- c(variables[i], vars)
  
  my_predictions <- predict.gam(gam0, newdata = newdat, exclude = "s(coast)", type = "response", se.fit = TRUE) # predict with new data excluding random effect
  my_fit <- my_predictions$fit
  my_se <- my_predictions$se.fit
  
  partial <- data.frame(my_predictions, variable = var1) 
  partial <- partial %>% mutate(variable = var1,
                                variable_id = variables[i])
  
  partial_plot_proportion = rbind(partial_plot_proportion, partial)
  
}


# Temperature:
p <- ggplot() +

  geom_smooth(data = stats_df3, aes(x = Temp_C, y = 100 * pred_proportion), color = "black", se = T, alpha = .4, size = .5) +         # Predicted regression
  
  geom_line(data = subset(partial_plot_proportion, variable_id == "Temp_C"), aes(x = variable, y = 100 * fit), color = "orange", size = 1) + 
  geom_ribbon(data = subset(partial_plot_proportion, variable_id == "Temp_C"), aes(y = 100 * fit, x = variable, ymin = 100 * (fit - se.fit), ymax = 100 * (fit + se.fit)),
              fill = "orange", alpha = .1) +
  
  geom_pointrange(data = stats_df2, aes(x = Temp_C, y = mean_cephalopoda_fish_proportion_me,                               # Average in ME +- se
                                        ymin = mean_cephalopoda_fish_proportion_me - se_cephalopoda_fish_proportion_me,
                                        ymax = mean_cephalopoda_fish_proportion_me + se_cephalopoda_fish_proportion_me,
                                        color = coast)) +
  geom_pointrange(data = stats_df2, aes(x = Temp_C, y = mean_cephalopoda_fish_proportion_me, 
                                        xmin = Temp_C - se_Temp_C,
                                        xmax = Temp_C + se_Temp_C,
                                        color = coast)) +
  scale_color_manual(values = my_colors) +
  ylim(c(-.5, NA)) +
  theme_base() +
  ylab("% ceph. biom.") +
  xlab("Temperature (ºC)") + 
  theme(plot.background = element_blank(),
        legend.title = element_blank()) 

p
# ggsave("plots/statistics/ceph_proportion/GAMM_ceph_proportion_predict_Temp_C_me.png", p, height = 40 , width = 40, units = "mm", scale = 4)

c1 <- p +
  theme(legend.position = "none",
        axis.text.x = element_blank()) +
  xlab("") 


# Depth:
p <- ggplot() +
 
  geom_smooth(data = stats_df3, aes(x = depth, y = 100 * pred_proportion), color = "black", se = T, alpha = .4, size = .5) +         # Predicted regression
  
  geom_line(data = subset(partial_plot_proportion, variable_id == "depth"), aes(x = variable, y = 100 * fit), color = "orange", size = 1) + 
  geom_ribbon(data = subset(partial_plot_proportion, variable_id == "depth"), aes(y = 100 * fit, x = variable, ymin = 100 * (fit - se.fit), ymax = 100 * (fit + se.fit)),
              fill = "orange", alpha = .1) +
  
  geom_pointrange(data = stats_df2, aes(x = depth, y = mean_cephalopoda_fish_proportion_me,                               # Average in ME +- se
                                        ymin = mean_cephalopoda_fish_proportion_me - se_cephalopoda_fish_proportion_me,
                                        ymax = mean_cephalopoda_fish_proportion_me + se_cephalopoda_fish_proportion_me,
                                        color = coast)) +
  geom_pointrange(data = stats_df2, aes(x = depth, y = mean_cephalopoda_fish_proportion_me, 
                                        xmin = depth - se_depth,
                                        xmax = depth + se_depth,
                                        color = coast)) +
  scale_color_manual(values = my_colors) +
  ylim(c(-.5, NA)) +
  theme_base() +
  ylab("% cephalopod biomass") +
  xlab("Depth (m)") + 
  theme(plot.background = element_blank(),
        legend.title = element_blank())

p
# ggsave("plots/statistics/ceph_proportion/GAMM_ceph_proportion_predict_depth_me.png", p, height = 40 , width = 40, units = "mm", scale = 4)

c2 <- p +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  xlab("") +
  ylab("")

# Zoo productivity:
p <- ggplot() +

  geom_smooth(data = stats_df3, aes(x = z_prod, y = 100 * pred_proportion), color = "black", se = T, alpha = .4, size = .5) +        # Predicted regression
  
  geom_line(data = subset(partial_plot_proportion, variable_id == "z_prod"), aes(x = variable, y = 100 * fit), color = "orange", size = 1) + 
  geom_ribbon(data = subset(partial_plot_proportion, variable_id == "z_prod"), aes(y = 100 * fit, x = variable, ymin = 100 * (fit - se.fit), ymax = 100 * (fit + se.fit)),
              fill = "orange", alpha = .1) +
  
  geom_pointrange(data = stats_df2, aes(x = z_prod, y = mean_cephalopoda_fish_proportion_me,                               # Average in ME +- se
                                        ymin = mean_cephalopoda_fish_proportion_me - se_cephalopoda_fish_proportion_me,
                                        ymax = mean_cephalopoda_fish_proportion_me + se_cephalopoda_fish_proportion_me,
                                        color = coast)) +
  geom_pointrange(data = stats_df2, aes(x = z_prod, y = mean_cephalopoda_fish_proportion_me, 
                                        xmin = z_prod - se_z_prod,
                                        xmax = z_prod + se_z_prod,
                                        color = coast)) +
  scale_color_manual(values = my_colors) +
  ylim(c(-.5, NA)) +
  theme_base() +
  ylab("% cephalopod biomass") +
  xlab(expression(paste("Zooplankton productivity (g ", m^-2, " ", yr^-1, ")"))) + 
  theme(plot.background = element_blank(),
        legend.title = element_blank()) 

p
# ggsave("plots/statistics/ceph_proportion/GAMM_ceph_proportion_predict_z_prod_me.png", p, height = 40 , width = 40, units = "mm", scale = 4)

c3 <- p +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  xlab("") +
  ylab("")


# Plot predictions using exagons as points:
# As a function of temperature:
p <- ggplot(data = stats_df3) +
  geom_point(aes(x = Temp_C, y = cephalopoda_fish_proportion, color = coast), shape = 21, alpha = .5)  + # Observed values
  geom_line(data = subset(partial_plot_proportion, variable_id == "Temp_C"), aes(x = variable, y = fit), color = "orange", size = 1) + 
  geom_ribbon(data = subset(partial_plot_proportion, variable_id == "Temp_C"), aes(y = fit, x = variable, ymin = (fit - se.fit), ymax = (fit + se.fit)),
              fill = "orange", alpha = .1) +
  # geom_point(aes(x = Temp_C, y = pred_proportion, color = coast)) +                                      # Predicted values
  geom_smooth(aes(x = Temp_C, y = pred_proportion), color = "black", se = F) +                           # Predicted regression
  theme_classic()

p
# ggsave("plots/statistics/ceph_proportion/GAMM_ceph_proportion_predict_Temp_C.png", p, height = 40 , width = 40, units = "mm", scale = 4)

# As a function of depth:
p <- ggplot(data = stats_df3) +
  geom_point(aes(x = depth, y = cephalopoda_fish_proportion, color = coast), shape = 21, alpha = .5)  +  # Observed values
  geom_line(data = subset(partial_plot_proportion, variable_id == "depth"), aes(x = variable, y = fit), color = "orange", size = 1) + 
  geom_ribbon(data = subset(partial_plot_proportion, variable_id == "depth"), aes(y = fit, x = variable, ymin = (fit - se.fit), ymax = (fit + se.fit)),
              fill = "orange", alpha = .1) +
  # geom_point(aes(x = depth, y = pred_proportion, color = coast)) +                                       # Predicted values
  geom_smooth(aes(x = depth, y = pred_proportion), color = "black", se = F) +                            # Predicted regression
  theme_classic()

p
# ggsave("plots/statistics/ceph_proportion/GAMM_ceph_proportion_predict_depth.png", p, height = 40 , width = 40, units = "mm", scale = 4)

# As a function of z_prod:
p <- ggplot(data = stats_df3) +
  geom_point(aes(x = z_prod, y = cephalopoda_fish_proportion, color = coast), shape = 21, alpha = .5)  +  # Observed values
  geom_line(data = subset(partial_plot_proportion, variable_id == "z_prod"), aes(x = variable, y = fit), color = "orange", size = 1) + 
  geom_ribbon(data = subset(partial_plot_proportion, variable_id == "z_prod"), aes(y = fit, x = variable, ymin = (fit - se.fit), ymax = (fit + se.fit)),
              fill = "orange", alpha = .1) +
  # geom_point(aes(x = z_prod, y = pred_proportion, color = coast)) +                                       # Predicted values
  geom_smooth(aes(x = z_prod, y = pred_proportion), color = "black", se = F) +                            # Predicted regression
  theme_classic()

p
# ggsave("plots/statistics/ceph_proportion/GAMM_ceph_proportion_predict_z_prod.png", p, height = 40 , width = 40, units = "mm", scale = 4)



#####################
# Cephalopod biomass:

# Make predictions of range of a variable keeping all other variables constant:
# Loop to obtain partial effect predictions for all the variables:
variables <- c("Temp_C", "depth", "z_prod")
partial_plot_ceph_bio <- data.frame()

test <- stats_df3 %>%
  dplyr::select(Temp_C, depth, z_prod, coast) %>%
  as.data.frame()

for(i in 1:length(variables)){
  var1 <- seq(min(test[,which(colnames(test) == variables[i])]),
              max(test[,which(colnames(test) == variables[i])]), length = nrow(test))
  
  vars <- variables[-i]
  
  var2 <- rep(mean(test[,which(colnames(test) == vars[1])]),length(var1)) 
  var3 <- rep(mean(test[,which(colnames(test) == vars[2])]),length(var1))
  
  coast <- rep("West coast", length(var1)) # set the random effect to the same coast
  
  newdat <- data.frame(var1, var2, var3, coast) 
  
  colnames(newdat)[1:3] <- c(variables[i], vars)
  
  my_predictions <- predict.gam(gam1, newdata = newdat, exclude = "s(coast)", type = "response", se.fit = TRUE) # predict with new data excluding random effect
  my_fit <- my_predictions$fit
  my_se <- my_predictions$se.fit
  
  partial <- data.frame(my_predictions, variable = var1) 
  partial <- partial %>% mutate(variable = var1,
                                variable_id = variables[i])
  
  partial_plot_ceph_bio = rbind(partial_plot_ceph_bio, partial)
  
}


# Temperature:
p <- ggplot() +
  
  geom_smooth(data = stats_df3, aes(x = Temp_C, y = pred_ceph_biomass), color = "black", se = T, alpha = .4, size = .5) +         # Predicted regression
  
  geom_line(data = subset(partial_plot_ceph_bio, variable_id == "Temp_C"), aes(x = variable, y = fit), color = "orange", size = 1) + 
  geom_ribbon(data = subset(partial_plot_ceph_bio, variable_id == "Temp_C"), aes(y = fit, x = variable, ymin = (fit - se.fit), ymax =  (fit + se.fit)),
              fill = "orange", alpha = .1) +
  
  
  geom_pointrange(data = stats_df2, aes(x = Temp_C, y = mean_wt_cpue_cephalopoda_me,                               # Average in ME +- se
                                        ymin = mean_wt_cpue_cephalopoda_me - se_wt_cpue_cephalopoda_me,
                                        ymax = mean_wt_cpue_cephalopoda_me + se_wt_cpue_cephalopoda_me,
                                        color = coast)) +
  geom_pointrange(data = stats_df2, aes(x = Temp_C, y = mean_wt_cpue_cephalopoda_me, 
                                        xmin = Temp_C - se_Temp_C,
                                        xmax = Temp_C + se_Temp_C,
                                        color = coast)) +
  scale_color_manual(values = my_colors) +
  theme_base() +
  ylab(expression(paste("Ceph. biom. (kg ", km^-2, ")"))) +
  xlab("Temperature (ºC)") + 
  theme(plot.background = element_blank(),
        legend.title = element_blank()) +
  scale_y_log10(limits = c(.2, 3700)) 
  # scale_y_continuous(trans = log_trans(), 
  #                    breaks = trans_breaks("log", function(x) exp(x)),
  #                    labels = label_number(accuracy = 1))


p
# ggsave("plots/statistics/ceph_biomass/GAMM_ceph_biomass_predict_Temp_C_me.png", p, height = 40 , width = 40, units = "mm", scale = 4)

a1 <- p +
  theme(legend.position = "none",
        axis.text.x = element_blank()) +
  xlab("") 


# Depth:
p <- ggplot() +

  geom_smooth(data = stats_df3, aes(x = depth, y = pred_ceph_biomass), color = "black", se = T, alpha = .4, size = .5) +         # Predicted regression
  
  geom_line(data = subset(partial_plot_ceph_bio, variable_id == "depth"), aes(x = variable, y = fit), color = "orange", size = 1) + 
  geom_ribbon(data = subset(partial_plot_ceph_bio, variable_id == "depth"), aes(y = fit, x = variable, ymin = (fit - se.fit), ymax =  (fit + se.fit)),
              fill = "orange", alpha = .1) +
  
  
  geom_pointrange(data = stats_df2, aes(x = depth, y = mean_wt_cpue_cephalopoda_me,                               # Average in ME +- se
                                        ymin = mean_wt_cpue_cephalopoda_me - se_wt_cpue_cephalopoda_me,
                                        ymax = mean_wt_cpue_cephalopoda_me + se_wt_cpue_cephalopoda_me,
                                        color = coast)) +
  geom_pointrange(data = stats_df2, aes(x = depth, y = mean_wt_cpue_cephalopoda_me, 
                                        xmin = depth - se_depth,
                                        xmax = depth + se_depth,
                                        color = coast)) +
  scale_color_manual(values = my_colors) +
  theme_base() +
  ylab(expression(paste("Ceph. biomass (kg ", km^-2, ")"))) +
  xlab("Depth (m)") + 
  theme(plot.background = element_blank(),
        legend.title = element_blank()) +
  scale_y_log10(limits = c(.2, 3700))
  # scale_y_continuous(trans = log_trans(), 
  #                    breaks = trans_breaks("log", function(x) exp(x)),
  #                    labels = label_number(accuracy = 1))

p
# ggsave("plots/statistics/ceph_biomass/GAMM_ceph_biomass_predict_depth_me.png", p, height = 40 , width = 40, units = "mm", scale = 4)

a2 <- p +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  xlab("") +
  ylab("")

# Zoo productivity:
p <- ggplot() +
 
  geom_smooth(data = stats_df3, aes(x = z_prod, y = pred_ceph_biomass), color = "black", se = T, alpha = .4, size = .5) +        # Predicted regression
  
  geom_line(data = subset(partial_plot_ceph_bio, variable_id == "z_prod"), aes(x = variable, y = fit), color = "orange", size = 1) + 
  geom_ribbon(data = subset(partial_plot_ceph_bio, variable_id == "z_prod"), aes(y = fit, x = variable, ymin = (fit - se.fit), ymax =  (fit + se.fit)),
              fill = "orange", alpha = .1) +
  
  geom_pointrange(data = stats_df2, aes(x = z_prod, y = mean_wt_cpue_cephalopoda_me,                               # Average in ME +- se
                                        ymin = mean_wt_cpue_cephalopoda_me - se_wt_cpue_cephalopoda_me,
                                        ymax = mean_wt_cpue_cephalopoda_me + se_wt_cpue_cephalopoda_me,
                                        color = coast)) +
  geom_pointrange(data = stats_df2, aes(x = z_prod, y = mean_wt_cpue_cephalopoda_me, 
                                        xmin = z_prod - se_z_prod,
                                        xmax = z_prod + se_z_prod,
                                        color = coast)) +
  scale_color_manual(values = my_colors) +
  theme_base() +
  ylab(expression(paste("Cephalopod biomass (kg ", km^-2, ")"))) +
  xlab(expression(paste("Zooplankton productivity (g ", m^-2, " ", yr^-1, ")"))) + 
  theme(plot.background = element_blank(),
        legend.title = element_blank()) +
  scale_y_log10(limits = c(.2, 3700))
# scale_y_continuous(trans = log_trans(), 
#                    breaks = trans_breaks("log", function(x) exp(x)),
#                    labels = label_number(accuracy = 1))

p
# ggsave("plots/statistics/ceph_biomass/GAMM_ceph_biomass_predict_z_prod_me.png", p, height = 40 , width = 40, units = "mm", scale = 4)

a3 <- p +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  xlab("") +
  ylab("")


# Plot all dot hexagons:
# As a function of temperature:
p <- ggplot(data = stats_df3) +
  geom_point(aes(x = Temp_C, y = mean_wt_cpue_cephalopoda, color = coast), shape = 21, alpha = .5)  + # Observed values
  geom_line(data = subset(partial_plot_ceph_bio, variable_id == "Temp_C"), aes(x = variable, y = fit), color = "orange", size = 1) + 
  geom_ribbon(data = subset(partial_plot_ceph_bio, variable_id == "Temp_C"), aes(y = fit, x = variable, ymin = (fit - se.fit), ymax = (fit + se.fit)),
              fill = "orange", alpha = .1) +
  # geom_point(aes(x = Temp_C, y = pred_proportion, color = coast)) +                                      # Predicted values
  geom_smooth(aes(x = Temp_C, y = pred_ceph_biomass), color = "black", se = F) +                           # Predicted regression
  theme_classic() +
  scale_y_log10()

p
# ggsave("plots/statistics/ceph_biomass/GAMM_ceph_biomass_predict_Temp_C.png", p, height = 40 , width = 40, units = "mm", scale = 4)

# As a function of depth:
p <- ggplot(data = stats_df3) +
  geom_point(aes(x = depth, y = mean_wt_cpue_cephalopoda, color = coast), shape = 21, alpha = .5)  +  # Observed values
  geom_line(data = subset(partial_plot_ceph_bio, variable_id == "depth"), aes(x = variable, y = fit), color = "orange", size = 1) + 
  geom_ribbon(data = subset(partial_plot_ceph_bio, variable_id == "depth"), aes(y = fit, x = variable, ymin = (fit - se.fit), ymax = (fit + se.fit)),
              fill = "orange", alpha = .1) +
  # geom_point(aes(x = depth, y = pred_proportion, color = coast)) +                                       # Predicted values
  geom_smooth(aes(x = depth, y = pred_ceph_biomass), color = "black", se = F) +                            # Predicted regression
  theme_classic() +
  scale_y_log10()

p
# ggsave("plots/statistics/ceph_biomass/GAMM_ceph_biomass_predict_depth.png", p, height = 40 , width = 40, units = "mm", scale = 4)

# As a function of z_prod:
p <- ggplot(data = stats_df3) +
  geom_point(aes(x = z_prod, y = mean_wt_cpue_cephalopoda, color = coast), shape = 21, alpha = .5)  +  # Observed values
  geom_line(data = subset(partial_plot_ceph_bio, variable_id == "z_prod"), aes(x = variable, y = fit), color = "orange", size = 1) + 
  geom_ribbon(data = subset(partial_plot_ceph_bio, variable_id == "z_prod"), aes(y = fit, x = variable, ymin = (fit - se.fit), ymax = (fit + se.fit)),
              fill = "orange", alpha = .1) +
  # geom_point(aes(x = z_prod, y = pred_proportion, color = coast)) +                                       # Predicted values
  geom_smooth(aes(x = z_prod, y = pred_ceph_biomass), color = "black", se = F) +                            # Predicted regression
  theme_classic() +
  scale_y_log10()

p
# ggsave("plots/statistics/ceph_biomass/GAMM_ceph_biomass_predict_z_prod.png", p, height = 40 , width = 40, units = "mm", scale = 4)


#####################
# Fish biomass:

# Make predictions of range of a variable keeping all other variables constant:
# Loop to obtain partial effect predictions for all the variables:
variables <- c("Temp_C", "depth", "z_prod")
partial_plot_fish_bio <- data.frame()

test <- stats_df3 %>%
  dplyr::select(Temp_C, depth, z_prod, coast) %>%
  as.data.frame()

for(i in 1:length(variables)){
  var1 <- seq(min(test[,which(colnames(test) == variables[i])]),
              max(test[,which(colnames(test) == variables[i])]), length = nrow(test))
  
  vars <- variables[-i]
  
  var2 <- rep(mean(test[,which(colnames(test) == vars[1])]),length(var1)) 
  var3 <- rep(mean(test[,which(colnames(test) == vars[2])]),length(var1))
  
  coast <- rep("West coast", length(var1)) # set the random effect to the same coast
  
  newdat <- data.frame(var1, var2, var3, coast) 
  
  colnames(newdat)[1:3] <- c(variables[i], vars)
  
  my_predictions <- predict.gam(gam2, newdata = newdat, exclude = "s(coast)", type = "response", se.fit = TRUE) # predict with new data excluding random effect
  my_fit <- my_predictions$fit
  my_se <- my_predictions$se.fit
  
  partial <- data.frame(my_predictions, variable = var1) 
  partial <- partial %>% mutate(variable = var1,
                                variable_id = variables[i])
  
  partial_plot_fish_bio = rbind(partial_plot_fish_bio, partial)
  
}

# Temperature:
p <- ggplot() +
  
  geom_smooth(data = stats_df3, aes(x = Temp_C, y = pred_fish_biomass), color = "black", se = T, alpha = .4, size = .5) +         # Predicted regression
  
  geom_line(data = subset(partial_plot_fish_bio, variable_id == "Temp_C"), aes(x = variable, y = fit), color = "orange", size = 1) + 
  geom_ribbon(data = subset(partial_plot_fish_bio, variable_id == "Temp_C"), aes(y = fit, x = variable, ymin = (fit - se.fit), ymax =  (fit + se.fit)),
              fill = "orange", alpha = .1) +
  
  geom_pointrange(data = stats_df2, aes(x = Temp_C, y = mean_wt_cpue_fish_me,                               # Average in ME +- se
                                        ymin = mean_wt_cpue_fish_me - se_wt_cpue_cephalopoda_me,
                                        ymax = mean_wt_cpue_fish_me + se_wt_cpue_cephalopoda_me,
                                        color = coast)) +
  geom_pointrange(data = stats_df2, aes(x = Temp_C, y = mean_wt_cpue_fish_me, 
                                        xmin = Temp_C - se_Temp_C,
                                        xmax = Temp_C + se_Temp_C,
                                        color = coast)) +
  scale_color_manual(values = my_colors) +
  theme_base() +
  ylab(expression(paste("Fish biom. (kg ", km^-2, ")"))) +
  xlab("Temperature (ºC)") + 
  theme(plot.background = element_blank(),
        legend.title = element_blank()) +
  scale_y_log10(limits = c(NA, 2*10^5))
# scale_y_continuous(trans = log_trans(), 
#                    breaks = trans_breaks("log", function(x) exp(x)),
#                    labels = label_number(accuracy = 1))
# +
#   ylim(0, NA)

p
# ggsave("plots/statistics/fish_biomass/GAMM_fish_biomass_predict_Temp_C_me.png", p, height = 40 , width = 40, units = "mm", scale = 4)

b1 <- p +
  theme(legend.position = "none",
        axis.text.x = element_blank()) +
  xlab("") 

# Depth:
p <- ggplot() +
  
  geom_smooth(data = stats_df3, aes(x = depth, y = pred_fish_biomass), color = "black", se = T, alpha = .4, size = .5) +         # Predicted regression
  
  geom_line(data = subset(partial_plot_fish_bio, variable_id == "depth"), aes(x = variable, y = fit), color = "orange", size = 1) + 
  geom_ribbon(data = subset(partial_plot_fish_bio, variable_id == "depth"), aes(y = fit, x = variable, ymin = (fit - se.fit), ymax =  (fit + se.fit)),
              fill = "orange", alpha = .1) +
  
  geom_pointrange(data = stats_df2, aes(x = depth, y = mean_wt_cpue_fish_me,                               # Average in ME +- se
                                        ymin = mean_wt_cpue_fish_me - se_wt_cpue_cephalopoda_me,
                                        ymax = mean_wt_cpue_fish_me + se_wt_cpue_cephalopoda_me,
                                        color = coast)) +
  geom_pointrange(data = stats_df2, aes(x = depth, y = mean_wt_cpue_fish_me, 
                                        xmin = depth - se_depth,
                                        xmax = depth + se_depth,
                                        color = coast)) +
  scale_color_manual(values = my_colors) +
  theme_base() +
  ylab(expression(paste("Fish biomass (kg ", km^-2, ")"))) +
  xlab("Depth (m)") + 
  theme(plot.background = element_blank(),
        legend.title = element_blank()) +
  scale_y_log10(limits = c(NA, 2*10^5))
# scale_y_continuous(trans = log_trans(), 
#                    breaks = trans_breaks("log", function(x) exp(x)),
#                    labels = label_number(accuracy = 1))

p
# ggsave("plots/statistics/fish_biomass/GAMM_fish_biomass_predict_depth_me.png", p, height = 40 , width = 40, units = "mm", scale = 4)

b2 <- p +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  xlab("") +
  ylab("")

# Zoo productivity:
p <- ggplot() +
 
  geom_smooth(data = stats_df3, aes(x = z_prod, y = pred_fish_biomass), color = "black", se = T, alpha = .4, size = .5) +        # Predicted regression
  
  geom_line(data = subset(partial_plot_fish_bio, variable_id == "z_prod"), aes(x = variable, y = fit), color = "orange", size = 1) + 
  geom_ribbon(data = subset(partial_plot_fish_bio, variable_id == "z_prod"), aes(y = fit, x = variable, ymin = (fit - se.fit), ymax =  (fit + se.fit)),
              fill = "orange", alpha = .1) +
  
  geom_pointrange(data = stats_df2, aes(x = z_prod, y = mean_wt_cpue_fish_me,                               # Average in ME +- se
                                        ymin = mean_wt_cpue_fish_me - se_wt_cpue_cephalopoda_me,
                                        ymax = mean_wt_cpue_fish_me + se_wt_cpue_cephalopoda_me,
                                        color = coast)) +
  geom_pointrange(data = stats_df2, aes(x = z_prod, y = mean_wt_cpue_fish_me, 
                                        xmin = z_prod - se_z_prod,
                                        xmax = z_prod + se_z_prod,
                                        color = coast)) +
  scale_color_manual(values = my_colors) +
  theme_base() +
  ylab(expression(paste("Fish biomass (kg ", km^-2, ")"))) +
  xlab(expression(paste("Zooplankton productivity (g ", m^-2, " ", yr^-1, ")"))) + 
  theme(plot.background = element_blank(),
        legend.title = element_blank()) +
  scale_y_log10(limits = c(NA, 2*10^5))
# scale_y_continuous(trans = log_trans(), 
#                    breaks = trans_breaks("log", function(x) exp(x)),
#                    labels = label_number(accuracy = 1))

p
# ggsave("plots/statistics/fish_biomass/GAMM_fish_biomass_predict_z_prod_me.png", p, height = 40 , width = 40, units = "mm", scale = 4)

b3 <- p +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  xlab("") +
  ylab("")


#####################
# Cephalopod asymptotic size:

# Make predictions of range of a variable keeping all other variables constant:
# Loop to obtain partial effect predictions for all the variables:
variables <- c("Temp_C", "depth", "z_prod")
partial_plot_a_size <- data.frame()

test <- stats_df4 %>%
  dplyr::select(Temp_C, depth, z_prod, coast) %>%
  as.data.frame()

for(i in 1:length(variables)){
  var1 <- seq(min(test[,which(colnames(test) == variables[i])]),
              max(test[,which(colnames(test) == variables[i])]), length = nrow(test))
  
  vars <- variables[-i]
  
  var2 <- rep(mean(test[,which(colnames(test) == vars[1])]),length(var1)) 
  var3 <- rep(mean(test[,which(colnames(test) == vars[2])]),length(var1))
  
  coast <- rep("West coast", length(var1)) # set the random effect to the same coast
  
  newdat <- data.frame(var1, var2, var3, coast) 
  
  colnames(newdat)[1:3] <- c(variables[i], vars)
  
  my_predictions <- predict.gam(gam3, newdata = newdat, exclude = "s(coast)", type = "response", se.fit = TRUE) # predict with new data excluding random effect
  my_fit <- my_predictions$fit
  my_se <- my_predictions$se.fit
  
  partial <- data.frame(my_predictions, variable = var1) 
  partial <- partial %>% mutate(variable = var1,
                                variable_id = variables[i])
  
  partial_plot_a_size = rbind(partial_plot_a_size, partial)
  
}

# Temperature:
p <- ggplot() +

  geom_smooth(data = stats_df4, aes(x = Temp_C, y = 10^pred_a_size), color = "black", se = T, alpha = .4, size = .5) +         # Predicted regression
  
  geom_line(data = subset(partial_plot_a_size, variable_id == "Temp_C"), aes(x = variable, y = 10^fit), color = "orange", size = 1) + 
  geom_ribbon(data = subset(partial_plot_a_size, variable_id == "Temp_C"), aes(y = 10^fit, x = variable, ymin = 10^(fit - se.fit), ymax =  10^(fit + se.fit)),
              fill = "orange", alpha = .1) +
  
  geom_pointrange(data = stats_df2, aes(x = Temp_C, y = mean_a_size_weighted_wt_cpue_cephalopoda_me,                               # Average in ME +- se
                                        ymin = mean_a_size_weighted_wt_cpue_cephalopoda_me - se_a_size_weighted_wt_cpue_cephalopoda_me,
                                        ymax = mean_a_size_weighted_wt_cpue_cephalopoda_me + se_a_size_weighted_wt_cpue_cephalopoda_me,
                                        color = coast)) +
  geom_pointrange(data = stats_df2, aes(x = Temp_C, y = mean_a_size_weighted_wt_cpue_cephalopoda_me, 
                                        xmin = Temp_C - se_Temp_C,
                                        xmax = Temp_C + se_Temp_C,
                                        color = coast)) +
  scale_color_manual(values = my_colors) +
  theme_base() +
  ylab("Asymptotic size (g)") +
  xlab("Temperature (ºC)") + 
  theme(plot.background = element_blank(),
        legend.title = element_blank()) +
  scale_y_log10(limits = c(100, NA), breaks = c(1000, 10000)) 

p
# ggsave("plots/statistics/a_size/GAMM_a_size_predict_Temp_C_me.png", p, height = 40 , width = 40, units = "mm", scale = 4)

d1 <- p +
  theme(legend.position = "none",
        axis.text.x = element_blank()) +
  xlab("") 

# Depth:
p <- ggplot() +
  
  geom_smooth(data = stats_df4, aes(x = depth, y = 10^pred_a_size), color = "black", se = T, alpha = .4, size = .5) +         # Predicted regression
  
  geom_line(data = subset(partial_plot_a_size, variable_id == "depth"), aes(x = variable, y = 10^fit), color = "orange", size = 1) + 
  geom_ribbon(data = subset(partial_plot_a_size, variable_id == "depth"), aes(y = 10^fit, x = variable, ymin = 10^(fit - se.fit), ymax =  10^(fit + se.fit)),
              fill = "orange", alpha = .1) +
  
  geom_pointrange(data = stats_df2, aes(x = depth, y = mean_a_size_weighted_wt_cpue_cephalopoda_me,                               # Average in ME +- se
                                        ymin = mean_a_size_weighted_wt_cpue_cephalopoda_me - se_a_size_weighted_wt_cpue_cephalopoda_me,
                                        ymax = mean_a_size_weighted_wt_cpue_cephalopoda_me + se_a_size_weighted_wt_cpue_cephalopoda_me,
                                        color = coast)) +
  geom_pointrange(data = stats_df2, aes(x = depth, y = mean_a_size_weighted_wt_cpue_cephalopoda_me, 
                                        xmin = depth - se_depth,
                                        xmax = depth + se_depth,
                                        color = coast)) +
  scale_color_manual(values = my_colors) +
  theme_base() +
  ylab("Asymptotic size (g)") +
  xlab("Depth (m)") + 
  theme(plot.background = element_blank(),
        legend.title = element_blank()) +
  scale_y_log10(limits = c(100, NA), breaks = c(1000, 10000)) 

p
# ggsave("plots/statistics/a_size/GAMM_a_size_predict_depth_me.png", p, height = 40 , width = 40, units = "mm", scale = 4)

d2 <- p +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  xlab("") +
  ylab("")

# Zoo productivity:
p <- ggplot() +

  geom_smooth(data = stats_df4, aes(x = z_prod, y = 10^pred_a_size), color = "black", se = T, alpha = .4, size = .5) +        # Predicted regression
  
  geom_line(data = subset(partial_plot_a_size, variable_id == "z_prod"), aes(x = variable, y = 10^fit), color = "orange", size = 1) + 
  geom_ribbon(data = subset(partial_plot_a_size, variable_id == "z_prod"), aes(y = 10^fit, x = variable, ymin = 10^(fit - se.fit), ymax =  10^(fit + se.fit)),
              fill = "orange", alpha = .1) +
  
  geom_pointrange(data = stats_df2, aes(x = z_prod, y = mean_a_size_weighted_wt_cpue_cephalopoda_me,                               # Average in ME +- se
                                        ymin = mean_a_size_weighted_wt_cpue_cephalopoda_me - se_a_size_weighted_wt_cpue_cephalopoda_me,
                                        ymax = mean_a_size_weighted_wt_cpue_cephalopoda_me + se_a_size_weighted_wt_cpue_cephalopoda_me,
                                        color = coast)) +
  geom_pointrange(data = stats_df2, aes(x = z_prod, y = mean_a_size_weighted_wt_cpue_cephalopoda_me, 
                                        xmin = z_prod - se_z_prod,
                                        xmax = z_prod + se_z_prod,
                                        color = coast)) +
  scale_color_manual(values = my_colors) +
  theme_base() +
  ylab("Asymptotic size (g)") +
  xlab(expression(paste("Zooplankton productivity (g ", m^-2, " ", yr^-1, ")"))) + 
  theme(plot.background = element_blank(),
        legend.title = element_blank()) +
  scale_y_log10(limits = c(100, NA), breaks = c(1000, 10000)) 

p
# ggsave("plots/statistics/a_size/GAMM_a_size_predict_z_prod_me.png", p, height = 40 , width = 40, units = "mm", scale = 4)

d3 <- p +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  xlab("") +
  ylab("")



# Plot all dot hexagons:
# As a function of temperature:
p <- ggplot(data = stats_df4) +
  geom_point(aes(x = Temp_C, y = mean_a_size_weighted_wt_cpue_cephalopoda, color = coast), shape = 21, alpha = .5)  + # Observed values
  geom_line(data = subset(partial_plot_a_size, variable_id == "Temp_C"), aes(x = variable, y = 10^fit), color = "orange", size = 1) + 
  geom_ribbon(data = subset(partial_plot_a_size, variable_id == "Temp_C"), aes(y =10^ fit, x = variable, ymin = (10^fit - 10^se.fit), ymax = (10^fit + 10^se.fit)),
              fill = "orange", alpha = .1) +
  # geom_point(aes(x = Temp_C, y = pred_proportion, color = coast)) +                                      # Predicted values
  geom_smooth(aes(x = Temp_C, y = 10^pred_a_size), color = "black", se = F) +                           # Predicted regression
  theme_classic()

p
# ggsave("plots/statistics/a_size/GAMM_a_size_predict_Temp_C.png", p, height = 40 , width = 40, units = "mm", scale = 4)

# As a function of depth:
p <- ggplot(data = stats_df4) +
  geom_point(aes(x = depth, y = mean_a_size_weighted_wt_cpue_cephalopoda, color = coast), shape = 21, alpha = .5)  +  # Observed values
  geom_line(data = subset(partial_plot_a_size, variable_id == "depth"), aes(x = variable, y = 10^fit), color = "orange", size = 1) + 
  geom_ribbon(data = subset(partial_plot_a_size, variable_id == "depth"), aes(y = 10^fit, x = variable, ymin = (10^fit - 10^se.fit), ymax = (10^fit + 10^se.fit)),
              fill = "orange", alpha = .1) +
  # geom_point(aes(x = depth, y = pred_proportion, color = coast)) +                                       # Predicted values
  geom_smooth(aes(x = depth, y = 10^pred_a_size), color = "black", se = F) +                            # Predicted regression
  theme_classic() 
p
# ggsave("plots/statistics/a_size/GAMM_a_size_predict_depth.png", p, height = 40 , width = 40, units = "mm", scale = 4)

# As a function of z_prod:
p <- ggplot(data = stats_df4) +
  geom_point(aes(x = z_prod, y = mean_a_size_weighted_wt_cpue_cephalopoda, color = coast), shape = 21, alpha = .5)  +  # Observed values
  geom_line(data = subset(partial_plot_a_size, variable_id == "z_prod"), aes(x = variable, y = 10^fit), color = "orange", size = 1) + 
  geom_ribbon(data = subset(partial_plot_a_size, variable_id == "z_prod"), aes(y = 10^fit, x = variable, ymin = (10^fit - 10^se.fit), ymax = (10^fit + 10^se.fit)),
              fill = "orange", alpha = .1) +
  # geom_point(aes(x = z_prod, y = pred_proportion, color = coast)) +                                       # Predicted values
  geom_smooth(aes(x = z_prod, y = 10^pred_a_size), color = "black", se = F) +                            # Predicted regression
  theme_classic() 

p
# ggsave("plots/statistics/a_size/GAMM_a_size_predict_z_prod.png", p, height = 40 , width = 40, units = "mm", scale = 4)


#####################
# Cephalopod lifespan:


# Make predictions of range of a variable keeping all other variables constant:
# Loop to obtain partial effect predictions for all the variables:
variables <- c("Temp_C", "depth", "z_prod")
partial_plot_lifespan <- data.frame()

test <- stats_df4 %>%
  dplyr::select(Temp_C, depth, z_prod, coast) %>%
  as.data.frame()

for(i in 1:length(variables)){
  var1 <- seq(min(test[,which(colnames(test) == variables[i])]),
              max(test[,which(colnames(test) == variables[i])]), length = nrow(test))
  
  vars <- variables[-i]
  
  var2 <- rep(mean(test[,which(colnames(test) == vars[1])]),length(var1)) 
  var3 <- rep(mean(test[,which(colnames(test) == vars[2])]),length(var1))
  
  coast <- rep("West coast", length(var1)) # set the random effect to the same coast
  
  newdat <- data.frame(var1, var2, var3, coast) 
  
  colnames(newdat)[1:3] <- c(variables[i], vars)
  
  my_predictions <- predict.gam(gam4, newdata = newdat, exclude = "s(coast)", type = "response", se.fit = TRUE) # predict with new data excluding random effect
  my_fit <- my_predictions$fit
  my_se <- my_predictions$se.fit
  
  partial <- data.frame(my_predictions, variable = var1) 
  partial <- partial %>% mutate(variable = var1,
                                variable_id = variables[i])
  
  partial_plot_lifespan = rbind(partial_plot_lifespan, partial)
  
}


# Temperature:
p <- ggplot() +

  geom_smooth(data = stats_df4, aes(x = Temp_C, y = pred_lifespan), color = "black", se = T, alpha = .4, size = .5) +         # Predicted regression
  
  geom_line(data = subset(partial_plot_lifespan, variable_id == "Temp_C"), aes(x = variable, y = fit), color = "orange", size = 1) + 
  geom_ribbon(data = subset(partial_plot_lifespan, variable_id == "Temp_C"), aes(y = fit, x = variable, ymin = (fit - se.fit), ymax =  (fit + se.fit)),
              fill = "orange", alpha = .1) +
  
  geom_pointrange(data = stats_df2, aes(x = Temp_C, y = mean_lifespan_weighted_wt_cpue_cephalopoda_me,                               # Average in ME +- se
                                        ymin = mean_lifespan_weighted_wt_cpue_cephalopoda_me - se_lifespan_weighted_wt_cpue_cephalopoda_me,
                                        ymax = mean_lifespan_weighted_wt_cpue_cephalopoda_me + se_lifespan_weighted_wt_cpue_cephalopoda_me,
                                        color = coast)) +
  geom_pointrange(data = stats_df2, aes(x = Temp_C, y = mean_lifespan_weighted_wt_cpue_cephalopoda_me, 
                                        xmin = Temp_C - se_Temp_C,
                                        xmax = Temp_C + se_Temp_C,
                                        color = coast)) +
  scale_color_manual(values = my_colors) +
  ylim(c(4, 27)) +
  theme_base() +
  ylab("Lifespan (months)") +
  xlab("Temperature (ºC)") + 
  theme(plot.background = element_blank(),
        legend.title = element_blank())
# +
#   ylim(0, NA)

p
# ggsave("plots/statistics/lifespan/GAMM_lifespan_predict_Temp_C_me.png", p, height = 40 , width = 40, units = "mm", scale = 4)

e1 <- p +
  theme(legend.position = "none",
        axis.text.x = element_blank()) +
  xlab("")

# Depth:
p <- ggplot() +
  
  geom_smooth(data = stats_df4, aes(x = depth, y = pred_lifespan), color = "black", se = T, alpha = .4, size = .5) +         # Predicted regression
  
  geom_line(data = subset(partial_plot_lifespan, variable_id == "depth"), aes(x = variable, y = fit), color = "orange", size = 1) + 
  geom_ribbon(data = subset(partial_plot_lifespan, variable_id == "depth"), aes(y = fit, x = variable, ymin = (fit - se.fit), ymax =  (fit + se.fit)),
              fill = "orange", alpha = .1) +
  
  geom_pointrange(data = stats_df2, aes(x = depth, y = mean_lifespan_weighted_wt_cpue_cephalopoda_me,                               # Average in ME +- se
                                        ymin = mean_lifespan_weighted_wt_cpue_cephalopoda_me - se_lifespan_weighted_wt_cpue_cephalopoda_me,
                                        ymax = mean_lifespan_weighted_wt_cpue_cephalopoda_me + se_lifespan_weighted_wt_cpue_cephalopoda_me,
                                        color = coast)) +
  geom_pointrange(data = stats_df2, aes(x = depth, y = mean_lifespan_weighted_wt_cpue_cephalopoda_me, 
                                        xmin = depth - se_depth,
                                        xmax = depth + se_depth,
                                        color = coast)) +
  scale_color_manual(values = my_colors) +
  ylim(c(4, 27)) +
  theme_base() +
  ylab("Lifespan (months)") +
  xlab("Depth (m)") + 
  theme(plot.background = element_blank(),
        legend.title = element_blank()) 

p
# ggsave("plots/statistics/lifespan/GAMM_lifespan_predict_depth_me.png", p, height = 40 , width = 40, units = "mm", scale = 4)

e2 <- p +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  xlab("") +
  ylab("")

# Zoo productivity:
p <- ggplot() +

  geom_smooth(data = stats_df4, aes(x = z_prod, y = pred_lifespan), color = "black", se = T, alpha = .4, size = .5) +        # Predicted regression
  
  geom_line(data = subset(partial_plot_lifespan, variable_id == "z_prod"), aes(x = variable, y = fit), color = "orange", size = 1) + 
  geom_ribbon(data = subset(partial_plot_lifespan, variable_id == "z_prod"), aes(y = fit, x = variable, ymin = (fit - se.fit), ymax =  (fit + se.fit)),
              fill = "orange", alpha = .1) +
  
  geom_pointrange(data = stats_df2, aes(x = z_prod, y = mean_lifespan_weighted_wt_cpue_cephalopoda_me,                               # Average in ME +- se
                                        ymin = mean_lifespan_weighted_wt_cpue_cephalopoda_me - se_lifespan_weighted_wt_cpue_cephalopoda_me,
                                        ymax = mean_lifespan_weighted_wt_cpue_cephalopoda_me + se_lifespan_weighted_wt_cpue_cephalopoda_me,
                                        color = coast)) +
  geom_pointrange(data = stats_df2, aes(x = z_prod, y = mean_lifespan_weighted_wt_cpue_cephalopoda_me, 
                                        xmin = z_prod - se_z_prod,
                                        xmax = z_prod + se_z_prod,
                                        color = coast)) +
  scale_color_manual(values = my_colors) +
  ylim(c(4, 27)) +
  theme_base() +
  ylab("Lifespan (months)") +
  xlab(expression(paste("Zooplankton productivity (g ", m^-2, " ", yr^-1, ")"))) + 
  theme(plot.background = element_blank(),
        legend.title = element_blank())

p
# ggsave("plots/statistics/lifespan/GAMM_lifespan_predict_z_prod_me.png", p, height = 40 , width = 40, units = "mm", scale = 4)

e3 <- p +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  xlab("") +
  ylab("")


#####################
# Cephalopod a growth:


# Make predictions of range of a variable keeping all other variables constant:
# Loop to obtain partial effect predictions for all the variables:
variables <- c("Temp_C", "depth", "z_prod")
partial_plot_a_growth <- data.frame()

test <- stats_df4 %>%
  dplyr::select(Temp_C, depth, z_prod, coast) %>%
  as.data.frame()

for(i in 1:length(variables)){
  var1 <- seq(min(test[,which(colnames(test) == variables[i])]),
              max(test[,which(colnames(test) == variables[i])]), length = nrow(test))
  
  vars <- variables[-i]
  
  var2 <- rep(mean(test[,which(colnames(test) == vars[1])]),length(var1)) 
  var3 <- rep(mean(test[,which(colnames(test) == vars[2])]),length(var1))
  
  coast <- rep("West coast", length(var1)) # set the random effect to the same coast
  
  newdat <- data.frame(var1, var2, var3, coast) 
  
  colnames(newdat)[1:3] <- c(variables[i], vars)
  
  my_predictions <- predict.gam(gam5, newdata = newdat, exclude = "s(coast)", type = "response", se.fit = TRUE) # predict with new data excluding random effect
  my_fit <- my_predictions$fit
  my_se <- my_predictions$se.fit
  
  partial <- data.frame(my_predictions, variable = var1) 
  partial <- partial %>% mutate(variable = var1,
                                variable_id = variables[i])
  
  partial_plot_a_growth = rbind(partial_plot_a_growth, partial)
  
}

# Temperature:
p <- ggplot() +
  
  geom_line(data = subset(partial_plot_a_growth, variable_id == "Temp_C"), aes(x = variable, y = fit), color = "orange", size = 1) + 
  geom_ribbon(data = subset(partial_plot_a_growth, variable_id == "Temp_C"), aes(y = fit, x = variable, ymin = (fit - se.fit), ymax =  (fit + se.fit)),
              fill = "orange", alpha = .1) +
  
  geom_smooth(data = stats_df4, aes(x = Temp_C, y = pred_a_growth), color = "black", se = T, alpha = .4, size = .5) +         # Predicted regression
  
  geom_pointrange(data = stats_df2, aes(x = Temp_C, y = mean_A_growth_weighted_wt_cpue_cephalopoda_me,                               # Average in ME +- se
                                        ymin = mean_A_growth_weighted_wt_cpue_cephalopoda_me - se_A_growth_weighted_wt_cpue_cephalopoda_me,
                                        ymax = mean_A_growth_weighted_wt_cpue_cephalopoda_me + se_A_growth_weighted_wt_cpue_cephalopoda_me,
                                        color = coast)) +
  geom_pointrange(data = stats_df2, aes(x = Temp_C, y = mean_A_growth_weighted_wt_cpue_cephalopoda_me, 
                                        xmin = Temp_C - se_Temp_C,
                                        xmax = Temp_C + se_Temp_C,
                                        color = coast)) +
  scale_color_manual(values = my_colors) +
  theme_base() +
  ylab("Growth coefficient A") +
  xlab("Temperature (ºC)") + 
  theme(plot.background = element_blank(),
        legend.title = element_blank())
# +
#   ylim(0, NA)

p
# ggsave("plots/statistics/a_growth/GAMM_a_growth_predict_Temp_C_me.png", p, height = 40 , width = 40, units = "mm", scale = 4)

f1 <- p +
  theme(legend.position = "none")  

# Depth:
p <- ggplot() +
  
  geom_line(data = subset(partial_plot_a_growth, variable_id == "depth"), aes(x = variable, y = fit), color = "orange", size = 1) + 
  geom_ribbon(data = subset(partial_plot_a_growth, variable_id == "depth"), aes(y = fit, x = variable, ymin = (fit - se.fit), ymax =  (fit + se.fit)),
              fill = "orange", alpha = .1) +
  
  geom_smooth(data = stats_df4, aes(x = depth, y = pred_a_growth), color = "black", se = T, alpha = .4, size = .5) +         # Predicted regression
  
  geom_pointrange(data = stats_df2, aes(x = depth, y = mean_A_growth_weighted_wt_cpue_cephalopoda_me,                               # Average in ME +- se
                                        ymin = mean_A_growth_weighted_wt_cpue_cephalopoda_me - se_A_growth_weighted_wt_cpue_cephalopoda_me,
                                        ymax = mean_A_growth_weighted_wt_cpue_cephalopoda_me + se_A_growth_weighted_wt_cpue_cephalopoda_me,
                                        color = coast)) +
  geom_pointrange(data = stats_df2, aes(x = depth, y = mean_A_growth_weighted_wt_cpue_cephalopoda_me, 
                                        xmin = depth - se_depth,
                                        xmax = depth + se_depth,
                                        color = coast)) +
  scale_color_manual(values = my_colors) +
  theme_base() +
  ylab("Growth coefficient A") +
  xlab("Depth (m)") + 
  theme(plot.background = element_blank(),
        legend.title = element_blank()) 

p
# ggsave("plots/statistics/a_growth/GAMM_a_growth_predict_depth_me.png", p, height = 40 , width = 40, units = "mm", scale = 4)

f2 <- p +
  theme(legend.position = "none",
        axis.text.y = element_blank()) +
  ylab("")

# Zoo productivity:
p <- ggplot() +
  
  geom_line(data = subset(partial_plot_a_growth, variable_id == "z_prod"), aes(x = variable, y = fit), color = "orange", size = 1) + 
  geom_ribbon(data = subset(partial_plot_a_growth, variable_id == "z_prod"), aes(y = fit, x = variable, ymin = (fit - se.fit), ymax =  (fit + se.fit)),
              fill = "orange", alpha = .1) +
  
  geom_smooth(data = stats_df4, aes(x = z_prod, y = pred_a_growth), color = "black", se = T, alpha = .4, size = .5) +        # Predicted regression
  
  geom_pointrange(data = stats_df2, aes(x = z_prod, y = mean_A_growth_weighted_wt_cpue_cephalopoda_me,                               # Average in ME +- se
                                        ymin = mean_A_growth_weighted_wt_cpue_cephalopoda_me - se_A_growth_weighted_wt_cpue_cephalopoda_me,
                                        ymax = mean_A_growth_weighted_wt_cpue_cephalopoda_me + se_A_growth_weighted_wt_cpue_cephalopoda_me,
                                        color = coast)) +
  geom_pointrange(data = stats_df2, aes(x = z_prod, y = mean_A_growth_weighted_wt_cpue_cephalopoda_me, 
                                        xmin = z_prod - se_z_prod,
                                        xmax = z_prod + se_z_prod,
                                        color = coast)) +
  scale_color_manual(values = my_colors) +
  theme_base() +
  ylab("Growth coefficient A") +
  xlab(expression(paste("Zoo. productivity (g ", m^-2, " ", yr^-1, ")"))) + 
  theme(#plot.background = element_blank(),
        legend.title = element_blank())

p
# ggsave("plots/statistics/a_growth/GAMM_a_growth_predict_z_prod_me.png", p, height = 40 , width = 40, units = "mm", scale = 4)

f3 <- p +
  theme(legend.position = "none",
        axis.text.y = element_blank()) +
  ylab("")

g <- p +
  theme(legend.position = "bottom") 

g <- cowplot::get_legend(g)


# Plot all dot hexagons:
# As a function of temperature:
p <- ggplot(data = stats_df4) +
  geom_point(aes(x = Temp_C, y = mean_A_growth_weighted_wt_cpue_cephalopoda, color = coast), shape = 21, alpha = .5)  + # Observed values
  geom_line(data = subset(partial_plot_a_growth, variable_id == "Temp_C"), aes(x = variable, y = fit), color = "orange", size = 1) + 
  geom_ribbon(data = subset(partial_plot_a_growth, variable_id == "Temp_C"), aes(y = fit, x = variable, ymin = (fit - se.fit), ymax = (fit + se.fit)),
              fill = "orange", alpha = .1) +
  # geom_point(aes(x = Temp_C, y = pred_proportion, color = coast)) +                                      # Predicted values
  geom_smooth(aes(x = Temp_C, y = pred_a_growth), color = "black", se = F) +                           # Predicted regression
  theme_classic()

p
# ggsave("plots/statistics/a_growth/GAMM_a_growth_predict_Temp_C.png", p, height = 40 , width = 40, units = "mm", scale = 4)

# As a function of depth:
p <- ggplot(data = stats_df4) +
  geom_point(aes(x = depth, y = mean_A_growth_weighted_wt_cpue_cephalopoda, color = coast), shape = 21, alpha = .5)  +  # Observed values
  geom_line(data = subset(partial_plot_a_growth, variable_id == "depth"), aes(x = variable, y = fit), color = "orange", size = 1) + 
  geom_ribbon(data = subset(partial_plot_a_growth, variable_id == "depth"), aes(y = fit, x = variable, ymin = (fit - se.fit), ymax = (fit + se.fit)),
              fill = "orange", alpha = .1) +
  # geom_point(aes(x = depth, y = pred_proportion, color = coast)) +                                       # Predicted values
  geom_smooth(aes(x = depth, y = pred_a_growth), color = "black", se = F) +                            # Predicted regression
  theme_classic() 
p
# ggsave("plots/statistics/a_growth/GAMM_a_growth_predict_depth.png", p, height = 40 , width = 40, units = "mm", scale = 4)

# As a function of z_prod:
p <- ggplot(data = stats_df4) +
  geom_point(aes(x = z_prod, y = mean_A_growth_weighted_wt_cpue_cephalopoda, color = coast), shape = 21, alpha = .5)  +  # Observed values
  geom_line(data = subset(partial_plot_a_growth, variable_id == "z_prod"), aes(x = variable, y = fit), color = "orange", size = 1) + 
  geom_ribbon(data = subset(partial_plot_a_growth, variable_id == "z_prod"), aes(y = fit, x = variable, ymin = (fit - se.fit), ymax = (fit + se.fit)),
              fill = "orange", alpha = .1) +
  # geom_point(aes(x = z_prod, y = pred_proportion, color = coast)) +                                       # Predicted values
  geom_smooth(aes(x = z_prod, y = pred_a_growth), color = "black", se = F) +                            # Predicted regression
  theme_classic() 

p
# ggsave("plots/statistics/a_growth/GAMM_a_growth_predict_z_prod.png", p, height = 40 , width = 40, units = "mm", scale = 4)


# Combine all smooths in a plot:
p <- (a1 | a2 | a3) /
  (b1 | b2 | b3) / 
  (c1 | c2 | c3) /
  (d1 | d2 | d3) /
  (e1 | e2 | e3) /
  (f1 | f2 | f3) +
  plot_annotation(tag_levels = 'a') +
  theme(plot.tag = element_text(size = 20)) & 
  theme(plot.background = element_blank())

p
ggsave("plots/statistics/regression_plots.png", p, height = 100 , width = 66, units = "mm", scale = 4)

p <- (a1 | a3) /
  (b1 | b3) / 
  (c1 | c3) /
  (d1 | d3) /
  (e1 | e3) /
  (f1 | f3) +
  plot_annotation(tag_levels = 'a') +
  theme(plot.tag = element_text(size = 20)) & 
  theme(plot.background = element_blank())

p
ggsave("plots/statistics/regression_plots2.png", p, height = 100 , width = 45, units = "mm", scale = 4)


rm(list = c("a1", "a2", "a3", "b1", "b2", "b3", "c1", "c2", "c3", "d1", "d2", "d3", "e1", "e2", "e3", "f1", "f2", "f3"))


# # Combine all partial effect data frames:
# 
# partial_plots_df <- rbind(partial_plot_proportion, 
#                           partial_plot_ceph_bio, 
#                           partial_plot_fish_bio, 
#                           partial_plot_a_size, 
#                           partial_plot_lifespan, 
#                           partial_plot_a_growth)
# 
# # save(partial_plots_df, file = "data/partial_plots_df.Rdata")
# 
# 
# rm(list = c("partial_plot_proportion", "partial_plot_ceph_bio", "partial_plot_fish_bio", 
#              "partial_plot_a_size", "partial_plot_lifespan", "partial_plot_a_growth"))


#-------------------------------------------------------------------------
#
# Random Forest -----------------------------------------------------------
#

# Create an empty list to store the plots:
rf_list <- list()

# Convert proportion to %:
stats_df$cephalopoda_fish_proportion <- 100 * stats_df$cephalopoda_fish_proportion 

# Cephalopod proportion:
rf0 <- randomForest(cephalopoda_fish_proportion ~ Temp_C + depth + z_prod + coast, #+ mean_wt_cpue_fish
                   data = stats_df,
                   ntree = 1000,
                   importance = T,
                   mtry = 2,
                   na.action = na.omit)

rf0

# Partial plots:
vars <- c("Temp_C", "depth", "z_prod", "coast")

for(i in 1:length(vars)) {
  p <- rf0 %>%  
    pdp::partial(pred.var = vars[i]) %>%
    autoplot(smooth = TRUE) + 
    ylab("Cephalopod proportion (%)") +
    # scale_y_continuous(position = "right") +
    theme_base() +
    ylim(0, 10)
    
  
  rf_list[[i]] <- p

  ggsave(paste0("plots/statistics/ceph_proportion/random_forest/ceph_proportion_", vars[i],".png"), p, height = 40 , width = 40, units = "mm", scale = 4)
}

p1 <- rf_list[[1]]
p2 <- rf_list[[2]]
p3 <- rf_list[[3]]

# Cephalopod biomass:
rf1 <- randomForest(mean_wt_cpue_cephalopoda ~ Temp_C + depth + z_prod + coast,
                    data = stats_df,
                    ntree = 1000,
                    importance = T,
                    mtry = 2,
                    na.action = na.omit)

rf1


for(i in 1:length(vars)) {
  p <- rf1 %>%  
    pdp::partial(pred.var = vars[i]) %>%
    autoplot(smooth = TRUE) + 
    ylab("Cephalopod biomass (kg km^-2)") +
    # scale_y_continuous(position = "right") +
    theme_base() +
    scale_y_log10(limits = c(80, 600))
  
  rf_list[[i]] <- p
  
  ggsave(paste0("plots/statistics/ceph_biomass/random_forest/ceph_biomass_", vars[i],".png"), p, height = 40 , width = 40, units = "mm", scale = 4)
}

p4 <- rf_list[[1]]
p5 <- rf_list[[2]]
p6 <- rf_list[[3]]

# Fish biomass:
rf2 <- randomForest(mean_wt_cpue_fish ~ Temp_C + depth + z_prod + coast,
                    data = stats_df,
                    ntree = 1000,
                    importance = T,
                    mtry = 2,
                    na.action = na.omit)

rf2


for(i in 1:length(vars)) {
  p <- rf2 %>%  
    pdp::partial(pred.var = vars[i]) %>%
    autoplot(smooth = TRUE) + 
    ylab("Fish biomass (kg km^-2)") +
    # scale_y_continuous(position = "right") +
    theme_base()  +
    scale_y_log10(limits = c(3.5*10^4, 1.5*10^5))
  
  rf_list[[i]] <- p
  
  ggsave(paste0("plots/statistics/fish_biomass/random_forest/fish_biomass_", vars[i],".png"), p, height = 40 , width = 40, units = "mm", scale = 4)
}

p7 <- rf_list[[1]]
p8 <- rf_list[[2]]
p9 <- rf_list[[3]]

# Asymptotic size:
rf3 <- randomForest(mean_a_size_weighted_wt_cpue_cephalopoda ~ Temp_C + depth + z_prod + coast,
                    data = stats_df4,
                    ntree = 1000,
                    importance = T,
                    mtry = 2,
                    na.action = na.omit)
rf3

# Partial plots:
for(i in 1:length(vars)) {
  p <- rf3 %>%  
    pdp::partial(pred.var = vars[i]) %>%
    autoplot(smooth = TRUE) + 
    ylab("Asymptotic size (g)") +
    # scale_y_continuous(position = "right") +
    theme_base()  +
    scale_y_log10(limits = c(1500, 5500))
  
  rf_list[[i]] <- p
  
  ggsave(paste0("plots/statistics/a_size/random_forest/a_size_", vars[i],".png"), p, height = 40 , width = 40, units = "mm", scale = 4)
}

p10 <- rf_list[[1]]
p11 <- rf_list[[2]]
p12 <- rf_list[[3]]


# Lifespan:
rf4 <- randomForest(mean_lifespan_weighted_wt_cpue_cephalopoda ~ Temp_C + depth + z_prod + coast,
                    data = stats_df4,
                    ntree = 1000,
                    importance = T,
                    mtry = 2,
                    na.action = na.omit)

rf4

# Partial plots:
for(i in 1:length(vars)) {
  p <- rf4 %>%  
    pdp::partial(pred.var = vars[i]) %>%
    autoplot(smooth = TRUE) + ylab("Lifespan (months)") +
    # scale_y_continuous(position = "right") +
    theme_base() +
    ylim(c(12, 17))
  
  rf_list[[i]] <- p
  
  ggsave(paste0("plots/statistics/lifespan/random_forest/lifespan_", vars[i],".png"), p, height = 40 , width = 40, units = "mm", scale = 4)
}

p13 <- rf_list[[1]]
p14 <- rf_list[[2]]
p15 <- rf_list[[3]]


# A growth:
rf5 <- randomForest(mean_A_growth_weighted_wt_cpue_cephalopoda ~ Temp_C + depth + z_prod + coast,
                    data = stats_df4,
                    ntree = 1000,
                    importance = T,
                    mtry = 2,
                    na.action = na.omit)


rf5

# Partial plots:
for(i in 1:length(vars)) {
  p <- rf5 %>%  
    pdp::partial(pred.var = vars[i]) %>%
    autoplot(smooth = TRUE) + ylab("Growth parameter A") +
    # scale_y_continuous(position = "right") +
    theme_base() 
  rf_list[[i]] <- p

  ggsave(paste0("plots/statistics/a_growth/random_forest/a_growth_", vars[i],".png"), p, height = 40 , width = 40, units = "mm", scale = 4)
}

p16 <- rf_list[[1]]
p17 <- rf_list[[2]]
p18 <- rf_list[[3]]


# Combine all smooths in a plot:
p <- (p4 | p5 | p6) /
  (p7 | p8 | p9) / 
  (p1 | p2 | p3) /
  (p10 | p11 | p12) /
  (p13 | p14 | p15) /
  (p16 | p17 | p18) +
  plot_annotation(tag_levels = 'a') +
  theme(plot.tag = element_text(size = 20)) & 
  theme(plot.background = element_blank())

p
# ggsave("plots/statistics/random_forest.png", p, height = 105 , width = 66, units = "mm", scale = 5)



######################################################################################
#-------------------------------------------------------------------------------------


#
# Power analysis
# Generate a random sequence and split data into four parts:
#

set.seed(123)
n_replicates <- 100

# Ceph proportion:
funRand <- function(){
  # Test and training data 
  d1 <- stats_df %>% sample_n(max(0.66 * n(), 1))
  p1 <- stats_df %>%  anti_join(d1)
  
  # Fit and predict
  fit1 <- gam(cephalopoda_fish_proportion ~ 
                s(Temp_C, k = 3) +
                s(depth, k = 3) +
                s(z_prod, k = 3) + 
                # s(mean_wt_cpue_fish, k = 3) +
                s(coast, bs = "re"),
              data = d1,
              method = "REML",
              family = betar())
  
  f1 <- summary(fit1)$dev.expl

    fit2 <- randomForest(cephalopoda_fish_proportion ~ Temp_C + depth + z_prod + coast,
                       data = d1,
                       ntree = 1000,
                       importance = T,
                       mtry = 2,
                       na.action = na.omit)
  
  f2 <- mean(fit2$rsq)
  
  # Predict the remaining part and compare with observations
  Pred1 <- predict.gam(fit1, newdata = as.data.frame(p1), type = "response")
  Pred2 <- predict(fit2, newdata = as.data.frame(p1))

  # Mean squared errors
  R1 <- (p1[ , 13] - Pred1)^2
  R1 <- sum(R1$cephalopoda_fish_proportion)/(nrow(R1) - 1)
  R2 <- (p1[ , 13] - Pred2)^2 
  R2 <- sum(R2$cephalopoda_fish_proportion)/(nrow(R2) - 1)
  
  return(c(f1, R1, f2, R2))
}
funRand()

# Return results for runs 
retlistS <- t(unlist(replicate(n_replicates, funRand())))
colnames(retlistS) <- c("GAM","GAM","RF","RF") 
summary(retlistS)
comparison <- as.data.frame(retlistS)
expl_var <- reshape2::melt(comparison[ , c(1, 3)])
MSE <- reshape2::melt(comparison[,c(2,4)])

ggplot() + geom_boxplot(data = expl_var, aes(x = variable, y = value)) + labs(y = "Variance", x = "") + theme_bw() +
  ggplot() + geom_boxplot(data = MSE, aes(x = variable, y = value)) + labs(y = "MSE", x = "") + theme_bw() 

# Check average output:
expl_var %>%
  group_by(variable) %>%
  summarise(mean = mean(value))

MSE %>%
  group_by(variable) %>%
  summarise(mean = mean(value))


# Try to predict one region based on the other two regions:
funRandCoast <- function(coast_name){
  # Test and training data 
  d1 <- stats_df %>% filter(coast != coast_name)
  p1 <- stats_df %>%  filter(coast == coast_name)
  
  # Fit and predict
  fit1 <- gam(cephalopoda_fish_proportion ~ 
                s(Temp_C, k = 3) +
                s(depth, k = 3) +
                s(z_prod, k = 3) + 
                # s(mean_wt_cpue_fish, k = 3) +
                s(coast, bs = "re"),
              data = d1,
              method = "REML",
              family = betar())
  
  f1 <- summary(fit1)$dev.expl
  
  fit2 <- randomForest(cephalopoda_fish_proportion ~ Temp_C + depth + z_prod + coast,
                       data = d1,
                       ntree = 1000,
                       importance = T,
                       mtry = 2,
                       na.action = na.omit)
  
  f2 <- mean(fit2$rsq)
  
  # Predict the remaining part and compare with observations
  Pred1 <- predict.gam(fit1, newdata = as.data.frame(p1), type = "response")
  Pred2 <- predict(fit2, newdata = as.data.frame(p1))
  
  
  R1 <- (p1[ , 13] - Pred1)^2
  R1 <- sum(R1$cephalopoda_fish_proportion)/(nrow(R1) - 1)
  R2 <- (p1[ , 13] - Pred2)^2 
  R2 <- sum(R2$cephalopoda_fish_proportion)/(nrow(R2) - 1)
  
  return(c(f1, R1, f2, R2))
}

funRandCoast("Europe")

coast_name <- as.vector(unique(stats_df$coast))

retlistS <- c()

for(i in 1:length(coast_name)) {
  retlistS <- rbind(retlistS, funRandCoast(coast_name[i]))
}


colnames(retlistS) <- c("GAM2","GAM2","RF2","RF2") 
summary(retlistS)
comparison <- as.data.frame(retlistS)
expl_var2 <- reshape2::melt(comparison[ , c(1, 3)])
MSE2 <- reshape2::melt(comparison[,c(2,4)])

# Check average output:
expl_var2 %>%
  group_by(variable) %>%
  summarise(mean = mean(value))

MSE2 %>%
  group_by(variable) %>%
  summarise(mean = mean(value))

expl_var3 <- rbind(expl_var, expl_var2)
MSE3 <- rbind(MSE, MSE2)


p <- ggplot() + 
  geom_boxplot(data = expl_var, aes(x = variable, y = value)) + 
  labs(y = "Variance", x = "") + 
  theme_base() +
  theme(plot.background = element_blank()) +
  ggtitle("Random bootsrap") +
  
  ggplot() + 
  geom_boxplot(data = MSE, aes(x = variable, y = value)) + 
  labs(y = "MSE", x = "") + 
  theme_base() +
  theme(plot.background = element_blank()) +
  ggtitle("Random bootsrap") +
  
  ggplot() + 
  geom_boxplot(data = expl_var2, aes(x = variable, y = value)) + 
  labs(y = "Variance", x = "") + 
  theme_base() + 
  theme(plot.background = element_blank()) +
  ggtitle("Regional bootsrap") +
  
  ggplot() + 
  geom_boxplot(data = MSE2, aes(x = variable, y = value)) + 
  labs(y = "MSE", x = "") + 
  theme_base() +
  theme(plot.background = element_blank()) +
  ggtitle("Regional bootsrap") 

p

ggsave("plots/statistics/ceph_proportion/power_analysis/ceph_proportion.png", p, height = 50 , width = 50, units = "mm", scale = 4)


# Run some t-tests:
# MSE
GAM <- subset(MSE3, variable == "GAM")[2]
RF <- subset(MSE3, variable == "RF")[2]
t.test(GAM, RF) # p-value = 5.642e-06
t1 <- t.test(GAM, RF)$p.val

GAM2 <- subset(MSE3, variable == "GAM2")[2]
RF2 <- subset(MSE3, variable == "RF2")[2]
t.test(GAM2, RF2) # p-value = 0.5212
t2 <- t.test(GAM2, RF2)$p.val

# Variance
GAM <- subset(expl_var3, variable == "GAM")[2]
RF <- subset(expl_var3, variable == "RF")[2]
t.test(GAM, RF) # p-value < 5.118065e-61
t3 <- t.test(GAM, RF)$p.val

GAM2 <- subset(expl_var3, variable == "GAM2")[2]
RF2 <- subset(expl_var3, variable == "RF2")[2]
t.test(GAM2, RF2) # p-value = 0.2158
t4 <- t.test(GAM2, RF2)$p.val


# Within GAM:
gam_random_var <- subset(expl_var3, variable == "GAM")[2]
gam_region_var <- subset(expl_var3, variable == "GAM2")[2]
t.test(gam_random_var, gam_region_var) # p-value = 0.5145
t5 <- t.test(gam_random_var, gam_region_var)$p.val

gam_random_mse <- subset(MSE3, variable == "GAM")[2]
gam_region_mse <- subset(MSE3, variable == "GAM2")[2]
t.test(gam_random_mse, gam_region_mse) # p-value = 0.7068
t6 <- t.test(gam_random_mse, gam_region_mse)$p.val


#--------------------------------------------------------------

# Ceph biomass:
funRand <- function(){
  # Test and training data 
  d1 <- stats_df %>% sample_n(max(0.66 * n(), 1))
  p1 <- stats_df %>%  anti_join(d1)
  
  # Fit and predict
  fit1 <- gam(mean_wt_cpue_cephalopoda ~ 
                s(Temp_C, k = 3) +
                s(depth, k = 3) +
                s(z_prod, k = 3) + 
                # s(mean_wt_cpue_fish, k = 3) +
                s(coast, bs = "re"),
              data = d1,
              method = "REML",
              family = nb())
  
  f1 <- summary(fit1)$dev.expl
  
  fit2 <- randomForest(mean_wt_cpue_cephalopoda ~ Temp_C + depth + z_prod + coast,
                       data = d1,
                       ntree = 1000,
                       importance = T,
                       mtry = 2,
                       na.action = na.omit)
  
  f2 <- mean(fit2$rsq)
  
  # Predict the remaining part and compare with observations
  Pred1 <- predict.gam(fit1, newdata = as.data.frame(p1), type = "response")
  Pred2 <- predict(fit2, newdata = as.data.frame(p1))
  
  # Mean squared errors
  R1 <- (p1[ , 6] - Pred1)^2
  R1 <- sum(R1$mean_wt_cpue_cephalopoda)/(nrow(R1) - 1)
  R2 <- (p1[ , 6] - Pred2)^2 
  R2 <- sum(R2$mean_wt_cpue_cephalopoda)/(nrow(R2) - 1)
  
  return(c(f1, R1, f2, R2))
}
funRand()

# Return results for runs 
retlistS <- t(unlist(replicate(n_replicates, funRand())))
colnames(retlistS) <- c("GAM","GAM","RF","RF") 
summary(retlistS)
comparison <- as.data.frame(retlistS)
expl_var <- reshape2::melt(comparison[ , c(1, 3)])
MSE <- reshape2::melt(comparison[,c(2,4)])

ggplot() + geom_boxplot(data = expl_var, aes(x = variable, y = value)) + labs(y = "Variance", x = "") + theme_bw() +
  ggplot() + geom_boxplot(data = MSE, aes(x = variable, y = value)) + labs(y = "MSE", x = "") + theme_bw() 


# Check average output:
expl_var %>%
  group_by(variable) %>%
  summarise(mean = mean(value))

MSE %>%
  group_by(variable) %>%
  summarise(mean = mean(value))


# Try to predict one region based on the other two regions:
funRandCoast <- function(coast_name){
  # Test and training data 
  d1 <- stats_df %>% filter(coast != coast_name)
  p1 <- stats_df %>%  filter(coast == coast_name)
  
  # Fit and predict
  fit1 <- gam(mean_wt_cpue_cephalopoda ~ 
                s(Temp_C, k = 3) +
                s(depth, k = 3) +
                s(z_prod, k = 3) + 
                # s(mean_wt_cpue_fish, k = 3) +
                s(coast, bs = "re"),
              data = d1,
              method = "REML",
              family = nb())
  
  f1 <- summary(fit1)$dev.expl
  
  fit2 <- randomForest(mean_wt_cpue_cephalopoda ~ Temp_C + depth + z_prod + coast,
                       data = d1,
                       ntree = 1000,
                       importance = T,
                       mtry = 2,
                       na.action = na.omit)
  
  f2 <- mean(fit2$rsq)
  
  # Predict the remaining part and compare with observations
  Pred1 <- predict.gam(fit1, newdata = as.data.frame(p1), type = "response")
  Pred2 <- predict(fit2, newdata = as.data.frame(p1))
  
  
  R1 <- (p1[ , 6] - Pred1)^2
  R1 <- sum(R1$mean_wt_cpue_cephalopoda)/(nrow(R1) - 1)
  R2 <- (p1[ , 6] - Pred2)^2 
  R2 <- sum(R2$mean_wt_cpue_cephalopoda)/(nrow(R2) - 1)
  
  return(c(f1, R1, f2, R2))
}

funRandCoast("Europe")

coast_name <- as.vector(unique(stats_df$coast))

retlistS <- c()

for(i in 1:length(coast_name)) {
  retlistS <- rbind(retlistS, funRandCoast(coast_name[i]))
}


colnames(retlistS) <- c("GAM2","GAM2","RF2","RF2") 
summary(retlistS)
comparison <- as.data.frame(retlistS)
expl_var2 <- reshape2::melt(comparison[ , c(1, 3)])
MSE2 <- reshape2::melt(comparison[,c(2,4)])

# Check average output:
expl_var2 %>%
  group_by(variable) %>%
  summarise(mean = mean(value))

MSE2 %>%
  group_by(variable) %>%
  summarise(mean = mean(value))


expl_var3 <- rbind(expl_var, expl_var2)
MSE3 <- rbind(MSE, MSE2)


p <- ggplot() + 
  geom_boxplot(data = expl_var, aes(x = variable, y = value)) + 
  labs(y = "Variance", x = "") + 
  theme_base() +
  theme(plot.background = element_blank()) +
  ggtitle("Random bootsrap") +
  
  ggplot() + 
  geom_boxplot(data = MSE, aes(x = variable, y = value)) + 
  labs(y = "MSE", x = "") + 
  theme_base() +
  theme(plot.background = element_blank()) +
  ggtitle("Random bootsrap") +
  
  ggplot() + 
  geom_boxplot(data = expl_var2, aes(x = variable, y = value)) + 
  labs(y = "Variance", x = "") + 
  theme_base() + 
  theme(plot.background = element_blank()) +
  ggtitle("Regional bootsrap") +
  
  ggplot() + 
  geom_boxplot(data = MSE2, aes(x = variable, y = value)) + 
  labs(y = "MSE", x = "") + 
  theme_base() +
  theme(plot.background = element_blank()) +
  ggtitle("Regional bootsrap") 

p

ggsave("plots/statistics/ceph_biomass/power_analysis/ceph_biomass.png", p, height = 50 , width = 50, units = "mm", scale = 4)


# Run some t-tests:
# MSE
GAM <- subset(MSE3, variable == "GAM")[2]
RF <- subset(MSE3, variable == "RF")[2]
t.test(GAM, RF) # p-value = 0.2160759
t7 <- t.test(GAM, RF)$p.val

GAM2 <- subset(MSE3, variable == "GAM2")[2]
RF2 <- subset(MSE3, variable == "RF2")[2]
t.test(GAM2, RF2) # p-value = 0.4226151
t8 <- t.test(GAM2, RF2)$p.val

# Variance
GAM <- subset(expl_var3, variable == "GAM")[2]
RF <- subset(expl_var3, variable == "RF")[2]
t.test(GAM, RF) # p-value = 4.161496e-67
t9 <- t.test(GAM, RF)$p.val

GAM2 <- subset(expl_var3, variable == "GAM2")[2]
RF2 <- subset(expl_var3, variable == "RF2")[2]
t.test(GAM2, RF2) # p-value = 0.06268907
t10 <- t.test(GAM2, RF2)$p.val


# Within GAM:
gam_random_var <- subset(expl_var3, variable == "GAM")[2]
gam_region_var <- subset(expl_var3, variable == "GAM2")[2]
t11 <- t.test(gam_random_var, gam_region_var)$p.val # p-value = 0.7463

gam_random_mse <- subset(MSE3, variable == "GAM")[2]
gam_region_mse <- subset(MSE3, variable == "GAM2")[2]
t12 <- t.test(gam_random_mse, gam_region_mse)$p.val # p-value = 0.4226


#--------------------------------------------------------------

# Fish biomass:
funRand <- function(){
  # Test and training data 
  d1 <- stats_df %>% sample_n(max(0.66 * n(), 1))
  p1 <- stats_df %>%  anti_join(d1)
  
  # Fit and predict
  fit1 <- gam(mean_wt_cpue_fish ~ 
                s(Temp_C, k = 3) +
                s(depth, k = 3) +
                s(z_prod, k = 3) + 
                # s(mean_wt_cpue_fish, k = 3) +
                s(coast, bs = "re"),
              data = d1,
              method = "REML",
              family = nb())
  
  f1 <- summary(fit1)$dev.expl
  
  fit2 <- randomForest(mean_wt_cpue_fish ~ Temp_C + depth + z_prod + coast,
                       data = d1,
                       ntree = 1000,
                       importance = T,
                       mtry = 2,
                       na.action = na.omit)
  
  f2 <- mean(fit2$rsq)
  
  # Predict the remaining part and compare with observations
  Pred1 <- predict.gam(fit1, newdata = as.data.frame(p1), type = "response")
  Pred2 <- predict(fit2, newdata = as.data.frame(p1))
  
  # Mean squared errors
  R1 <- (p1[ , 7] - Pred1)^2
  R1 <- sum(R1$mean_wt_cpue_fish)/(nrow(R1) - 1)
  R2 <- (p1[ , 7] - Pred2)^2 
  R2 <- sum(R2$mean_wt_cpue_fish)/(nrow(R2) - 1)
  
  return(c(f1, R1, f2, R2))
}
funRand()

# Return results for runs 
retlistS <- t(unlist(replicate(n_replicates, funRand())))
colnames(retlistS) <- c("GAM","GAM","RF","RF") 
summary(retlistS)
comparison <- as.data.frame(retlistS)
expl_var <- reshape2::melt(comparison[ , c(1, 3)])
MSE <- reshape2::melt(comparison[,c(2,4)])

ggplot() + geom_boxplot(data = expl_var, aes(x = variable, y = value)) + labs(y = "Variance", x = "") + theme_bw() +
  ggplot() + geom_boxplot(data = MSE, aes(x = variable, y = value)) + labs(y = "MSE", x = "") + theme_bw() 


# Check average output:
expl_var %>%
  group_by(variable) %>%
  summarise(mean = mean(value))

MSE %>%
  group_by(variable) %>%
  summarise(mean = mean(value))


# Try to predict one region based on the other two regions:
funRandCoast <- function(coast_name){
  # Test and training data 
  d1 <- stats_df %>% filter(coast != coast_name)
  p1 <- stats_df %>%  filter(coast == coast_name)
  
  # Fit and predict
  fit1 <- gam(mean_wt_cpue_fish ~ 
                s(Temp_C, k = 3) +
                s(depth, k = 3) +
                s(z_prod, k = 3) + 
                # s(mean_wt_cpue_fish, k = 3) +
                s(coast, bs = "re"),
              data = d1,
              method = "REML",
              family = nb())
  
  f1 <- summary(fit1)$dev.expl
  
  fit2 <- randomForest(mean_wt_cpue_fish ~ Temp_C + depth + z_prod + coast,
                       data = d1,
                       ntree = 1000,
                       importance = T,
                       mtry = 2,
                       na.action = na.omit)
  
  f2 <- mean(fit2$rsq)
  
  # Predict the remaining part and compare with observations
  Pred1 <- predict.gam(fit1, newdata = as.data.frame(p1), type = "response")
  Pred2 <- predict(fit2, newdata = as.data.frame(p1))
  
  
  R1 <- (p1[ , 7] - Pred1)^2
  R1 <- sum(R1$mean_wt_cpue_fish)/(nrow(R1) - 1)
  R2 <- (p1[ , 7] - Pred2)^2 
  R2 <- sum(R2$mean_wt_cpue_fish)/(nrow(R2) - 1)
  
  return(c(f1, R1, f2, R2))
}

funRandCoast("Europe")

coast_name <- as.vector(unique(stats_df$coast))

retlistS <- c()

for(i in 1:length(coast_name)) {
  retlistS <- rbind(retlistS, funRandCoast(coast_name[i]))
}


colnames(retlistS) <- c("GAM2","GAM2","RF2","RF2") 
summary(retlistS)
comparison <- as.data.frame(retlistS)
expl_var2 <- reshape2::melt(comparison[ , c(1, 3)])
MSE2 <- reshape2::melt(comparison[,c(2,4)])

# Check average output:
expl_var2 %>%
  group_by(variable) %>%
  summarise(mean = mean(value))

MSE2 %>%
  group_by(variable) %>%
  summarise(mean = mean(value))


expl_var3 <- rbind(expl_var, expl_var2)
MSE3 <- rbind(MSE, MSE2)


p <- ggplot() + 
  geom_boxplot(data = expl_var, aes(x = variable, y = value)) + 
  labs(y = "Variance", x = "") + 
  theme_base() +
  theme(plot.background = element_blank()) +
  ggtitle("Random bootsrap") +
  
  ggplot() + 
  geom_boxplot(data = MSE, aes(x = variable, y = value)) + 
  labs(y = "MSE", x = "") + 
  theme_base() +
  theme(plot.background = element_blank()) +
  ggtitle("Random bootsrap") +
  
  ggplot() + 
  geom_boxplot(data = expl_var2, aes(x = variable, y = value)) + 
  labs(y = "Variance", x = "") + 
  theme_base() + 
  theme(plot.background = element_blank()) +
  ggtitle("Regional bootsrap") +
  
  ggplot() + 
  geom_boxplot(data = MSE2, aes(x = variable, y = value)) + 
  labs(y = "MSE", x = "") + 
  theme_base() +
  theme(plot.background = element_blank()) +
  ggtitle("Regional bootsrap") 

p

ggsave("plots/statistics/fish_biomass/power_analysis/fish_biomass.png", p, height = 50 , width = 50, units = "mm", scale = 4)


# Run some t-tests:
# MSE
GAM <- subset(MSE3, variable == "GAM")[2]
RF <- subset(MSE3, variable == "RF")[2]
t.test(GAM, RF) # p-value = 0.8953138
t13 <- t.test(GAM, RF)$p.val

GAM2 <- subset(MSE3, variable == "GAM2")[2]
RF2 <- subset(MSE3, variable == "RF2")[2]
t.test(GAM2, RF2) # p-value = 0.1793559
t14 <- t.test(GAM2, RF2)$p.val

# Variance
GAM <- subset(expl_var3, variable == "GAM")[2]
RF <- subset(expl_var3, variable == "RF")[2]
t.test(GAM, RF) # p-value = 1.700959e-79
t15 <- t.test(GAM, RF)$p.val

GAM2 <- subset(expl_var3, variable == "GAM2")[2]
RF2 <- subset(expl_var3, variable == "RF2")[2]
t.test(GAM2, RF2) # p-value = 0.3153718
t16 <- t.test(GAM2, RF2)$p.val


# Within GAM:
gam_random_var <- subset(expl_var3, variable == "GAM")[2]
gam_region_var <- subset(expl_var3, variable == "GAM2")[2]
t17 <- t.test(gam_random_var, gam_region_var)$p.val # p-value = 0.8867

gam_random_mse <- subset(MSE3, variable == "GAM")[2]
gam_region_mse <- subset(MSE3, variable == "GAM2")[2]
t18 <- t.test(gam_random_mse, gam_region_mse)$p.val # p-value = 0.9562


#--------------------------------------------------------------

# Asymptotic size:
funRand <- function(){
  # Test and training data 
  d1 <- stats_df4 %>% sample_n(max(0.66 * n(), 1))
  p1 <- stats_df4 %>%  anti_join(d1)
  
  # Fit and predict
  fit1 <- gam(log10(mean_a_size_weighted_wt_cpue_cephalopoda) ~ 
                s(Temp_C, k = 3) +
                s(depth, k = 3) +
                s(z_prod, k = 3) + 
                # s(mean_wt_cpue_fish, k = 3) +
                s(coast, bs = "re"),
              data = d1,
              method = "REML",
              family = gaussian())
  
  f1 <- summary(fit1)$dev.expl
  
  fit2 <- randomForest(log10(mean_a_size_weighted_wt_cpue_cephalopoda) ~ Temp_C + depth + z_prod + coast,
                       data = d1,
                       ntree = 1000,
                       importance = T,
                       mtry = 2,
                       na.action = na.omit)
  
  f2 <- mean(fit2$rsq)
  
  # Predict the remaining part and compare with observations
  Pred1 <- 10^(predict.gam(fit1, newdata = as.data.frame(p1), type = "response"))
  Pred2 <- 10^(predict(fit2, newdata = as.data.frame(p1)))
  
  # Mean squared errors
  R1 <- (p1[ , 10] - Pred1)^2
  R1 <- sum(R1$mean_a_size_weighted_wt_cpue_cephalopoda)/(nrow(R1) - 1)
  R2 <- (p1[ , 10] - Pred2)^2 
  R2 <- sum(R2$mean_a_size_weighted_wt_cpue_cephalopoda)/(nrow(R2) - 1)
  
  return(c(f1, R1, f2, R2))
}
funRand()

# Return results for runs 
retlistS <- t(unlist(replicate(n_replicates, funRand())))
colnames(retlistS) <- c("GAM","GAM","RF","RF") 
summary(retlistS)
comparison <- as.data.frame(retlistS)
expl_var <- reshape2::melt(comparison[ , c(1, 3)])
MSE <- reshape2::melt(comparison[,c(2,4)])

ggplot() + geom_boxplot(data = expl_var, aes(x = variable, y = value)) + labs(y = "Variance", x = "") + theme_bw() +
  ggplot() + geom_boxplot(data = MSE, aes(x = variable, y = value)) + labs(y = "MSE", x = "") + theme_bw() 


# Check average output:
expl_var %>%
  group_by(variable) %>%
  summarise(mean = mean(value))

MSE %>%
  group_by(variable) %>%
  summarise(mean = mean(value))


# Try to predict one region based on the other two regions:
funRandCoast <- function(coast_name){
  # Test and training data 
  d1 <- stats_df4 %>% filter(coast != coast_name)
  p1 <- stats_df4 %>%  filter(coast == coast_name)
  
  # Fit and predict
  fit1 <- gam(log10(mean_a_size_weighted_wt_cpue_cephalopoda) ~ 
                s(Temp_C, k = 3) +
                s(depth, k = 3) +
                s(z_prod, k = 3) + 
                # s(mean_wt_cpue_fish, k = 3) +
                s(coast, bs = "re"),
              data = d1,
              method = "REML",
              family = gaussian())
  
  f1 <- summary(fit1)$dev.expl
  
  fit2 <- randomForest(log10(mean_a_size_weighted_wt_cpue_cephalopoda) ~ Temp_C + depth + z_prod + coast,
                       data = d1,
                       ntree = 1000,
                       importance = T,
                       mtry = 2,
                       na.action = na.omit)
  
  f2 <- mean(fit2$rsq)
  
  # Predict the remaining part and compare with observations
  Pred1 <- 10^(predict.gam(fit1, newdata = as.data.frame(p1), type = "response"))
  Pred2 <- 10^(predict(fit2, newdata = as.data.frame(p1)))
  
  
  R1 <- (p1[ , 10] - Pred1)^2
  R1 <- sum(R1$mean_a_size_weighted_wt_cpue_cephalopoda)/(nrow(R1) - 1)
  R2 <- (p1[ , 10] - Pred2)^2 
  R2 <- sum(R2$mean_a_size_weighted_wt_cpue_cephalopoda)/(nrow(R2) - 1)
  
  return(c(f1, R1, f2, R2))
}


funRandCoast("Europe")

coast_name <- as.vector(unique(stats_df$coast))

retlistS <- c()

for(i in 1:length(coast_name)) {
  retlistS <- rbind(retlistS, funRandCoast(coast_name[i]))
}


colnames(retlistS) <- c("GAM2","GAM2","RF2","RF2") 
summary(retlistS)
comparison <- as.data.frame(retlistS)
expl_var2 <- reshape2::melt(comparison[ , c(1, 3)])
MSE2 <- reshape2::melt(comparison[,c(2,4)])

# Check average output:
expl_var2 %>%
  group_by(variable) %>%
  summarise(mean = mean(value))

MSE2 %>%
  group_by(variable) %>%
  summarise(mean = mean(value))


expl_var3 <- rbind(expl_var, expl_var2)
MSE3 <- rbind(MSE, MSE2)


p <- ggplot() + 
  geom_boxplot(data = expl_var, aes(x = variable, y = value)) + 
  labs(y = "Variance", x = "") + 
  theme_base() +
  theme(plot.background = element_blank()) +
  ggtitle("Random bootsrap") +
  
  ggplot() + 
  geom_boxplot(data = MSE, aes(x = variable, y = value)) + 
  labs(y = "MSE", x = "") + 
  theme_base() +
  theme(plot.background = element_blank()) +
  ggtitle("Random bootsrap") +
  
  ggplot() + 
  geom_boxplot(data = expl_var2, aes(x = variable, y = value)) + 
  labs(y = "Variance", x = "") + 
  theme_base() + 
  theme(plot.background = element_blank()) +
  ggtitle("Regional bootsrap") +
  
  ggplot() + 
  geom_boxplot(data = MSE2, aes(x = variable, y = value)) + 
  labs(y = "MSE", x = "") + 
  theme_base() +
  theme(plot.background = element_blank()) +
  ggtitle("Regional bootsrap") 

p

ggsave("plots/statistics/a_size/power_analysis/a_size.png", p, height = 50 , width = 50, units = "mm", scale = 4)


# Run some t-tests:
# MSE
GAM <- subset(MSE3, variable == "GAM")[2]
RF <- subset(MSE3, variable == "RF")[2]
t.test(GAM, RF) # p-value = 4.608149e-41
t19 <- t.test(GAM, RF)$p.val

GAM2 <- subset(MSE3, variable == "GAM2")[2]
RF2 <- subset(MSE3, variable == "RF2")[2]
t.test(GAM2, RF2) # p-value = 
t20 <- t.test(GAM2, RF2)$p.val

# Variance
GAM <- subset(expl_var3, variable == "GAM")[2]
RF <- subset(expl_var3, variable == "RF")[2]
t.test(GAM, RF) # p-value = 3.127053e-121
t21 <- t.test(GAM, RF)$p.val

GAM2 <- subset(expl_var3, variable == "GAM2")[2]
RF2 <- subset(expl_var3, variable == "RF2")[2]
t.test(GAM2, RF2) # p-value = 0.07548468
t22 <- t.test(GAM2, RF2)$p.val



# Within GAM:
gam_random_var <- subset(expl_var3, variable == "GAM")[2]
gam_region_var <- subset(expl_var3, variable == "GAM2")[2]
t23 <- t.test(gam_random_var, gam_region_var)$p.val # p-value =  0.4903

gam_random_mse <- subset(MSE3, variable == "GAM")[2]
gam_region_mse <- subset(MSE3, variable == "GAM2")[2]
t24 <- t.test(gam_random_mse, gam_region_mse)$p.val # p-value = 0.5676


#--------------------------------------------------------------

# Lifespan:
funRand <- function(){
  # Test and training data 
  d1 <- stats_df4 %>% sample_n(max(0.66 * n(), 1))
  p1 <- stats_df4 %>%  anti_join(d1)
  
  # Fit and predict
  fit1 <- gam(mean_lifespan_weighted_wt_cpue_cephalopoda ~ 
                s(Temp_C, k = 3) +
                s(depth, k = 3) +
                s(z_prod, k = 3) + 
                # s(mean_wt_cpue_fish, k = 3) +
                s(coast, bs = "re"),
              data = d1,
              method = "REML",
              family = gaussian())
  
  f1 <- summary(fit1)$dev.expl
  
  fit2 <- randomForest(mean_lifespan_weighted_wt_cpue_cephalopoda ~ Temp_C + depth + z_prod + coast,
                       data = d1,
                       ntree = 1000,
                       importance = T,
                       mtry = 2,
                       na.action = na.omit)
  
  f2 <- mean(fit2$rsq)
  
  # Predict the remaining part and compare with observations
  Pred1 <- predict.gam(fit1, newdata = as.data.frame(p1), type = "response")
  Pred2 <- predict(fit2, newdata = as.data.frame(p1))
  
  # Mean squared errors
  R1 <- (p1[ , 11] - Pred1)^2
  R1 <- sum(R1$mean_lifespan_weighted_wt_cpue_cephalopoda)/(nrow(R1) - 1)
  R2 <- (p1[ , 11] - Pred2)^2 
  R2 <- sum(R2$mean_lifespan_weighted_wt_cpue_cephalopoda)/(nrow(R2) - 1)
  
  return(c(f1, R1, f2, R2))
}
funRand()

# Return results for runs 
retlistS <- t(unlist(replicate(n_replicates, funRand())))
colnames(retlistS) <- c("GAM","GAM","RF","RF") 
summary(retlistS)
comparison <- as.data.frame(retlistS)
expl_var <- reshape2::melt(comparison[ , c(1, 3)])
MSE <- reshape2::melt(comparison[,c(2,4)])

ggplot() + geom_boxplot(data = expl_var, aes(x = variable, y = value)) + labs(y = "Variance", x = "") + theme_bw() +
  ggplot() + geom_boxplot(data = MSE, aes(x = variable, y = value)) + labs(y = "MSE", x = "") + theme_bw() 


# Check average output:
expl_var %>%
  group_by(variable) %>%
  summarise(mean = mean(value))

MSE %>%
  group_by(variable) %>%
  summarise(mean = mean(value))


# Try to predict one region based on the other two regions:
funRandCoast <- function(coast_name){
  # Test and training data 
  d1 <- stats_df4 %>% filter(coast != coast_name)
  p1 <- stats_df4 %>%  filter(coast == coast_name)
  
  # Fit and predict
  fit1 <- gam(mean_lifespan_weighted_wt_cpue_cephalopoda ~ 
                s(Temp_C, k = 3) +
                s(depth, k = 3) +
                s(z_prod, k = 3) + 
                # s(mean_wt_cpue_fish, k = 3) +
                s(coast, bs = "re"),
              data = d1,
              method = "REML",
              family = gaussian())
  
  f1 <- summary(fit1)$dev.expl
  
  fit2 <- randomForest(mean_lifespan_weighted_wt_cpue_cephalopoda ~ Temp_C + depth + z_prod + coast,
                       data = d1,
                       ntree = 1000,
                       importance = T,
                       mtry = 2,
                       na.action = na.omit)
  
  f2 <- mean(fit2$rsq)
  
  # Predict the remaining part and compare with observations
  Pred1 <- predict.gam(fit1, newdata = as.data.frame(p1), type = "response")
  Pred2 <- predict(fit2, newdata = as.data.frame(p1))
  
  
  R1 <- (p1[ , 11] - Pred1)^2
  R1 <- sum(R1$mean_lifespan_weighted_wt_cpue_cephalopoda)/(nrow(R1) - 1)
  R2 <- (p1[ , 11] - Pred2)^2 
  R2 <- sum(R2$mean_lifespan_weighted_wt_cpue_cephalopoda)/(nrow(R2) - 1)
  
  return(c(f1, R1, f2, R2))
}


funRandCoast("Europe")

coast_name <- as.vector(unique(stats_df$coast))

retlistS <- c()

for(i in 1:length(coast_name)) {
  retlistS <- rbind(retlistS, funRandCoast(coast_name[i]))
}


colnames(retlistS) <- c("GAM2","GAM2","RF2","RF2") 
summary(retlistS)
comparison <- as.data.frame(retlistS)
expl_var2 <- reshape2::melt(comparison[ , c(1, 3)])
MSE2 <- reshape2::melt(comparison[,c(2,4)])

# Check average output:
expl_var2 %>%
  group_by(variable) %>%
  summarise(mean = mean(value))

MSE2 %>%
  group_by(variable) %>%
  summarise(mean = mean(value))


expl_var3 <- rbind(expl_var, expl_var2)
MSE3 <- rbind(MSE, MSE2)


p <- ggplot() + 
  geom_boxplot(data = expl_var, aes(x = variable, y = value)) + 
  labs(y = "Variance", x = "") + 
  theme_base() +
  theme(plot.background = element_blank()) +
  ggtitle("Random bootsrap") +
  
  ggplot() + 
  geom_boxplot(data = MSE, aes(x = variable, y = value)) + 
  labs(y = "MSE", x = "") + 
  theme_base() +
  theme(plot.background = element_blank()) +
  ggtitle("Random bootsrap") +
  
  ggplot() + 
  geom_boxplot(data = expl_var2, aes(x = variable, y = value)) + 
  labs(y = "Variance", x = "") + 
  theme_base() + 
  theme(plot.background = element_blank()) +
  ggtitle("Regional bootsrap") +
  
  ggplot() + 
  geom_boxplot(data = MSE2, aes(x = variable, y = value)) + 
  labs(y = "MSE", x = "") + 
  theme_base() +
  theme(plot.background = element_blank()) +
  ggtitle("Regional bootsrap") 

p

ggsave("plots/statistics/lifespan/power_analysis/lifespan.png", p, height = 50 , width = 50, units = "mm", scale = 4)


# Run some t-tests:
# MSE
GAM <- subset(MSE3, variable == "GAM")[2]
RF <- subset(MSE3, variable == "RF")[2]
t.test(GAM, RF) # p-value = 3.508845e-117
t25 <- t.test(GAM, RF)$p.val

GAM2 <- subset(MSE3, variable == "GAM2")[2]
RF2 <- subset(MSE3, variable == "RF2")[2]
t.test(GAM2, RF2) # p-value = 0.2723017
t26 <- t.test(GAM2, RF2)$p.val

# Variance
GAM <- subset(expl_var3, variable == "GAM")[2]
RF <- subset(expl_var3, variable == "RF")[2]
t.test(GAM, RF) # p-value = 9.065513e-160
t27 <- t.test(GAM, RF)$p.val

GAM2 <- subset(expl_var3, variable == "GAM2")[2]
RF2 <- subset(expl_var3, variable == "RF2")[2]
t.test(GAM2, RF2) # p-value = 0.3913116
t28 <- t.test(GAM2, RF2)$p.val


# Within GAM:
gam_random_var <- subset(expl_var3, variable == "GAM")[2]
gam_region_var <- subset(expl_var3, variable == "GAM2")[2]
t29 <- t.test(gam_random_var, gam_region_var)$p.val # p-value =  0.8559

gam_random_mse <- subset(MSE3, variable == "GAM")[2]
gam_region_mse <- subset(MSE3, variable == "GAM2")[2]
t30 <- t.test(gam_random_mse, gam_region_mse)$p.val # p-value = 0.03473


#--------------------------------------------------------------



#--------------------------------------------------------------

# A growth:
funRand <- function(){
  # Test and training data 
  d1 <- stats_df4 %>% sample_n(max(0.66 * n(), 1))
  p1 <- stats_df4 %>%  anti_join(d1)
  
  # Fit and predict
  fit1 <- gam(mean_A_growth_weighted_wt_cpue_cephalopoda ~ 
                s(Temp_C, k = 3) +
                s(depth, k = 3) +
                s(z_prod, k = 3) + 
                # s(mean_wt_cpue_fish, k = 3) +
                s(coast, bs = "re"),
              data = d1,
              method = "REML",
              family = gaussian())
  
  f1 <- summary(fit1)$dev.expl
  
  fit2 <- randomForest(mean_A_growth_weighted_wt_cpue_cephalopoda ~ Temp_C + depth + z_prod + coast,
                       data = d1,
                       ntree = 1000,
                       importance = T,
                       mtry = 2,
                       na.action = na.omit)
  
  f2 <- mean(fit2$rsq)
  
  # Predict the remaining part and compare with observations
  Pred1 <- predict.gam(fit1, newdata = as.data.frame(p1), type = "response")
  Pred2 <- predict(fit2, newdata = as.data.frame(p1))
  
  # Mean squared errors
  R1 <- (p1[ , 12] - Pred1)^2
  R1 <- sum(R1$mean_A_growth_weighted_wt_cpue_cephalopoda)/(nrow(R1) - 1)
  R2 <- (p1[ , 12] - Pred2)^2 
  R2 <- sum(R2$mean_A_growth_weighted_wt_cpue_cephalopoda)/(nrow(R2) - 1)
  
  return(c(f1, R1, f2, R2))
}
funRand()

# Return results for runs 
retlistS <- t(unlist(replicate(n_replicates, funRand())))
colnames(retlistS) <- c("GAM","GAM","RF","RF") 
summary(retlistS)
comparison <- as.data.frame(retlistS)
expl_var <- reshape2::melt(comparison[ , c(1, 3)])
MSE <- reshape2::melt(comparison[,c(2,4)])

p <- ggplot() + geom_boxplot(data = expl_var, aes(x = variable, y = value)) + labs(y = "Variance", x = "") + theme_bw() +
  ggplot() + geom_boxplot(data = MSE, aes(x = variable, y = value)) + labs(y = "MSE", x = "") + theme_bw() 


# Check average output:
expl_var %>%
  group_by(variable) %>%
  summarise(mean = mean(value))

MSE %>%
  group_by(variable) %>%
  summarise(mean = mean(value))


# Try to predict one region based on the other two regions:
funRandCoast <- function(coast_name){
  # Test and training data 
  d1 <- stats_df4 %>% filter(coast != coast_name)
  p1 <- stats_df4 %>%  filter(coast == coast_name)
  
  # Fit and predict
  fit1 <- gam(mean_A_growth_weighted_wt_cpue_cephalopoda ~ 
                s(Temp_C, k = 3) +
                s(depth, k = 3) +
                s(z_prod, k = 3) + 
                # s(mean_wt_cpue_fish, k = 3) +
                s(coast, bs = "re"),
              data = d1,
              method = "REML",
              family = gaussian())
  
  f1 <- summary(fit1)$dev.expl
  
  fit2 <- randomForest(mean_A_growth_weighted_wt_cpue_cephalopoda ~ Temp_C + depth + z_prod + coast,
                       data = d1,
                       ntree = 1000,
                       importance = T,
                       mtry = 2,
                       na.action = na.omit)
  
  f2 <- mean(fit2$rsq)
  
  # Predict the remaining part and compare with observations
  Pred1 <- predict.gam(fit1, newdata = as.data.frame(p1), type = "response")
  Pred2 <- predict(fit2, newdata = as.data.frame(p1))
  
  
  R1 <- (p1[ , 12] - Pred1)^2
  R1 <- sum(R1$mean_A_growth_weighted_wt_cpue_cephalopoda)/(nrow(R1) - 1)
  R2 <- (p1[ , 12] - Pred2)^2 
  R2 <- sum(R2$mean_A_growth_weighted_wt_cpue_cephalopoda)/(nrow(R2) - 1)
  
  return(c(f1, R1, f2, R2))
}


funRandCoast("Europe")

coast_name <- as.vector(unique(stats_df$coast))

retlistS <- c()

for(i in 1:length(coast_name)) {
  retlistS <- rbind(retlistS, funRandCoast(coast_name[i]))
}


colnames(retlistS) <- c("GAM2","GAM2","RF2","RF2") 
summary(retlistS)
comparison <- as.data.frame(retlistS)
expl_var2 <- reshape2::melt(comparison[ , c(1, 3)])
MSE2 <- reshape2::melt(comparison[,c(2,4)])

# Check average output:
expl_var2 %>%
  group_by(variable) %>%
  summarise(mean = mean(value))

MSE2 %>%
  group_by(variable) %>%
  summarise(mean = mean(value))


expl_var3 <- rbind(expl_var, expl_var2)
MSE3 <- rbind(MSE, MSE2)


p <- ggplot() + 
  geom_boxplot(data = expl_var, aes(x = variable, y = value)) + 
  labs(y = "Variance", x = "") + 
  theme_base() +
  theme(plot.background = element_blank()) +
  ggtitle("Random bootsrap") +
  
  ggplot() + 
  geom_boxplot(data = MSE, aes(x = variable, y = value)) + 
  labs(y = "MSE", x = "") + 
  theme_base() +
  theme(plot.background = element_blank()) +
  ggtitle("Random bootsrap") +
  
  ggplot() + 
  geom_boxplot(data = expl_var2, aes(x = variable, y = value)) + 
  labs(y = "Variance", x = "") + 
  theme_base() + 
  theme(plot.background = element_blank()) +
  ggtitle("Regional bootsrap") +
  
  ggplot() + 
  geom_boxplot(data = MSE2, aes(x = variable, y = value)) + 
  labs(y = "MSE", x = "") + 
  theme_base() +
  theme(plot.background = element_blank()) +
  ggtitle("Regional bootsrap") 

p

ggsave("plots/statistics/a_growth/power_analysis/a_growth.png", p, height = 50 , width = 50, units = "mm", scale = 4)

expl_var4 <- rbind(expl_var, expl_var2) %>%
  filter(variable %in% c("GAM", "GAM2"))

p <- ggplot() + 
  geom_boxplot(data = expl_var4, aes(x = variable, y = value, color = variable)) + 
  labs(y = "Variance", x = "") + 
  theme_base() +
  theme(plot.background = element_blank()) +
  ggtitle("R2 on growth parameter A") +
  scale_x_discrete(labels = c("Random","Regional")) +
  theme(legend.position = "none")

p

# Run some t-tests:
# MSE
GAM <- subset(MSE3, variable == "GAM")[2]
RF <- subset(MSE3, variable == "RF")[2]
t.test(GAM, RF) # p-value = 1.123668e-85
t31 <- t.test(GAM, RF)$p.val

GAM2 <- subset(MSE3, variable == "GAM2")[2]
RF2 <- subset(MSE3, variable == "RF2")[2]
t.test(GAM2, RF2) # p-value = 0.7888193
t32 <- t.test(GAM2, RF2)$p.val

# Variance
GAM <- subset(expl_var3, variable == "GAM")[2]
RF <- subset(expl_var3, variable == "RF")[2]
t.test(GAM, RF) # p-value = 1.297618e-140
t33 <- t.test(GAM, RF)$p.val

GAM2 <- subset(expl_var3, variable == "GAM2")[2]
RF2 <- subset(expl_var3, variable == "RF2")[2]
t.test(GAM2, RF2) # p-value = 0.2484645
t34 <- t.test(GAM2, RF2)$p.val

# Within GAM:
gam_random_var <- subset(expl_var3, variable == "GAM")[2]
gam_region_var <- subset(expl_var3, variable == "GAM2")[2]
t35 <- t.test(gam_random_var, gam_region_var)$p.val # p-value =  0.9668

gam_random_mse <- subset(MSE3, variable == "GAM")[2]
gam_region_mse <- subset(MSE3, variable == "GAM2")[2]
t36 <- t.test(gam_random_mse, gam_region_mse)$p.val # p-value = 0.1636

my_t_tests <- c(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10,
                t11, t12, t13, t14, t15, t16, t27, t28, t29, t20,
                t21, t22, t23, t24, t25, t26, t27, t28, t29, t30,
                t31, t32, t33, t34, t35, t36)




#                                 END OF SCRIPT
#################################################################################