
## Script name: All survey data statistics with size-specific catchability coefficient
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
## This script is part of the sensitivity analysis to evaluate how the catchability coefficient affects the statistical results.
## It uses biomass estimates of cephalopods calculated with a size-specific catchability coefficient. 
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
load(file = "data/ices_cephalopods_stats_data_size_q.Rdata")
stats_df1 <- stats_df %>%
  mutate(coast = "Europe")

load(file = "data/norway_cephalopods_stats_data_size_q.Rdata") 
stats_df2 <- stats_df %>%
  mutate(coast = "Europe")

load(file = "data/us_cephalopods_stats_data_size_q.Rdata")
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
            mean_prop = 100 * mean(cephalopoda_fish_proportion))


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
              s(coast, bs = "re"),
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
              s(coast, bs = "re"),
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
              s(coast, bs = "re"),
            data = stats_df4,
            method = "REML",
            family = gaussian())

# Time out:
print(Sys.time() - time0)

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

stats_df4 <- stats_df4 %>% 
  mutate(pred_a_size = fitted(gam3),
         pred_lifespan = fitted(gam4),
         pred_a_growth = fitted(gam5))

#-------------------------------------------------------------------------

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


#                                 END OF SCRIPT
#################################################################################