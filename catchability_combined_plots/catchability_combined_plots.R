
## Script name: Combined smooth plots with different catchability coefficient
##
## Authors: Daniel Ottmann 
## Email: daniel.ottmann.riera@gmail.com
##
## Date created: June 2023
## Last update:  July 2023
##
## ---------------------------
##
## Readme:
##
## This script combines the smooth plots of the GAMMs produced with different catchability coefficients to facilitate their comparison.
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


#-----------------------------------------------------------------------

# Catchability coefficient = 0.3

#-----------------------------------------------------------------------

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

# Save panels for joint plot:
q03_c1 <-  draw(gam0,          
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

q03_c2 <-  draw(gam0,      
            ci_alpha = 0.2,         
            ci_col = "orange",          
            smooth_col = "orange")[[2]] +
  ylim(c(-1.5, 1.5)) +
  ylab(expression("P"["Ceph"])) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")

q03_c3 <-  draw(gam0,       
            ci_alpha = 0.2,    
            ci_col = "orange",    
            smooth_col = "orange")[[4]] +
  ylim(c(-2.1, 1.5)) +
  ylab(expression("P"["Ceph"])) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank())  +
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

# Save panels for joint plot:
q03_a1 <-  draw(gam1,        
            ci_alpha = 0.2,     
            ci_col = "orange",      
            smooth_col = "orange")[[3]] +
  ylim(c(-6.5, 3))  +
  ylab(expression(paste("B"["Ceph"], " (ln(kg ", km^-2,"))"))) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        # plot.title = element_blank(),
        axis.title.x = element_blank())  +
  geom_hline(yintercept = 0, color = "gray") +
  ggtitle("Same catchability")

q03_a2 <-  draw(gam1,      
            ci_alpha = 0.2, 
            ci_col = "orange",      
            smooth_col = "orange")[[2]] +
  ylim(c(-6.5, 3))  +
  ylab(expression(paste("B"["Ceph"], " (ln(kg ", km^-2,"))"))) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        # plot.title = element_blank(),
        axis.title.x = element_blank())  +
  geom_hline(yintercept = 0, color = "gray") +
  ggtitle("Same catchability") 

q03_a3 <-  draw(gam1,     
            ci_alpha = 0.2,   
            ci_col = "orange",   
            smooth_col = "orange")[[4]] +
  ylim(c(-6.5, 3)) +
  ylab(expression(paste("B"["Ceph"], " (ln(kg ", km^-2,"))"))) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        # plot.title = element_blank(),
        axis.title.x = element_blank())  +
  geom_hline(yintercept = 0, color = "gray") +
  ggtitle("Same catchability")



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


# Save panels for joint plot:
q03_b1 <-  draw(gam2, 
            ci_alpha = 0.2,  
            ci_col = "orange",         
            smooth_col = "orange")[[3]] +
  ylim(c(-1.5, 2)) +
  ylab(expression(paste("B"["Fish"], " (ln(kg ", km^-2,"))"))) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")

q03_b2 <-  draw(gam2,
            ci_alpha = 0.2,
            ci_col = "orange", 
            smooth_col = "orange")[[2]]  +
  ylim(c(-1.5, 2)) +
  ylab(expression(paste("B"["Fish"], " (ln(kg ", km^-2,"))"))) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")

q03_b3 <-  draw(gam2,
            ci_alpha = 0.2,
            ci_col = "orange",
            smooth_col = "orange")[[4]] +
  ylim(c(-5, 2)) +
  ylab(expression(paste("B"["Fish"], " (ln(kg ", km^-2,"))"))) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank())  +
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


# Save panels for joint plot:
q03_d1 <-  draw(gam3,
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

q03_d2 <-  draw(gam3,
            ci_alpha = 0.2,
            ci_col = "orange",
            smooth_col = "orange")[[2]] +
  ylim(c(-1.2, .5))  +
  ylim(c(-1.2, .5))  +
  ylab(expression(paste("W"[infinity], " log(g)"))) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")

q03_d3 <-  draw(gam3,
            ci_alpha = 0.2,
            ci_col = "orange",
            smooth_col = "orange")[[4]] +
  ylim(c(-1.2, .5)) +
  ylim(c(-1.2, .5))  +
  ylab(expression(paste("W"[infinity], " log(g)"))) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank())  +
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


# Save panels for joint plot:
q03_e1 <-  draw(gam4,
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

q03_e2 <-  draw(gam4,
            ci_alpha = 0.2,
            ci_col = "orange",
            smooth_col = "orange")[[2]] +
  ylim(c(-10, 8.5)) +
  ylab(expression(paste("R"[0]," (months)"))) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")

q03_e3 <-  draw(gam4,
            ci_alpha = 0.2,
            ci_col = "orange",
            smooth_col = "orange")[[4]] +
  ylim(c(-10, 8.5)) +
  ylab(expression(paste("R"[0]," (months)"))) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank())  +
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


# Save panels for joint plot:
q03_f1 <-  draw(gam5,
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

q03_f2 <-  draw(gam5,
            ci_alpha = 0.2,
            ci_col = "orange",
            smooth_col = "orange")[[2]] +
  ylim(c(-7, 10.5)) +
  ylab(expression(paste("A (", g^(1/3), yr^-1,")"))) +
  xlab("Depth (m)") +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")

q03_f3 <-  draw(gam5,
            ci_alpha = 0.2,
            ci_col = "orange",
            smooth_col = "orange")[[4]] +
  ylim(c(-7, 10.5)) +
  ylab(expression(paste("A (", g^(1/3), yr^-1,")"))) +
  xlab(expression(paste("Z. Product (g ", m^-2, " ", yr^-1, ")"))) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")




#-----------------------------------------------------------

# No catchability coefficient

#-----------------------------------------------------------

# Load data:
load(file = "data/ices_cephalopods_stats_data_sensitivity_noq.Rdata")
stats_df1 <- stats_df %>%
  mutate(coast = "Europe")

load(file = "data/norway_cephalopods_stats_data_sensitivity_noq.Rdata") 
stats_df2 <- stats_df %>%
  mutate(coast = "Europe")

load(file = "data/us_cephalopods_stats_data_sensitivity_noq.Rdata")
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


# Save panels for joint plot:
nq_c1 <-  draw(gam0,          
            ci_alpha = 0.2,       
            ci_col = "orange",         
            smooth_col = "orange")[[3]] +
  ylim(c(-1.5, 1.5)) +
  ylab(expression("P"["Ceph"])) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_hline(yintercept = 0, color = "gray")

nq_c2 <-  draw(gam0,      
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

nq_c3 <-  draw(gam0,       
            ci_alpha = 0.2,    
            ci_col = "orange",    
            smooth_col = "orange")[[4]] +
  ylim(c(-2.1, 1.5)) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())  +
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


# Save panels for joint plot:
nq_a1 <-  draw(gam1,        
            ci_alpha = 0.2,     
            ci_col = "orange",      
            smooth_col = "orange")[[3]] +
  ylim(c(-6.5, 3))  +
  ylab(expression(paste("B"["Ceph"], " (ln(kg ", km^-2,"))"))) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        # plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray") +
  ggtitle("No catchability")

nq_a2 <-  draw(gam1,      
            ci_alpha = 0.2, 
            ci_col = "orange",      
            smooth_col = "orange")[[2]] +
  ylim(c(-6.5, 3))  +
  theme_classic()  +
  theme(plot.background = element_blank(),
        # plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray") +
  ggtitle("No catchability")

nq_a3 <-  draw(gam1,     
            ci_alpha = 0.2,   
            ci_col = "orange",   
            smooth_col = "orange")[[4]] +
  ylim(c(-6.5, 3)) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        # plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray") +
  ggtitle("No catchability")



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


# Save panels for joint plot:
nq_b1 <-  draw(gam2, 
            ci_alpha = 0.2,  
            ci_col = "orange",         
            smooth_col = "orange")[[3]] +
  ylim(c(-1.5, 2)) +
  ylab(expression(paste("B"["Fish"], " (ln(kg ", km^-2,"))"))) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")

nq_b2 <-  draw(gam2,
            ci_alpha = 0.2,
            ci_col = "orange", 
            smooth_col = "orange")[[2]]  +
  ylim(c(-1.5, 2)) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")

nq_b3 <-  draw(gam2,
            ci_alpha = 0.2,
            ci_col = "orange",
            smooth_col = "orange")[[4]] +
  ylim(c(-5, 2)) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())  +
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


# Save panels for joint plot:
nq_d1 <-  draw(gam3,
            ci_alpha = 0.2,
            ci_col = "orange",
            smooth_col = "orange")[[3]] +
  ylim(c(-1.2, .5))  +
  ylab(expression(paste("W"[infinity], " log(g)"))) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")

nq_d2 <-  draw(gam3,
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

nq_d3 <-  draw(gam3,
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

# Save panels for joint plot:
nq_e1 <-  draw(gam4,
            ci_alpha = 0.2,
            ci_col = "orange",
            smooth_col = "orange")[[3]] +
  ylim(c(-10, 8.5)) +
  ylab(expression(paste("R"[0]," (months)"))) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")

nq_e2 <-  draw(gam4,
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

nq_e3 <-  draw(gam4,
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


# Save panels for joint plot:
nq_f1 <-  draw(gam5,
            ci_alpha = 0.2,
            ci_col = "orange",
            smooth_col = "orange")[[3]] +
  ylim(c(-7, 10.5)) +
  ylab(expression(paste("A (", g^(1/3), yr^-1,")"))) +
  xlab("Temperature (ºC)") +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")

nq_f2 <-  draw(gam5,
            ci_alpha = 0.2,
            ci_col = "orange",
            smooth_col = "orange")[[2]] +
  ylim(c(-7, 10.5)) +
  ylab("") +
  xlab("Depth (m)") +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")

nq_f3 <-  draw(gam5,
            ci_alpha = 0.2,
            ci_col = "orange",
            smooth_col = "orange")[[4]] +
  ylim(c(-7, 10.5)) +
  ylab("") +
  xlab(expression(paste("Z. Product (g ", m^-2, " ", yr^-1, ")"))) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")



#--------------------------------------------------------------------------

# Size-speciffic catchability

#--------------------------------------------------------------------------


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


# Save panels for joint plot:
sizeq_c1 <-  draw(gam0,          
            ci_alpha = 0.2,       
            ci_col = "orange",         
            smooth_col = "orange")[[3]] +
  ylim(c(-1.5, 1.5)) +
  ylab(expression("P"["Ceph"])) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_hline(yintercept = 0, color = "gray")

sizeq_c2 <-  draw(gam0,      
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

sizeq_c3 <-  draw(gam0,       
            ci_alpha = 0.2,    
            ci_col = "orange",    
            smooth_col = "orange")[[4]] +
  ylim(c(-2.15, 1.5)) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())  +
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


# Save panels for joint plot:
sizeq_a1 <-  draw(gam1,        
            ci_alpha = 0.2,     
            ci_col = "orange",      
            smooth_col = "orange")[[3]] +
  ylim(c(-6.5, 3))  +
  ylab(expression(paste("B"["Ceph"], " (ln(kg ", km^-2,"))"))) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        # plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray") +
  ggtitle("Size-based catchability")

sizeq_a2 <-  draw(gam1,      
            ci_alpha = 0.2, 
            ci_col = "orange",      
            smooth_col = "orange")[[2]] +
  ylim(c(-6.5, 3))  +
  theme_classic()  +
  theme(plot.background = element_blank(),
        # plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray") +
  ggtitle("Size-based catchability")

sizeq_a3 <-  draw(gam1,     
            ci_alpha = 0.2,   
            ci_col = "orange",   
            smooth_col = "orange")[[4]] +
  ylim(c(-6.5, 3)) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        # plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray") +
  ggtitle("Size-based catchability")



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


# Save panels for joint plot:
sizeq_b1 <-  draw(gam2, 
            ci_alpha = 0.2,  
            ci_col = "orange",         
            smooth_col = "orange")[[3]] +
  ylim(c(-1.5, 2)) +
  ylab(expression(paste("B"["Fish"], " (ln(kg ", km^-2,"))"))) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")

sizeq_b2 <-  draw(gam2,
            ci_alpha = 0.2,
            ci_col = "orange", 
            smooth_col = "orange")[[2]]  +
  ylim(c(-1.5, 2)) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")

sizeq_b3 <-  draw(gam2,
            ci_alpha = 0.2,
            ci_col = "orange",
            smooth_col = "orange")[[4]] +
  ylim(c(-5, 2)) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())  +
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

# Save panels for joint plot:
sizeq_d1 <-  draw(gam3,
            ci_alpha = 0.2,
            ci_col = "orange",
            smooth_col = "orange")[[3]] +
  ylim(c(-1.2, .5))  +
  ylab(expression(paste("W"[infinity], " log(g)"))) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")

sizeq_d2 <-  draw(gam3,
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

sizeq_d3 <-  draw(gam3,
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

# Save panels for joint plot:
sizeq_e1 <-  draw(gam4,
            ci_alpha = 0.2,
            ci_col = "orange",
            smooth_col = "orange")[[3]] +
  ylim(c(-10, 8.5)) +
  ylab(expression(paste("R"[0]," (months)"))) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")

sizeq_e2 <-  draw(gam4,
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

sizeq_e3 <-  draw(gam4,
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

# Save panels for joint plot:
sizeq_f1 <-  draw(gam5,
            ci_alpha = 0.2,
            ci_col = "orange",
            smooth_col = "orange")[[3]] +
  ylim(c(-7, 10.5)) +
  ylab(expression(paste("A (", g^(1/3), yr^-1,")"))) +
  xlab("Temperature (ºC)") +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")

sizeq_f2 <-  draw(gam5,
            ci_alpha = 0.2,
            ci_col = "orange",
            smooth_col = "orange")[[2]] +
  ylim(c(-7, 10.5)) +
  ylab("") +
  xlab("Depth (m)") +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")

sizeq_f3 <-  draw(gam5,
            ci_alpha = 0.2,
            ci_col = "orange",
            smooth_col = "orange")[[4]] +
  ylim(c(-7, 10.5)) +
  ylab("") +
  xlab(expression(paste("Z. Product (g ", m^-2, " ", yr^-1, ")"))) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")



#------------------------------

p <- (q03_a1 | nq_a1 | sizeq_a1) /
  (q03_b1 | nq_b1 | sizeq_b1) / 
  (q03_c1 | nq_c1 | sizeq_c1) /
  (q03_d1 | nq_d1 | sizeq_d1) /
  (q03_e1 | nq_e1 | sizeq_e1) /
  (q03_f1 | nq_f1 | sizeq_f1) +
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 10)) & 
  theme(plot.background = element_blank())

p
ggsave("plots/smooths_temp.png", p, height = 80 , width = 60, units = "mm", scale = 3)



p <- (q03_a2 | nq_a2 | sizeq_a2) /
  (q03_b2 | nq_b2 | sizeq_b2) / 
  (q03_c2 | nq_c2 | sizeq_c2) /
  (q03_d2 | nq_d2 | sizeq_d2) /
  (q03_e2 | nq_e2 | sizeq_e2) /
  (q03_f2 | nq_f2 | sizeq_f2) +
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 10)) & 
  theme(plot.background = element_blank())

p
# ggsave("sensitivity_q/plots/smooths_depth.png", p, height = 80 , width = 60, units = "mm", scale = 3)




p <- (q03_a3 | nq_a3 | sizeq_a3) /
  (q03_b3 | nq_b3 | sizeq_b3) / 
  (q03_c3 | nq_c3 | sizeq_c3) /
  (q03_d3 | nq_d3 | sizeq_d3) /
  (q03_e3 | nq_e3 | sizeq_e3) /
  (q03_f3 | nq_f3 | sizeq_f3) +
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 10)) & 
  theme(plot.background = element_blank())

p
# ggsave("sensitivity_q/plots/smooths_z_prod.png", p, height = 80 , width = 60, units = "mm", scale = 3)


#                    END OF SCRIPT
############################################################