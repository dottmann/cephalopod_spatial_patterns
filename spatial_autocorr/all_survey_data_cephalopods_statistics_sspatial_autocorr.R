
## Script name: All survey data statistics with catchability coefficient = 0.3 and spatial dimension
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
## To evaluate ho sensitive is the analysis to spatial variables, in this script we repeat the analysis using spatial variables
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

# Check variance inflation factors of the variables:
vif <- lm(mean_wt_cpue_cephalopoda ~ depth + z_prod + Temp_C + mean_wt_cpue_fish + lat + lon,
          data = stats_df)

vif(vif)

# VIF for lat is too high due to correlation with temperature and z_prod.
# if we remove it:
vif <- lm(mean_wt_cpue_cephalopoda ~ depth + z_prod + Temp_C + mean_wt_cpue_fish + lon,
          data = stats_df)

vif(vif)

# All VIF values are lower than 2

# Check correlation of variables:
cor_table <- stats_df %>%
  filter(coast == "Europe") %>%
  dplyr::select(depth, z_prod, Temp_C, mean_wt_cpue_fish, lon, lat)

cor(cor_table)

# Same issue with latitude:


# Including lon in North American waters is pointless because the data are essentially 2 vertical lines,
# and in European waters, lat and lon are correlated:
cor(subset(stats_df, coast == "Europe")$lon, subset(stats_df, coast == "Europe")$lat)

# So, including spatial structures would supress the effect of explanatory variables.

# If we remove lat and lon: 
# Check correlation between variables:

vif <- lm(mean_wt_cpue_cephalopoda ~ depth + z_prod + Temp_C + mean_wt_cpue_fish,
          data = stats_df)

vif(vif)

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


#-------------------------------------------------------------------------------------------


############
# Run stats!


#
# GAM models -----------------------------------------------------------
#

###########################
# Cephalopod proportion:
###########################

# Time in:
time0 <- Sys.time()

gam0 <- gamm(cephalopoda_fish_proportion ~ 
              s(Temp_C, k = 3) +
              s(depth, k = 3) +
              s(z_prod, k = 3) + 
              s(coast, bs = "re"),
            correlation = corExp(form = ~ lat + lon),
            data = stats_df,
            method = "REML") # spatial correlation structure does not work with betar() family

# Time out:
print(Sys.time() - time0)

# save(gam0, file = "models/gam0.Rdata")

summary(gam0$lme)

# Check smooths:
p <- draw(gam0) &
  theme_classic()

p
# ggsave("plots/statistics/GAMM_ceph_proportion_smooths.png", p, height = 40 , width = 40, units = "mm", scale = 4)

# Save panels for joint plot:
c1 <-  draw(gam0,          
            ci_alpha = 0.2,       
            ci_col = "orange",         
            smooth_col = "orange")[[3]] +
  ylim(c(-.025, .035)) +
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
  ylim(c(-.025, .035)) +
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
  ylim(c(-.025, .035)) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")


#######################
# Cephalopod biomass:
#######################

# Time in:
time0 <- Sys.time()

gam1 <- gamm(mean_wt_cpue_cephalopoda ~ 
              s(Temp_C, k = 3) +
              s(depth, k = 3) +
              s(z_prod, k = 3) + 
              # s(log10(mean_wt_cpue_fish), k = 3) +
              s(coast, bs = "re"),
            correlation = corExp(form = ~ lat + lon),
            data = stats_df,
            method = "REML",
            family = nb())

# Time out:
print(Sys.time() - time0)

# save(gam1, file = "models/gam1.Rdata")

summary(gam1$gam)

# Check smooths:
p <- draw(gam1) &
  theme_classic()

p
# ggsave("plots/statistics/GAMM_ceph_biomass_smooths.png", p, height = 40 , width = 40, units = "mm", scale = 4)

# Save panels for joint plot:
a1 <-  draw(gam1,        
            ci_alpha = 0.2,     
            ci_col = "orange",      
            smooth_col = "orange")[[3]] +
  ylim(c(-4.5, 3))  +
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
  ylim(c(-4.5, 3))  +
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
  ylim(c(-4.5, 3)) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")






#######################
# Fish biomass:
#######################

# Time in:
time0 <- Sys.time()

gam2 <- gamm(mean_wt_cpue_fish ~ 
              s(Temp_C, k = 3) +
              s(depth, k = 3) +
              s(z_prod, k = 3) + 
              # s(mean_wt_cpue_cephalopoda, k = 3) +
              s(coast, bs = "re"),
            correlation = corExp(form = ~ lat + lon),
            data = stats_df,
            method = "REML",
            family = nb())

# Time out:
print(Sys.time() - time0)

# save(gam2, file = "models/gam2.Rdata")

summary(gam2$gam)

# Check smooths:
p <- draw(gam2) &
  theme_classic()

p
# ggsave("plots/statistics/GAMM_fish_biomass_smooths.png", p, height = 40 , width = 40, units = "mm", scale = 4)

# Save panels for joint plot:
b1 <-  draw(gam2, 
            ci_alpha = 0.2,  
            ci_col = "orange",         
            smooth_col = "orange")[[3]] +
  ylim(c(-2, 1.2)) +
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
  ylim(c(-2, 1.2)) +
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
  ylim(c(-2, 1.2)) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")




#######################
# Asymptotic size:
#######################

# Time in:
time0 <- Sys.time()

gam3 <- gamm(log10(mean_a_size_weighted_wt_cpue_cephalopoda) ~ 
              s(Temp_C, k = 3) +
              s(depth, k = 3) +
              s(z_prod, k = 3) + 
              # s(log10(mean_wt_cpue_fish), k = 3) +
              s(coast, bs = "re"),
            correlation = corExp(form = ~ lat + lon),
            data = stats_df4,
            method = "REML",
            family = gaussian())

# Time out:
print(Sys.time() - time0)

# save(gam3, file = "models/gam3.Rdata")

summary(gam3$gam)

# Check smooths:
p <- draw(gam3) &
  theme_classic()

p
# ggsave("plots/statistics/GAMM_ceph_a_size_smooths.png", p, height = 40 , width = 40, units = "mm", scale = 4)

# Save panels for joint plot:
d1 <-  draw(gam3,
            ci_alpha = 0.2,
            ci_col = "orange",
            smooth_col = "orange")[[3]] +
  ylim(c(-.5, .4))  +
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
  ylim(c(-.5, .4))  +
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
  ylim(c(-.5, .4)) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")





#######################
# Lifespan:
#######################

# Time in:
time0 <- Sys.time()

gam4 <- gamm(mean_lifespan_weighted_wt_cpue_cephalopoda ~ 
              s(Temp_C, k = 3) +
              s(depth, k = 3) +
              s(z_prod, k = 3) + 
              # s(mean_wt_cpue_fish, k = 3) +
              s(coast, bs = "re"),
            correlation = corExp(form = ~ lat + lon),
            data = stats_df4,
            method = "REML",
            family = gaussian())

# Time out:
print(Sys.time() - time0)

# save(gam4, file = "models/gam4.Rdata")

summary(gam4$gam)

# Check smooths:
p <- draw(gam4) &
  theme_classic()

p
# ggsave("plots/statistics/GAMM_ceph_lifespan_smooths.png", p, height = 40 , width = 40, units = "mm", scale = 4)

# Save panels for joint plot:
e1 <-  draw(gam4,
            ci_alpha = 0.2,
            ci_col = "orange",
            smooth_col = "orange")[[3]] +
  ylim(c(-3, 2.5)) +
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
  ylim(c(-3, 2.5)) +
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
  ylim(c(-3, 2.5)) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")




#######################
# A growth:
#######################

# Time in:
time0 <- Sys.time()

gam5 <- gamm(mean_A_growth_weighted_wt_cpue_cephalopoda ~ 
              s(Temp_C, k = 3) +
              s(depth, k = 3) +
              s(z_prod, k = 3) + 
              # s(log10(mean_wt_cpue_fish), k = 3) +
              s(coast, bs = "re"),
            correlation = corExp(form = ~ lat + lon),
            data = stats_df4,
            method = "REML",
            family = gaussian())

# Time out:
print(Sys.time() - time0)

# save(gam5, file = "models/gam5.Rdata")

summary(gam5)

# Check smooths:
p <- draw(gam5) &
  theme_classic()

p
# ggsave("plots/statistics/GAMM_a_growth_gam0_smooths.png", p, height = 40 , width = 40, units = "mm", scale = 4)

# Save panels for joint plot:
f1 <-  draw(gam5,
            ci_alpha = 0.2,
            ci_col = "orange",
            smooth_col = "orange")[[3]] +
  ylim(c(-7, 9)) +
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
  ylim(c(-7, 9)) +
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
  ylim(c(-7, 9)) +
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




#-------------------------------------------------------------------------

# Add predictions to the data frame:
stats_df3 <- stats_df %>% 
  mutate(pred_proportion = predict(gam0$gam, type = "response"), # Same output as fitted()
         pred_ceph_biomass = fitted(gam1$gam),
         pred_fish_biomass = fitted(gam2$gam))

stats_df4 <- stats_df4 %>% 
  mutate(pred_a_size = fitted(gam3$gam) ,
         pred_lifespan = fitted(gam4$gam) ,
         pred_a_growth = fitted(gam5$gam) )


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
  
  my_predictions <- predict(gam0$gam, newdata = newdat, exclude = "s(coast)", type = "response", se.fit = TRUE) # predict with new data excluding random effect
  my_fit <- my_predictions$fitmy_fit <- my_predictions$fit
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
  ylim(c(-1, NA)) +
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
  ylim(c(-1, NA)) +
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
  ylim(c(-1, NA)) +
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
  
  my_predictions <- predict(gam1$gam, newdata = newdat, exclude = "s(coast)", type = "response", se.fit = TRUE) # predict with new data excluding random effect
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
  
  my_predictions <- predict(gam2$gam, newdata = newdat, exclude = "s(coast)", type = "response", se.fit = TRUE) # predict with new data excluding random effect
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
  
  my_predictions <- predict(gam3$gam, newdata = newdat, exclude = "s(coast)", type = "response", se.fit = TRUE) # predict with new data excluding random effect
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
  
  my_predictions <- predict(gam4$gam, newdata = newdat, exclude = "s(coast)", type = "response", se.fit = TRUE) # predict with new data excluding random effect
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
  
  my_predictions <- predict(gam5$gam, newdata = newdat, exclude = "s(coast)", type = "response", se.fit = TRUE) # predict with new data excluding random effect
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


# Combine all partial effect data frames:

partial_plots_df <- rbind(partial_plot_proportion, 
                          partial_plot_ceph_bio, 
                          partial_plot_fish_bio, 
                          partial_plot_a_size, 
                          partial_plot_lifespan, 
                          partial_plot_a_growth)

# save(partial_plots_df, file = "data/partial_plots_df.Rdata")



#                                 END OF SCRIPT
#################################################################################