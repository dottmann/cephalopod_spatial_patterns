
## Script name: SAU stats
##
## Authors: Daniel Ottmann 
## Email: daniel.ottmann.riera@gmail.com
##
## Date created: March 2023
## Last update:  March 2023
##
## ---------------------------
##
## Readme:
##
## This script runs statistics on the global fisheries data fromSea Around Us
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
# library(sjstats)
library(gratia)
library(car)
library(randomForest)
library(MixRF)
library(pdp)
library(scales)

# Clear environment:
rm(list = ls())

# Create "exclude" function:
"%ni%" <- Negate("%in%")

# Load data:
load(file = "data/sau_data_stats.Rdata")

stats_df <- sau_data_stats

##########################################################################

# Edit data:
stats_df <- stats_df %>%
  mutate(cephalopod_biomass = case_when(is.na(cephalopod_biomass) ~ 0,
                                        T ~ cephalopod_biomass),
         cephalopod_fish_proportion = case_when(is.na(cephalopod_fish_proportion) ~ 0,
                                                T ~ cephalopod_fish_proportion)) %>%
  as.data.frame()



#-------------------------------------------------------------------------------------------

# Check variance inflation factors:
vif <- lm(cephalopod_fish_proportion ~ depth + z_prod + Temp_C + fish_biomass,
          data = stats_df)

vif(vif)

# All VIF values are lower than 1.2 

# Check correlation between variables:
cor_table <- stats_df %>% 
  dplyr::select(depth, z_prod, Temp_C, fish_biomass)

cor(cor_table)

# Correlation =< 0.3


#-------------------------------------------------------------------------------------------


############
# Run stats!

#
# GAM models -----------------------------------------------------------
#

###########################
# Cephalopod proportion:
###########################
gam0 <- gam(cephalopod_fish_proportion ~ 
              s(Temp_C, k = 3) +
              s(depth, k = 3) +
              s(z_prod, k = 3),
            data = stats_df,
            method = "REML",
            family = betar())

summary(gam0)
gam.check(gam0)

# Check smooths:
p <- draw(gam0) &
  theme_classic()

p
# ggsave("plots/GAM_ceph_proportion_smooths.png", p, height = 40 , width = 40, units = "mm", scale = 4)

# Save panels for joint plot:
a1 <-  draw(gam0,          
            ci_alpha = 0.2,       
            ci_col = "orange",         
            smooth_col = "orange")[[2]] +
  xlab("Temperature (ºC)") + 
  ylim(c(-2.1, .7)) +
  ylab(expression("P"["Ceph"])) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank()
        # ,
        # axis.title.x = element_blank()
  ) +
  geom_hline(yintercept = 0, color = "gray")

a2 <-  draw(gam0,      
            ci_alpha = 0.2,         
            ci_col = "orange",          
            smooth_col = "orange")[[1]] +
  xlab("Depth (m)") + 
  ylim(c(-2.1, .7)) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")

a3 <-  draw(gam0,       
            ci_alpha = 0.2,    
            ci_col = "orange",    
            smooth_col = "orange")[[3]] +
  xlab(expression(paste("Zooplankton productivity (g ", m^-2, " ", yr^-1, ")"))) + 
  ylim(c(-2.1, .7)) +
  theme_classic()  +
  theme(plot.background = element_blank(),
        plot.title = element_blank(),
        axis.title.y = element_blank())  +
  geom_hline(yintercept = 0, color = "gray")


# Combine all smooths in a plot:
p <- (a1  / a2 / a3) +
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 10)) & 
  theme(plot.background = element_blank())

p
# ggsave("plots/all_smooths.png", p, height = 60 , width = 28, units = "mm", scale = 3)

rm(list = c("a1", "a2", "a3"))


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
  ggplot(aes(x = Temp_C, y = resid)) +
  theme_classic() +
  geom_point(size = 3, shape = 1) + 
  geom_hline(yintercept = 0, color = "gray")

# Residuals vs z_prod:
stats_df %>% 
  mutate(resid = residuals(gam0, type = "response")) %>% 
  ggplot(aes(x = z_prod, y = resid)) +
  theme_classic() +
  geom_point(size = 3, shape = 1) + 
  geom_hline(yintercept = 0, color = "gray")

# Residuals vs depth:
stats_df %>% 
  mutate(resid = residuals(gam0, type = "response")) %>% 
  ggplot(aes(x = depth, y = resid)) +
  theme_classic() +
  geom_point(size = 3, shape = 1) + 
  geom_hline(yintercept = 0, color = "gray")


#-----------------------------------------------------------------

# Add predictions to the data frame:
stats_df2 <- stats_df %>% 
  mutate(pred_proportion = predict(gam0, type = "response"),
         trunc_vals = case_when(cephalopod_fish_proportion > .2 ~ .2,
                                T ~ cephalopod_fish_proportion),
         model_residuals = cephalopod_fish_proportion - pred_proportion)

# Save predictions:
save(stats_df2, file = "data/sau_predictions.Rdata")


#-------------------------------------------------------------------------

# Plot predictions:

####################
# Cephalopod proportion:

# Make predictions of range of a variable keeping all other variables constant:
# Loop to obtain partial effect predictions for all the variables:
variables <- c("Temp_C", "depth", "z_prod")
partial_plot_proportion <- data.frame()

test <- stats_df2 %>%
  dplyr::select(Temp_C, depth, z_prod) %>%
  as.data.frame()

for(i in 1:length(variables)){
  var1 <- seq(min(test[,which(colnames(test) == variables[i])]),
              max(test[,which(colnames(test) == variables[i])]), length = nrow(test))
  
  vars <- variables[-i]
  
  var2 <- rep(mean(test[,which(colnames(test) == vars[1])]),length(var1)) 
  var3 <- rep(mean(test[,which(colnames(test) == vars[2])]),length(var1))
  
  
  newdat <- data.frame(var1, var2, var3) 
  
  colnames(newdat)[1:3] <- c(variables[i], vars)
  
  my_predictions <- predict.gam(gam0, newdata = newdat, type = "response", se.fit = TRUE) # predict with new data excluding random effect
  my_fit <- my_predictions$fit
  my_se <- my_predictions$se.fit
  
  partial <- data.frame(my_predictions, variable = var1) 
  partial <- partial %>% mutate(variable = var1,
                                variable_id = variables[i])
  
  partial_plot_proportion = rbind(partial_plot_proportion, partial)
  
}


# Plot as pointranges:
stats_df3 <- stats_df2 %>%
  group_by(lme) %>%
  summarise(n = n(),
            sd_Temp_C = sd(Temp_C, na.rm = T),
            sd_cephalopod_fish_proportion = sd(cephalopod_fish_proportion, na.rm = T),
            se_Temp_C = sd_Temp_C / sqrt(n),
            se_cephalopod_fish_proportion = 100 * sd_cephalopod_fish_proportion / sqrt(n),
            cephalopod_fish_proportion = 100 * mean(cephalopod_fish_proportion, na.rm = T),
            Temp_C = mean(Temp_C, na.rm = T))

#---------------------------------------------------------------------------------

# Temperature:
p <- ggplot() +
  
  geom_point(data = stats_df2, aes(x = Temp_C, y = 100 * trunc_vals), alpha = .5, size = 1) +
  geom_point(data = subset(stats_df2, trunc_vals == .2), aes(x = Temp_C, y = 100 * trunc_vals), alpha = .5, size = 1, color = "red") +
  
  geom_smooth(data = stats_df2, aes(x = Temp_C, y = 100 * pred_proportion), color = "black", se = T, alpha = .3) +         # Predicted regression
  
  geom_line(data = subset(partial_plot_proportion, variable_id == "Temp_C"), aes(x = variable, y = 100 * fit), color = "orange", size = 1) + 
  geom_ribbon(data = subset(partial_plot_proportion, variable_id == "Temp_C"), aes(y = 100 * fit, x = variable, ymin = 100 * (fit - se.fit), ymax = 100 * (fit + se.fit)),
              fill = "orange", alpha = .1) +

  ylim(c(0, 20)) +
  theme_base() +
  ylab("% ceph. biom.") +
  xlab("Temperature (ºC)") + 
  theme(plot.background = element_blank(),
        legend.title = element_blank()) 

p
a1 <- p  


# Depth:
p <- ggplot() +
  
  geom_point(data = stats_df2, aes(x = depth, y = 100 * trunc_vals), alpha = .5, size = 1) +
  geom_point(data = subset(stats_df2, trunc_vals == .2), aes(x = depth, y = 100 * trunc_vals), alpha = .5, size = 1, color = "red") +
  
  geom_smooth(data = stats_df2, aes(x = depth, y = 100 * pred_proportion), color = "black", se = T, alpha = .3) +         # Predicted regression
  
  geom_line(data = subset(partial_plot_proportion, variable_id == "depth"), aes(x = variable, y = 100 * fit), color = "orange", size = 1) + 
  geom_ribbon(data = subset(partial_plot_proportion, variable_id == "depth"), aes(y = 100 * fit, x = variable, ymin = 100 * (fit - se.fit), ymax = 100 * (fit + se.fit)),
              fill = "orange", alpha = .1) +
  
  ylim(c(0, 20)) +
  theme_base() +
  ylab("") +
  xlab("Depth (m)") + 
  theme(plot.background = element_blank(),
        legend.title = element_blank()) 

p
a2 <- p 


# Zooplankton productivity:
p <- ggplot() +
  
  geom_point(data = stats_df2, aes(x = z_prod, y = 100 * trunc_vals), alpha = 0.5, size = 1) +
  geom_point(data = subset(stats_df2, trunc_vals == .2), aes(x = z_prod, y = 100 * trunc_vals), alpha = .5, size = 1, color = "red") +
  
  geom_smooth(data = stats_df2, aes(x = z_prod, y = 100 * pred_proportion), color = "black", se = T, alpha = .3) +         # Predicted regression
  
  geom_line(data = subset(partial_plot_proportion, variable_id == "z_prod"), aes(x = variable, y = 100 * fit), color = "orange", size = 1) + 
  geom_ribbon(data = subset(partial_plot_proportion, variable_id == "z_prod"), aes(y = 100 * fit, x = variable, ymin = 100 * (fit - se.fit), ymax = 100 * (fit + se.fit)),
              fill = "orange", alpha = .1) +

  ylim(c(0, 20)) +
  theme_base() +
  ylab("") +
  xlab(expression(paste("Zoo. productivity (g ", m^-2, " ", yr^-1, ")"))) + 
  theme(plot.background = element_blank(),
        legend.title = element_blank()) 

p
a3 <- p 

# Combine all smooths in a plot:
p <- (a1 | a2 | a3) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 15)) &
  theme(plot.background = element_blank())

p
ggsave("plots/regression_plots.png", p, height = 20 , width = 66, units = "mm", scale = 4)

rm(list = c("a1", "a2", "a3"))

#-------------------------------------------------------------------------
#
# Random Forest -----------------------------------------------------------
#

# Cephalopod proportion:
rf0 <- randomForest(cephalopod_fish_proportion ~ Temp_C + depth + z_prod, #+ mean_wt_cpue_fish
                    data = stats_df,
                    ntree = 1000,
                    importance = T,
                    mtry = 2,
                    na.action = na.omit)

rf0

# Partial plots:
vars <- c("Temp_C", "depth", "z_prod")

for(i in 1:length(vars)) {
  p <- rf0 %>%  
    pdp::partial(pred.var = vars[i]) %>%
    autoplot(smooth = TRUE) + 
    theme_base() + 
    theme(axis.text = element_text(size = 30),
          axis.title.x = element_text(size = 30),
          axis.title.y = element_text(size = 30))
  
  ggsave(paste0("plots/ceph_proportion_", vars[i],".png"), p, height = 40 , width = 40, units = "mm", scale = 4)
}

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
  fit1 <- gam(cephalopod_fish_proportion ~ 
                s(Temp_C, k = 3) +
                s(depth, k = 3) +
                s(z_prod, k = 3), 
                # s(mean_wt_cpue_fish, k = 3) +
              data = d1,
              method = "REML",
              family = betar())
  
  f1 <- summary(fit1)$dev.expl
  
  fit2 <- randomForest(cephalopod_fish_proportion ~ Temp_C + depth + z_prod,
                       data = d1,
                       ntree = 1000,
                       importance = T,
                       mtry = 2,
                       na.action = na.omit)
  
  f2 <- mean(fit2$rsq)
  
  # Predict the remaining part and compare with observations
  Pred1 <- predict.gam(fit1, newdata = as.data.frame(p1), type = "response")
  Pred2 <- c(predict(fit2, newdata = as.data.frame(p1)))
  
  # Mean squared errors
  R1 <- (p1[ , 4] - Pred1)^2
  R1 <- sum(R1)/(nrow(Pred1) - 1)
  R2 <- c((p1[ , 4] - Pred2)^2)
  R2 <- sum(R2)/(length(Pred2) - 1)
  
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

p <-  ggplot() + 
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
  ggtitle("Random bootsrap") 

p

expl_var %>%
  group_by(variable) %>%
  summarise(mean = mean(value))

MSE %>%
  group_by(variable) %>%
  summarise(mean = mean(value))


ggsave("plots/ceph_proportion.png", p, height = 50 , width = 50, units = "mm", scale = 4)

# Check signifficance:
# Run some t-tests:
# MSE
GAM <- subset(MSE, variable == "GAM")[2]
RF <- subset(MSE, variable == "RF")[2]
t.test(GAM, RF) 

# Variance
GAM <- subset(expl_var, variable == "GAM")[2]
RF <- subset(expl_var, variable == "RF")[2]
t.test(GAM, RF) 
