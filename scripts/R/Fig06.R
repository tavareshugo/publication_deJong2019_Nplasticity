#-----------------------------#
#     de Jong et al 2019      #
# Code to reproduce Figure 06 #
#-----------------------------#

#### setup ####

# Load packages
library(tidyverse)
library(patchwork)
library(lme4)
source("./scripts/R/functions/plotLmmDiag.R") # function to plot some diagnostics from lme4 object

# Change ggplot2 default aesthetics
theme_set(theme_bw() + 
            theme(panel.grid = element_blank(), text = element_text(size = 10)))



#### Read data ####

graft <- read_csv("./data/grafting/grafting_experiments.csv") %>% 
  mutate(experiment = as.character(experiment))

# Tidy data
graft <- graft %>% 
  # Remove plants that were too small
  filter(small != "small" | is.na(small)) %>% 
  # Add a unique ID for each individual
  # Make nitrate a factor
  mutate(id = paste(shoot_id, root_id, graft_type, sep = "\n"),
         nitrate = factor(nitrate, levels = c("HN", "LN")),
         type = ifelse(grepl("M11|M345", id), "MAGIC", "Accessions")) %>% 
  # Order ID by the shoot type (convenient for consistent plotting)
  mutate(id = reorder(as.factor(id), as.numeric(as.factor(paste0(shoot_type, graft_type)))),
         type = factor(type, levels = c("MAGIC", "Accessions")))



#### Inference about rooth/shoot on branches ####

# We use mixed model because design is unbalanced
graft %>% 
  count(experiment, id, nitrate) %>% 
  ggplot(aes(x = id, y = nitrate, label = n)) +
  geom_tile(fill = "white", colour = "black") +
  geom_text() +
  theme_minimal() +
  facet_grid(experiment ~ .)


# Fit model with random term accounting for graft replicates
# main effects we are interested in are the root and shoot types (and their interaction with N)
lmm1 <- lmer(totalbr ~ nitrate*shoot_type + nitrate*root_type + 
               experiment + (nitrate|id), 
             data = graft, REML = TRUE)

# More complex model - three-way interaction
lmm2 <- lmer(totalbr ~ nitrate*shoot_type*root_type + 
               experiment + (nitrate|id), 
             data = graft, REML = TRUE)

# Plot diagnostics
plotLmmDiag(lmm1)

# Information criteria suggest working with simpler model 
BIC(lmm2) - BIC(lmm1) # suggests loss of information
AIC(lmm2) - AIC(lmm1) # the two are quite close

# ANOVA from the model - results reported in the paper
car::Anova(lmm1, type = 2, test.statistic = "F")


#### Inference about rooth/shoot on flowering ####

# Fit similar model as above
bolt_lmm1 <- lmer(days_to_bolt ~ nitrate*shoot_type + nitrate*root_type + 
                    experiment + (nitrate|id), 
                  data = graft, REML = TRUE)

plotLmmDiag(bolt_lmm1)

# Run ANOVA
car::Anova(bolt_lmm1, test.statistic = "F", type = 2)



#### Make plot ####

# Sample size (to be added to legend)
graft %>% 
  count(experiment, id, nitrate) %>% 
  group_by(experiment) %>% 
  summarise(min(n), max(n), median(n))

p1 <- graft %>% 
  group_by(id, nitrate, experiment, type) %>% 
  summarise(mean = mean(totalbr, na.rm = TRUE),
            sd = sd(totalbr, na.rm = TRUE),
            n = sum(!is.na(totalbr)),
            se = 1.96*sd/sqrt(n)) %>% 
  mutate(id2 = paste(id, experiment)) %>% 
  ggplot(aes(interaction(nitrate, id), mean, ymin = mean - se, ymax = mean + se, colour = nitrate)) +
  geom_line(aes(group = id2, linetype = experiment), colour = "black") + 
  geom_point() +
  geom_errorbar(width = 0) +
  facet_grid( ~ type, scales = "free_x") +
  scale_linetype_discrete(guide = FALSE) +
  scale_colour_manual(guide = FALSE, values = c("black", "grey")) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(size = 3)) +
  labs(x = "", y = "Total Branches", colour = "Nitrate")

p2 <- graft %>% 
  group_by(id, nitrate, experiment, type) %>% 
  summarise(mean = mean(days_to_bolt, na.rm = TRUE),
            sd = sd(days_to_bolt, na.rm = TRUE),
            n = sum(!is.na(days_to_bolt)),
            se = 1.96*sd/sqrt(n)) %>% 
  mutate(id2 = paste(id, experiment)) %>% 
  ggplot(aes(interaction(nitrate, id), mean, ymin = mean - se, ymax = mean + se, colour = nitrate)) +
  geom_line(aes(group = id2, linetype = experiment), colour = "black") + 
  geom_point() +
  geom_errorbar(width = 0) +
  facet_grid( ~ type, scales = "free_x") +
  scale_linetype_discrete(guide = FALSE) +
  scale_colour_manual(values = c("black", "grey")) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(size = 3)) +
  labs(x = "", y = "Days to Flowering", colour = "Nitrate") +
  ylim(0, 24)

# Figure further edited in inkscape
pdf("./figures/figure6.pdf", width = 7, height = 3.5)
p2 + p1 + plot_layout(ncol = 1)
dev.off()
