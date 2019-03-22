#-----------------------------#
#     de Jong et al 2019      #
# Code to reproduce Figure 03 #
#-----------------------------#


#### setup ####

# Load packages
library(tidyverse)
library(patchwork)

# Change ggplot2 default aesthetics
theme_set(theme_classic() + 
            theme(text = element_text(size = 8)))

# Custom function 
source("./scripts/R/functions/corrLabel.R")


#
# Read data ----
#
# Read MAGIC line data
magic_plas <- read_rds("./data_processed/phenotypes/magic_plasticity.rds")
magic_sum <- read_rds("./data_processed/phenotypes/magic_summarised.rds")


#
# Identify plasticity extremes ----
#
# Make subset of extremes using earlier flowering plants 
# (25 on each end of the distribution)
magic_extremes <- magic_plas %>% 
  filter(bolt_mean_ln < 25) %>% 
  arrange(abs(totalbr_mean_D)) %>% 
  slice(c(1:25, (n()-24):n())) %>% 
  mutate(plastic_category = rep(c("Least plastic", "Most plastic"), each = 25)) %>% 
  select(id_kover, plastic_category, totalbr_mean_D)

# Join with summary table
magic_extremes <- left_join(magic_extremes, magic_sum, by = "id_kover")



#
# plot ----
#
# total branches against plasticity. 
p1 <- ggplot(magic_plas, aes(x = totalbr_mean_D)) +
  geom_point(aes(y = totalbr_mean_hn, colour = "HN"), size = 0.8) +
  geom_point(aes(y = totalbr_mean_ln, colour = "LN"), size = 0.8) +
  annotate(geom = "text", x = 6, y = 4, size = 2.5, hjust = 0,
           label = corLabel(magic_plas$totalbr_mean_D, magic_plas$totalbr_mean_ln)) +
  annotate(geom = "text", x = 6, y = 6, size = 2.5, hjust = 0,
           label = corLabel(magic_plas$totalbr_mean_D, magic_plas$totalbr_mean_hn)) +
  scale_colour_manual(values = c("black", "grey69")) +
  labs(y = "Total Branches", x = "Branching Plasticity", colour = "Nitrate", tag = "A") +
  theme(legend.position = "top")


# Boxplot
p2 <- magic_extremes %>% 
  group_by(id_leyser) %>% 
  mutate(bolt = mean(bolt_mean)) %>% 
  ggplot(aes(nitrate, totalbr_mean, 
             group = id_kover, colour = bolt)) + 
  geom_line() +
  geom_point() +
  scale_colour_viridis_c(name = "Days to Flowering", breaks = seq(12, 28, 3)) +
  theme(legend.position = "top") +
  labs(x = "Nitrate", y = "Total Branches", tag = "B")


# Correlation with flowering time
# Also add the correlation for different sets of the data
p3 <- ggplot(magic_plas, aes(bolt_mean_ln, totalbr_mean_D)) +
  geom_point(size = 0.8) +
  geom_vline(xintercept = 25, linetype = "dotted") +
  annotate(geom = "text", x = 26, y = 8, size = 2.5, hjust = 0, vjust = 1,
           label = corLabel(magic_plas$bolt_mean_ln, magic_plas$totalbr_mean_D)) +
  annotate(geom = "text", x = 12, y = 8, size = 2.5, hjust = 0, vjust = 1,
           label = with(filter(magic_plas, bolt_mean_ln < 25),
                        corLabel(bolt_mean_ln, totalbr_mean_D))) +
  stat_smooth(se = FALSE, method = "loess") +
  labs(x = "Days to Flowering (LN)", y = "Branching Plasticity") +
  scale_x_continuous(breaks = seq(10, 100, 10), trans = "log2")


# Save plot
#pdf("./figures/figure3.pdf", width = 5.2, height = 4.5)
{p1 + labs(tag = "A") +
    p2 + labs(tag = "B") + plot_layout(widths = c(3/5, 2/5))} /
  p3 + labs(tag = "C")
#dev.off()
