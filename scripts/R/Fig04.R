#-----------------------------#
#     de Jong et al 2019      #
# Code to reproduce Figure 04 #
#-----------------------------#


#### setup ####

# Load packages
library(tidyverse)
library(patchwork)

# Change ggplot2 default aesthetics
theme_set(theme_classic() + 
            theme(text = element_text(size = 10)))

# Custom function 
source("./scripts/R/functions/corrLabel.R")


#
# Read data ----
#
# Read Accession data - silique dataset
acc_plas <- read_rds("./data_processed/phenotypes/accessions_plasticity_silique.rds")
acc_sum <- read_rds("./data_processed/phenotypes/accessions_summarised_silique.rds")

# Read senescence data
acc_plas_sen <- read_rds("./data_processed/phenotypes/accessions_plasticity_senescence.rds")


#
# Identify plasticity extremes ----
#
# Make subset of extremes (25 on each end of the distribution)
acc_extremes <- acc_plas %>% 
  filter(bolt_mean_ln < 25) %>% 
  arrange(totalbr_mean_D) %>% 
  slice(c(1:25, (n()-24):n())) %>% 
  mutate(plastic_category = rep(c("Least plastic", "Most plastic"), each = 25)) %>% 
  select(id_leyser, plastic_category, totalbr_mean_D)

# Join with summary table
acc_extremes <- left_join(acc_extremes, acc_sum, by = "id_leyser")


#
# plot ----
#
# total branches against plasticity - add Pearson's correlation
p1 <- ggplot(acc_plas, aes(x = totalbr_mean_D)) +
  geom_point(aes(y = totalbr_mean_hn, colour = "HN"), size = 0.8) +
  geom_point(aes(y = totalbr_mean_ln, colour = "LN"), size = 0.8) +
  annotate(geom = "text", x = 4.2, y = 2.4, size = 3, hjust = 0,
           label = corLabel(acc_plas$totalbr_mean_D, acc_plas$totalbr_mean_ln)) +
  annotate(geom = "text", x = 4.2, y = 4, size = 3, hjust = 0,
           label = corLabel(acc_plas$totalbr_mean_D, acc_plas$totalbr_mean_hn)) +
  scale_colour_manual(values = c("black", "grey69")) +
  labs(y = "Total Branches", x = "Branching Plasticity", colour = "Nitrate", tag = "A") +
  theme(legend.position = "top") +
  guides(colour = guide_legend(title.position="top", direction = "vertical"))


# Boxplot
p2 <- acc_extremes %>% 
  group_by(id_leyser) %>% 
  mutate(bolt = mean(bolt_mean)) %>% 
  ggplot(aes(nitrate, totalbr_mean, 
             group = id_leyser, colour = bolt)) + 
  geom_line() +
  geom_point() +
  scale_colour_viridis_c(name = "Days to Flowering", breaks = seq(12, 28, 3)) +
  #theme(legend.position = c(1, 1), legend.justification = c(1, 1)) +
  theme(legend.position = "top") +
  labs(x = "Nitrate", y = "Total Branches", tag = "B") +
  guides(colour = guide_colourbar(title.position="top"))


# Plot of total siliques against plasticity. 
p3 <- ggplot(acc_plas_sen, aes(x = totalbr_mean_D)) +
  geom_point(aes(y = totalsil_mean_hn, colour = "HN"), size = 0.8) +
  geom_point(aes(y = totalsil_mean_ln, colour = "LN"), size = 0.8) +
  annotate(geom = "text", x = 4, y = 70, size = 3, hjust = 0, vjust = 1,
           label = corLabel(acc_plas_sen$totalbr_mean_D, acc_plas_sen$totalsil_mean_ln)) +
  annotate(geom = "text", x = 4, y = 220, size = 3, hjust = 0, vjust = 1,
           label = corLabel(acc_plas_sen$totalbr_mean_D, acc_plas_sen$totalsil_mean_hn)) +
  scale_colour_manual(values = c("black", "grey69")) +
  labs(y = "Total Siliques", x = "Branching Plasticity", colour = "Nitrate") +
  theme(legend.position = "none")

# Correlation with flowering time
p4 <- ggplot(acc_plas, aes(bolt_mean_ln, totalbr_mean_D)) +
  geom_point(size = 0.8) +
  geom_vline(xintercept = 25, linetype = "dotted") +
  annotate(geom = "text", x = 26, y = 8, size = 3, hjust = 0, vjust = 1,
           label = corLabel(acc_plas$bolt_mean_ln, acc_plas$totalbr_mean_D)) +
  annotate(geom = "text", x = 15, y = 8, size = 3, hjust = 0, vjust = 1,
           label = with(filter(acc_plas, bolt_mean_ln < 25),
                        corLabel(bolt_mean_ln, totalbr_mean_D))) +
  stat_smooth(se = FALSE, method = "loess") +
  labs(x = "Days to Flowering (LN)", y = "Branching Plasticity") +
  scale_x_continuous(breaks = seq(10, 100, 10), trans = "log2")

# Assemble figure
pdf("./figures/figure4.pdf", width = 5.2, height = 4.5)
{p1 + labs(tag = "A") +
    p2 + labs(tag = "B")} /
    {p3 + labs(tag = "C") +
        p4 + labs(tag = "D")}
dev.off()
