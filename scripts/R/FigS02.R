#-----------------------------#
#     de Jong et al 2019      #
# Code to reproduce Figure S1 #
#-----------------------------#

#
# setup ----
#
# Load packages
library(tidyverse)
library(patchwork)

# Change ggplot2 default aesthetics
theme_set(theme_classic() + 
            theme(text = element_text(size = 10)))

# Custom function 
source("./scripts/R/functions/corrLabel.R")


#
# read data ----
#

# Read seed yield data
seeds <- read_csv("./data/seed_yield/seed_yield.csv")

# Summarise data with means
seeds_sum <- seeds %>% 
  mutate(seedwgt = seedwgt*1000) %>% # convert weight to mg
  group_by(id_leyser, nitrate) %>% 
  summarise_at(vars(shootwgt, seedwgt, totalsil, totalbr), 
               list(mean = ~ mean(.), sd = ~ sd(.), se = ~ sd(.)/sqrt(n()))) %>% 
  ungroup() %>% 
  mutate(set = ifelse(str_detect(id_leyser, "accession"), "accession", "MAGIC"))

# Summarised accession data
acc_sum_sen <- read_rds("./data_processed/phenotypes/accessions_summarised_senescence.rds")

#
# plot ----
#
# Relationship between branch number and seed weight
p1 <- ggplot(seeds_sum, aes(totalbr_mean, seedwgt_mean, colour = nitrate)) +
  geom_point(aes(shape = set)) +
  geom_errorbar(aes(ymin = seedwgt_mean - 2*seedwgt_se, ymax = seedwgt_mean + 2*seedwgt_se), width = 0) +
  geom_errorbarh(aes(xmin = totalbr_mean - 2*totalbr_se, xmax = totalbr_mean + 2*totalbr_se), height = 0) +
  annotate(geom = "text", x = 0, y = 70, size = 3,
           label = corLabel(seeds_sum$totalbr_mean, seeds_sum$seedwgt_mean, TRUE),
           hjust = 0) +
  scale_colour_manual(values = c("black", "grey")) +
  labs(x = "Total Branches", y = "Seed weight (mg)", tag = "A") +
  theme(legend.position = "none") +
  scale_shape_manual(values = c(19, 4))

# Relationship between silique number and seed weight
p2 <- ggplot(seeds_sum, aes(totalsil_mean, seedwgt_mean, colour = nitrate)) +
  geom_point(aes(shape = set)) +
  geom_errorbar(aes(ymin = seedwgt_mean - 2*seedwgt_se, ymax = seedwgt_mean + 2*seedwgt_se)) +
  geom_errorbarh(aes(xmin = totalsil_mean - 2*totalsil_se, xmax = totalsil_mean + 2*totalsil_se)) +
  annotate(geom = "text", x = 25, y = 70, size = 3,
           label = corLabel(seeds_sum$totalsil_mean, seeds_sum$seedwgt_mean, TRUE),
           hjust = 0) +
  scale_colour_manual(values = c("black", "grey")) +
  labs(x = "Total Siliques", y = "Seed weight (mg)", colour = "Nitrate", tag = "B") +
  theme(legend.position = c(1, 0.1), legend.justification = c(1, 0.1)) +
  scale_shape_manual(values = c(19, 4), guide = FALSE)

# Relationship between branch number and silique number at senescence stage
p3 <- ggplot(acc_sum_sen, aes(totalbr_mean, totalsil_mean, colour = nitrate)) +
  geom_point() +
  annotate(geom = "text", x = 0, y = 200, size = 3, hjust = 0, 
           label = corLabel(acc_sum_sen$totalbr_mean, acc_sum_sen$totalsil_mean, TRUE)) +
  labs(x = "Total Branches", y = "Total Siliques", tag = "C") +
  scale_colour_manual(values = c("black", "grey")) +
  theme(legend.position = "none")

# Putting plots together
pdf("./figures/figureS2.pdf", width = 5.2, height = 5.2)
p1 + 
  p2 + 
  p3 + 
  plot_layout(ncol = 2)
dev.off()
