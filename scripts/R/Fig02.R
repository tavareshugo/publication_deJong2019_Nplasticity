#-----------------------------#
#     de Jong et al 2019      #
# Code to reproduce Figure 02 #
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
# We use three sets of data: 
# raw (individual)
# summarised (means per nitrate)
# summarised with plasticity traits (similar to above but in "wide" format)

# Read MAGIC line data
magic_ind <- read_rds("./data_processed/phenotypes/magic_individual.rds")
magic_plas <- read_rds("./data_processed/phenotypes/magic_plasticity.rds")
magic_sum <- read_rds("./data_processed/phenotypes/magic_summarised.rds")

# Read Accession data (note there is also senescence data for this set)
acc_ind <- read_rds("./data_processed/phenotypes/accessions_individual_silique.rds")
acc_plas <- read_rds("./data_processed/phenotypes/accessions_plasticity_silique.rds")
acc_sum <- read_rds("./data_processed/phenotypes/accessions_summarised_silique.rds")


#
# Choose lines to highlight ----
#
# Note, this is not intented as a formal inference test. It's more to highlight  
# that there is evidence against flat reaction norms for branching compared to flowering

# Rank-based test to highlight lines unlikely to have flat reaction norms
magic_highlight <- magic_ind %>% 
  group_by(id_leyser) %>% 
  summarise(bolt_pval = wilcox.test(bolt ~ nitrate, exact = FALSE) %>% broom::tidy() %>% pull(p.value),
            totalbr_pval = wilcox.test(totalbr ~ nitrate, exact = FALSE) %>% broom::tidy() %>% pull(p.value)) %>% 
  ungroup() %>% 
  mutate(bolt_padj = p.adjust(bolt_pval, method = "fdr"),
         totalbr_padj = p.adjust(totalbr_pval, method = "fdr")) %>% 
  filter(totalbr_padj < 0.05) %>% 
  pull(id_leyser)

acc_highlight <- acc_ind %>% 
  group_by(id_leyser) %>% 
  summarise(bolt_pval = wilcox.test(bolt ~ nitrate, exact = FALSE) %>% broom::tidy() %>% pull(p.value),
            totalbr_pval = wilcox.test(totalbr ~ nitrate, exact = FALSE) %>% broom::tidy() %>% pull(p.value)) %>% 
  ungroup() %>% 
  mutate(bolt_padj = p.adjust(bolt_pval, method = "fdr"),
         totalbr_padj = p.adjust(totalbr_pval, method = "fdr")) %>% 
  filter(totalbr_padj < 0.05) %>% 
  pull(id_leyser)


#
# plot ----
#
### MAGIC lines ###
# branches
p1 <- ggplot(magic_plas, aes(totalbr_mean_hn, totalbr_mean_ln)) +
  geom_point(size = 1) +
  geom_point(data = filter(magic_plas, id_leyser %in% magic_highlight), colour = "brown", size = 0.5) +
  annotate(geom = "text", x = 0, y = 10, 
           label = corLabel(magic_plas$totalbr_mean_hn, magic_plas$totalbr_mean_ln, TRUE),
           hjust = 0, vjust = 1, size = 2.5) +
  geom_abline(linetype = "dashed") +
  labs(x = "Total Branches (HN)", y = "Total Branches (LN)") +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  scale_y_continuous(breaks = seq(0, 10, 2)) +
  coord_fixed(xlim = c(0, 10), ylim = c(0, 10))

p1.1 <- magic_plas %>% 
  select(id_leyser, totalbr_rdpi) %>% 
  full_join(magic_sum, by = "id_leyser") %>% 
  ggplot(aes(nitrate, totalbr_mean)) +
  geom_line(aes(group = id_leyser, colour = totalbr_rdpi), alpha = 0.3, size = 1) +
  coord_cartesian(ylim = c(0, 10)) +
  scale_colour_viridis_c(limits = c(-1, 1)) +
  theme(legend.position = "none") +
  labs(x = "Nitrate", y = "Total Branches")

p1.2 <- magic_plas %>% 
  ggplot(aes(totalbr_rdpi)) + 
  geom_histogram(fill = "black", binwidth = 0.1) +
  #geom_histogram(data = filter(magic_plas, id_leyser %in% magic_highlight), fill = "brown", binwidth = 0.1) +
  geom_vline(xintercept = 0, colour = "grey", linetype = 2) +
  geom_rug(sides = "b", aes(colour = totalbr_rdpi), show.legend = FALSE) +
  scale_colour_viridis_c(limits = c(-1, 1)) +
  scale_x_continuous(limits = c(-1.2, 1.2)) +
  labs(x = "Relative Plasticity Index")


# Bolt
p2 <- ggplot(magic_plas, aes(2^boltLog_mean_hn, 2^boltLog_mean_ln)) +
  geom_point(size = 1) +
  annotate(geom = "text", x = 10, y = 108, 
           label = corLabel(magic_plas$boltLog_mean_hn, magic_plas$boltLog_mean_ln, TRUE),
           hjust = 0, vjust = 1, size = 2.5) +
  geom_abline(linetype = "dashed") +
  labs(x = "Days to Flowering (HN)", y = "Days to Flowering (LN)") +
  scale_x_continuous(breaks = seq(10, 110, 20), trans = "log2") +
  scale_y_continuous(breaks = seq(10, 110, 20), trans = "log2") +
  coord_fixed(xlim = c(10, 110), ylim = c(10, 110))

p2.1 <- magic_plas %>% 
  select(id_leyser, bolt_rdpi) %>% 
  full_join(magic_sum, by = "id_leyser") %>% 
  ggplot(aes(nitrate, 2^boltLog_mean)) +
  geom_line(aes(group = id_leyser, colour = bolt_rdpi), 
            alpha = 0.3, size = 1, show.legend = FALSE) +
  scale_y_continuous(breaks = seq(10, 110, 20), trans = "log2") +
  coord_cartesian(ylim = c(10, 110)) +
  scale_colour_viridis_c(limits = c(-1, 1)) +
  labs(x = "Nitrate", y = "Days to Flowering")

p2.2 <- magic_plas %>% 
  ggplot(aes(bolt_rdpi)) + 
  geom_histogram(fill = "black", binwidth = 0.1) +
  geom_vline(xintercept = 0, colour = "grey", linetype = 2) +
  geom_rug(sides = "b", aes(colour = bolt_rdpi), show.legend = FALSE) +
  scale_colour_viridis_c(limits = c(-1, 1)) +
  scale_x_continuous(limits = c(-1.2, 1.2)) +
  labs(x = "Relative Plasticity Index")


### Accessions ###
# branches
p3 <- ggplot(acc_plas, aes(totalbr_mean_hn, totalbr_mean_ln)) +
  geom_point(size = 1) +
  geom_point(data = filter(acc_plas, id_leyser %in% acc_highlight), colour = "brown", size = 0.5) +
  annotate(geom = "text", x = 0, y = 10, 
           label = corLabel(acc_plas$totalbr_mean_hn, acc_plas$totalbr_mean_ln, TRUE),
           hjust = 0, vjust = 1, size = 2.5) +
  geom_abline(linetype = "dashed") +
  labs(x = "Total Branches (HN)", y = "Total Branches (LN)") +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  scale_y_continuous(breaks = seq(0, 10, 2)) +
  coord_fixed(xlim = c(0, 10), ylim = c(0, 10))

p3.1 <- acc_plas %>% 
  select(id_leyser, totalbr_rdpi) %>% 
  full_join(acc_sum, by = "id_leyser") %>% 
  ggplot(aes(nitrate, totalbr_mean)) +
  geom_line(aes(group = id_leyser, colour = totalbr_rdpi), alpha = 0.3, size = 1) +
  coord_cartesian(ylim = c(0, 10)) +
  scale_colour_viridis_c(limits = c(-1, 1)) +
  theme(legend.position = "none") +
  labs(x = "Nitrate", y = "Total Branches")

p3.2 <- acc_plas %>% 
  ggplot(aes(totalbr_rdpi)) + 
  geom_histogram(fill = "black", binwidth = 0.1) +
  #geom_histogram(data = filter(magic_plas, id_leyser %in% magic_highlight), fill = "brown", binwidth = 0.1) +
  geom_vline(xintercept = 0, colour = "grey", linetype = 2) +
  geom_rug(sides = "b", aes(colour = totalbr_rdpi), show.legend = FALSE) +
  scale_colour_viridis_c(limits = c(-1, 1)) +
  scale_x_continuous(limits = c(-1.2, 1.2)) +
  labs(x = "Relative Plasticity Index")


# Bolt
p4 <- ggplot(acc_plas, aes(2^boltLog_mean_hn, 2^boltLog_mean_ln)) +
  geom_point(size = 1) +
  annotate(geom = "text", x = 10, y = 38, 
           label = corLabel(acc_plas$boltLog_mean_hn, acc_plas$boltLog_mean_ln, TRUE),
           hjust = 0, vjust = 1, size = 2.5) +
  geom_abline(linetype = "dashed") +
  labs(x = "Days to Flowering (HN)", y = "Days to Flowering (LN)") +
  scale_x_continuous(breaks = seq(10, 40, 5), trans = "log2") +
  scale_y_continuous(breaks = seq(10, 40, 5), trans = "log2") +
  coord_fixed(xlim = c(10, 40), ylim = c(10, 40))

p4.1 <- acc_plas %>% 
  select(id_leyser, bolt_rdpi) %>% 
  full_join(acc_sum, by = "id_leyser") %>% 
  ggplot(aes(nitrate, 2^boltLog_mean)) +
  geom_line(aes(group = id_leyser, colour = bolt_rdpi), 
            alpha = 0.3, size = 1, show.legend = FALSE) +
  scale_y_continuous(breaks = seq(10, 40, 5), trans = "log2") +
  coord_cartesian(ylim = c(10, 40)) +
  scale_colour_viridis_c(limits = c(-1, 1)) +
  labs(x = "Nitrate", y = "Days to Flowering")


p4.2 <- acc_plas %>% 
  ggplot(aes(bolt_rdpi)) + 
  geom_histogram(fill = "black", binwidth = 0.1) +
  geom_vline(xintercept = 0, colour = "grey", linetype = 2) +
  geom_rug(sides = "b", aes(colour = bolt_rdpi), show.legend = FALSE) +
  scale_colour_viridis_c(limits = c(-1, 1)) +
  scale_x_continuous(limits = c(-1.2, 1.2)) +
  labs(x = "Relative Plasticity Index")


# Put the plots together and save to pdf
#pdf("./figures/figure2.pdf", width = 7.5, height = 8.75)
p1 + labs(tag = "A") + p1.1 + p1.2 +
  p2 + labs(tag = "B") + p2.1 + p2.2 +
  p3 + labs(tag = "C") + p3.1 + p3.2 +
  p4 + labs(tag = "D") + p4.1 + p4.2 +
  plot_layout(ncol = 3, widths = c(0.4, 0.3, 0.3))
#dev.off()

