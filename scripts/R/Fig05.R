#-----------------------------#
#     de Jong et al 2019      #
# Code to reproduce Figure 06 #
#-----------------------------#

#### setup ####

# Load packages
library(tidyverse)
library(patchwork)

# Change ggplot2 default aesthetics
theme_set(theme_bw() + 
            theme(panel.grid = element_blank(), text = element_text(size = 8)))


#### Read data ####

# Read and tidy data
qpcr <- read_csv("./data/qpcr/nitrate_response_qpcr.csv") %>% 
  select(-line) %>%
  rename(line = line_name) %>% 
  mutate(line = factor(line, levels = c("Col-0", "MAGIC.11", "Hi-0", "Sha", "MAGIC.345", "Rsch-4", "Tsu-0")))


#### Normalise qPCR data ####
# Following Livak & Schmittgen 2001

# 1. Calculate average of technical replicates
qpcr_means <- qpcr %>% 
  group_by(sample_nr, gene, rep, line, line_type, nitrate) %>% 
  summarise(Cp_mean = mean(Cp, na.rm = TRUE), n = n()) %>% 
  ungroup()

# Check that all have 3 technical replicates
stopifnot(all(qpcr_means$n == 3))


# 2. Calculate the mean expression of the two reference genes
qpcr_means <- qpcr_means %>% 
  mutate(gene = ifelse(gene %in% c("APX3", "ubc9"), "reference", gene)) %>% 
  group_by(sample_nr, gene, rep, line, line_type, nitrate) %>% 
  summarise(Cp_mean = mean(Cp_mean)) %>% 
  ungroup()

# 3. Calculate the first delta-Cp: Cp_target - Cp_reference
# This normalises the gene expression in each sample by the reference genes
dcp <- qpcr_means %>% 
  group_by(sample_nr) %>% 
  # Add Cp of reference genes as a new column to calculate dCp
  mutate(Cp_ref = Cp_mean[gene == "reference"]) %>% 
  mutate(dCp = Cp_mean - Cp_ref) %>%
  ungroup() %>% 
  filter(gene != "reference") %>% 
  rename(Cp_target = Cp_mean)

# 4. Calculate the mean Cp in control samples per gene and genotype
dcp <- dcp %>% 
  group_by(line, gene) %>% 
  mutate(Cp_target_control = mean(Cp_target[nitrate == "control"]),
         Cp_ref_control = mean(Cp_ref[nitrate == "control"])) %>% 
  ungroup()

# 5. Calculate the delta-delta-Cp, i.e. the log2(fold-change)
dcp <- dcp %>% 
  mutate(ddCp = dCp - (Cp_target_control - Cp_ref_control),
         fold_change = 2^-ddCp)


#### Add branch information ####

line_ids <- read_csv("./data/qpcr/nitrate_response_qpcr.csv") %>% count(line, line_name)

# Read branching plasticity data
magics <- read_rds("./data_processed/phenotypes/magic_plasticity.rds") %>% 
  filter(id_leyser %in% c("11", "471")) %>% 
  select(id_kover, totalbr_mean_D) %>% 
  rename(line = id_kover)

accessions <- read_rds("./data_processed/phenotypes/accessions_plasticity_silique.rds") %>% 
  filter(id_leyser %in% c(36, 83, 86, 93, 95)) %>% 
  select(name, totalbr_mean_D) %>% 
  rename(line = name)

# Join with dcp table
dcp <- bind_rows(magics, accessions) %>% 
  full_join(dcp, by = "line")

# Turn variables to factor for plotting clarity
# Order line by branching plasticity
# Make GSR1 the first gene as it's the only one down-regulated
dcp <- dcp %>% 
  mutate(gene = fct_relevel(gene, "GSR1"),
         line = fct_reorder(line, totalbr_mean_D))



#### Make plot ####

p1 <- dcp %>% 
  filter(nitrate == "treat") %>% 
  mutate(plasticity = ifelse(totalbr_mean_D < 2, "< 2", "> 3")) %>% 
  ggplot(aes(line, log2(fold_change), colour = plasticity)) +
  geom_point(size = 0.8) +
  facet_grid( ~ gene, scales = "free_x") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_colour_manual(values = c("grey", "gray48")) +
  labs(x = "Line", y = expression(log[2]("KNO"[3]/"KCl")),
       colour = "Plasticity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = c(0, 1), legend.justification = c(-0.1, 1.1),
        legend.background = element_rect(color = "black", size = 0.5))

pdf("./figures/figure6.pdf", width = 7, height = 2.5)
p1
dev.off()


# For figure legend, add plasticity from experiments
dcp %>% 
  distinct(line, totalbr_mean_D) %>% 
  mutate(temp = paste0(line, " = ", round(totalbr_mean_D, 1))) %>% 
  mutate(temp = fct_reorder(temp, totalbr_mean_D)) %>% 
  pull(temp) %>% 
  levels() %>% 
  paste(collapse = ", ")


#### Effect of genotype ####

# Run a Kruskall-Wallis test per gene
# it's quite low-powered (two reps only) so unsurprising that p-values are not that low
dcp %>% 
  filter(nitrate == "treat") %>% 
  group_by(gene) %>% 
  nest() %>% 
  mutate(lm_fit = map(data, ~ kruskal.test(log2(fold_change) ~ line, data = .))) %>% 
  mutate(lm_fit = map(lm_fit, ~ broom::glance(.))) %>% 
  unnest(lm_fit)

# Check correlation between branch number and expression
dcp %>% 
  filter(nitrate == "treat") %>% 
  mutate(fold_change = log2(fold_change)) %>% 
  group_by(gene) %>% 
  nest() %>% 
  mutate(cor = map(data, ~ cor.test(.$fold_change, .$totalbr_mean_D, method = "spearman"))) %>% 
  mutate(cor = map(cor, ~ broom::glance(.))) %>% 
  unnest(cor) %>% 
  mutate(padj = p.adjust(p.value))
