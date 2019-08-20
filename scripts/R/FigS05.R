#-----------------------------#
#     de Jong et al 2019      #
# Code to reproduce Figure S5 #
#-----------------------------#

#### setup ####

# Load packages
library(tidyverse)
library(patchwork)

# Change ggplot2 default aesthetics
theme_set(theme_bw() + 
            theme(panel.grid = element_blank(), text = element_text(size = 10)))


#### read data ####

# Results from LIMIX multi-trait model
# ./scripts/python/limix_multitrait.py
gwas_multi <- read_csv("./data_processed/gwas_accessions/multi_trait_silique_early.csv") %>% 
  select(-X1) %>% 
  # retain only SNPS with MAF >= 5%
  filter(maf >= 0.05)


#### make plot ####

# Bonferroni threshold
thr <- gwas_multi %>% distinct(chrom, pos) %>% nrow()
thr <- 0.05/thr

# figure
pdf("./figures/figureS5.pdf", width = 7, height = 4)
gwas_multi %>% 
  select(chrom, pos, GxE = specific, G = common, maf) %>% 
  gather("component", "p", GxE:G) %>% 
  filter(p < 0.01 & maf > 0.05) %>% 
  ggplot(aes(pos/1e6, -log10(p))) +
  geom_point(size = 0.5) +
  facet_grid(component ~ chrom, scales = "free_x", space = "free_x") +
  geom_hline(yintercept = -log10(thr), linetype = "dashed", colour = "grey") +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(breaks = seq(0, 30, 10)) +
  labs(x = "Position (Mb)")
dev.off()


# One SNP gets above threshold (low MAF, noise?)
gwas_multi %>% 
  select(chrom, pos, GxE = specific, G = common, maf) %>% 
  gather("component", "p", GxE:G) %>% 
  filter(p < thr)

