#-----------------------------#
#     de Jong et al 2019      #
# Code to reproduce Figure S6 #
#-----------------------------#

#### setup ####

# Load packages
library(tidyverse)
library(patchwork)

# Change ggplot2 default aesthetics
theme_set(theme_bw() + 
            theme(panel.grid = element_blank(), text = element_text(size = 10)))


#### read data ####

# produced with GCTA software with LOCO method for relatedness correction
# ./scripts/shell/04_gwas_gcta_250k.sh
gwas_gene <- read_csv("./data_processed/gwas_accessions/compiled_gwas_gene_level.csv")

# Set genome-wide threshold using bonferroni correction
thr_gene <- gwas_gene %>% filter(trait == "bolt" & nitrate == "hn" & set == "snp250k") %>% nrow
thr_gene <- 0.05/thr_gene

# Distribution of number of SNPs per window
gwas_gene %>% 
  filter(trait %in% c("bolt", "totalbr")) %>% 
  ggplot(aes(No.SNPs)) +
  geom_density() +
  facet_grid(trait + nitrate ~ set, scales = "free_x", space = "free_x")



#### make plots ####

p1 <- gwas_gene %>% 
  filter(trait == "bolt" & set == "snp250k" & Pvalue < 0.01) %>% 
  mutate(nitrate = ifelse(nitrate == "D", "Plasticity", toupper(nitrate))) %>% 
  ggplot(aes(Start/1e6, -log10(Pvalue))) +
  geom_point(size = 0.5) +
  facet_grid(nitrate ~ Chr, scales = "free_x", space = "free_x") +
  geom_hline(yintercept = -log10(thr_gene), linetype = "dashed", colour = "grey") +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(breaks = seq(0, 30, 10)) +
  scale_y_continuous(expand = c(0, 1)) +
  labs(x = "Position (Mb)", tag = "A")

p2 <- gwas_gene %>% 
  filter(trait == "totalbr" & set == "snp250k" & Pvalue < 0.01) %>% 
  mutate(nitrate = ifelse(nitrate == "D", "Plasticity", toupper(nitrate))) %>% 
  ggplot(aes(Start/1e6, -log10(Pvalue))) +
  geom_point(size = 0.5) +
  facet_grid(nitrate ~ Chr, scales = "free_x", space = "free_x") +
  geom_hline(yintercept = -log10(thr_gene), linetype = "dashed", colour = "grey") +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(breaks = seq(0, 30, 10)) +
  scale_y_continuous(expand = c(0, 1)) +
  labs(x = "Position (Mb)", tag = "B")

pdf("./figures/figureS6.pdf", width = 7, height = 6)
p1 + p2 + plot_layout(ncol = 1, heights = c(2/5, 3/5))
dev.off()


# Even with imputed set there's nothing striking
gwas_gene %>% 
  filter(trait %in% c("bolt", "totalbr") & set == "imputed" & Pvalue < 0.01) %>% 
  mutate(nitrate = ifelse(nitrate == "D", "Plasticity", toupper(nitrate))) %>% 
  ggplot(aes(Start/1e6, -log10(Pvalue))) +
  geom_point(size = 0.5) +
  facet_grid(trait + nitrate ~ Chr, scales = "free_x", space = "free_x") +
  geom_hline(yintercept = -log10(thr_gene), linetype = "dashed", colour = "grey") +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(breaks = seq(0, 30, 10)) +
  labs(x = "Position (Mb)")
