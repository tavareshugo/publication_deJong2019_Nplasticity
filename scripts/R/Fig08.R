#-----------------------------#
#     de Jong et al 2019      #
# Code to reproduce Figure 08 #
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
gwas <- read_csv("./data_processed/gwas_accessions/compiled_gwas.csv")

# Set genome-wide threshold using bonferroni correction
# the number of SNPs from the 250K snp set is used
thr <- gwas %>% filter(trait == "bolt" & nitrate == "hn" & set == "snp250k") %>% nrow
thr <- 0.05/thr


#### set relaxed threshold ####

# set a threshold that includes two flowering QTL
thr2 <- 1e-5

#### make plot ####

# panel A
p1 <- gwas %>% 
  filter(trait == "bolt" & set == "snp250k" & p < 0.01 & stage == "silique") %>% 
  mutate(nitrate = ifelse(nitrate == "D", "Plasticity", toupper(nitrate))) %>% 
  mutate(sig = case_when(p < thr ~ "2", p < thr2 ~ "1", TRUE ~ "0")) %>% 
  ggplot(aes(bp/1e6, -log10(p))) +
  geom_point(aes(colour = sig, size = sig), show.legend = FALSE) +
  facet_grid(nitrate ~ Chr, scales = "free_x", space = "free_x") +
  geom_hline(yintercept = c(-log10(thr), -log10(thr2)), linetype = "dashed", colour = "grey") +
  theme(panel.grid = element_blank()) +
  scale_colour_manual(values = c("black", "orange", "brown"))  +
  scale_size_manual(values = c(0.5, 1, 1)) +
  scale_x_continuous(breaks = seq(0, 30, 10)) +
  labs(x = "Position (Mb)", tag = "A")

# panel B
p2 <- gwas %>% 
  filter(trait == "totalbr" & set == "snp250k" & p < 0.01 & stage == "silique") %>% 
  mutate(nitrate = ifelse(nitrate == "D", "Plasticity", toupper(nitrate))) %>% 
  mutate(sig = case_when(p < thr ~ "2", p < thr2 ~ "1", TRUE ~ "0")) %>% 
  ggplot(aes(bp/1e6, -log10(p))) +
  geom_point(aes(colour = sig, size = sig), show.legend = FALSE) +
  facet_grid(nitrate ~ Chr, scales = "free_x", space = "free_x") +
  geom_hline(yintercept = c(-log10(thr), -log10(thr2)), linetype = "dashed", colour = "grey") +
  theme(panel.grid = element_blank()) +
  scale_colour_manual(values = c("black", "orange", "brown"))  +
  scale_size_manual(values = c(0.5, 1, 1)) +
  scale_x_continuous(breaks = seq(0, 30, 10)) +
  labs(x = "Position (Mb)", tag = "B")

# Save figure
pdf("./figures/figure8.pdf", width = 7, height = 6)
p1 + p2 + plot_layout(ncol = 1)
dev.off()


# Using imputed SNPs (but showing the 250k threshold)
gwas %>% 
  filter(trait %in% c("bolt", "totalbr") & set == "imputed" & p < 0.01 & stage == "silique") %>% 
  mutate(nitrate = ifelse(nitrate == "D", "Plasticity", toupper(nitrate))) %>% 
  ggplot(aes(bp/1e6, -log10(p))) +
  geom_point(size = 0.5) +
  facet_grid(trait + nitrate ~ Chr, scales = "free_x", space = "free_x") +
  geom_hline(yintercept = -log10(thr), linetype = "dashed", colour = "grey") +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(breaks = seq(0, 30, 10)) +
  labs(x = "Position (Mb)")

# And senescence traits
gwas %>% 
  filter(trait %in% c("totalsil", "totalbr") & set == "snp250k" & p < 0.01 & stage == "senescence") %>% 
  mutate(nitrate = ifelse(nitrate == "D", "Plasticity", toupper(nitrate))) %>% 
  ggplot(aes(bp/1e6, -log10(p))) +
  geom_point(size = 0.5) +
  facet_grid(trait + nitrate ~ Chr, scales = "free_x", space = "free_x") +
  geom_hline(yintercept = -log10(thr), linetype = "dashed", colour = "grey") +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(breaks = seq(0, 30, 10)) +
  labs(x = "Position (Mb)")
