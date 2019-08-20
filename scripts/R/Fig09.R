#-----------------------------#
#     de Jong et al 2019      #
# Code to reproduce Figure 09 #
#-----------------------------#

#### setup ####

library(tidyverse)
library(patchwork)
library(tidytext) # for reorder_within()

# Change ggplot2 default aesthetics
theme_set(theme_bw() + 
            theme(panel.grid = element_blank(), 
                  text = element_text(size = 10),
                  plot.title = element_text(hjust = -0.07, size = 10, face = "bold", 
                                            margin = margin(t = -5, b = 1)),
                  plot.subtitle = element_text(size = 10, hjust = 0.5)))

# Source custom functions
source("./scripts/R/functions/qtl2_tidiers.R")


#### Read data ####

# Genotype objects
magic_gen <- read_rds("./data_processed/qtl_magic/magic_qtl2_cross2.rds")
magic_prob <- read_rds("./data_processed/qtl_magic/magic_qtl2_genoprob.rds")

# QTL scan results
qtl_results <- read_csv("./data_processed/qtl_magic/qtl2_scans_lm.csv")
qtl_thresh <- read_csv("./data_processed/qtl_magic/qtl2_scans_lm_perm.csv") %>% 
  filter(level == 0.05)

# Linear mixed model results
qtl_lmm <- read_csv("./data_processed/qtl_magic/qtl2_scans_lmm.csv") %>% 
  select(-LOD_full) %>% 
  gather("test", "LOD", matches("LOD")) %>% 
  mutate(test = str_remove(test, "LOD_"))
qtl_lmm_thresh <- read_csv("./data_processed/qtl_magic/qtl2_scans_lmm_perm.csv") %>% 
  select(-LOD_full, -seed) %>% 
  summarise_all(list(threshold = ~ quantile(., probs = 0.95))) %>% 
  gather("test", "threshold") %>% 
  mutate(test = str_remove(test, "LOD_")) %>% 
  separate(test, c("test", "level"), sep = "_") 


#### Estimate BLUPs ####

#### Height LN Chr2
hgt_ln_2 <- qtl_results %>% 
  filter(trait == "height" & nitrate == "hn") %>% 
  arrange(desc(LOD)) %>% slice(1) %>% .$marker

hgt_ln_2_blup <- qtl2::scan1blup(list("1" = magic_prob[[2]][, , hgt_ln_2, drop = FALSE]), 
                                 scale(magic_gen$pheno[, "height_mean_hn"]),
                                 se = TRUE) %>% 
  tidy() %>% 
  filter(coef != "intercept") %>% 
  mutate(coef = str_remove(coef, magic_gen$alleles),
         lo = estimate - 1.96*SE,
         hi = estimate + 1.96*SE)


##### flowering HN Chr5
ft_hn_5 <- qtl_results %>% 
  filter(trait == "bolt" & nitrate == "hn" & chromosome == 5) %>% 
  arrange(desc(LOD)) %>% slice(1) %>% .$marker

ft_hn_5_blup <- qtl2::scan1blup(list("1" = magic_prob[[5]][, , ft_hn_5, drop = FALSE]), 
                                 scale(magic_gen$pheno[, "bolt_mean_hn"]),
                                 se = TRUE) %>% 
  tidy() %>% 
  filter(coef != "intercept") %>% 
  mutate(coef = str_remove(coef, magic_gen$alleles),
         lo = estimate - 1.96*SE,
         hi = estimate + 1.96*SE)


#### Branching plasticity Chr2
totalbr_d_2 <- qtl_results %>% 
  filter(trait == "totalbr" & nitrate == "D" & chromosome == 2) %>% 
  arrange(desc(LOD)) %>% slice(1) %>% .$marker

totalbr_d_2_blup <- qtl2::scan1blup(list("1" = magic_prob[[2]][, , totalbr_d_2, drop = FALSE]), 
                                scale(magic_gen$pheno[, "totalbr_mean_D"]),
                                se = TRUE) %>% 
  tidy() %>% 
  filter(coef != "intercept") %>% 
  mutate(coef = str_remove(coef, magic_gen$alleles),
         lo = estimate - 1.96*SE,
         hi = estimate + 1.96*SE)


#### Plot ####

# join all effects together
all_blups <- hgt_ln_2_blup %>% 
  bind_rows(ft_hn_5_blup) %>% 
  bind_rows(totalbr_d_2_blup) %>% 
  # rename markers more sensibly
  mutate(marker = case_when(marker == "MN2_11300378" ~ "HGT.LN.2",
                            marker == "MN5_4327715" ~ "FT.HN.5",
                            marker == "SOC1_461" ~ "SB.Pl.2")) %>% 
  # order them
  mutate(marker = factor(marker, levels = c("HGT.LN.2", "FT.HN.5", "SB.Pl.2")))

p1 <- all_blups %>% 
  mutate(coef = reorder_within(coef, estimate, marker)) %>% 
  ggplot(aes(coef, estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(ymax = hi, ymin = lo)) +
  labs(x = "Accession", y = "Standardised effect") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_line(colour = "grey"),
        strip.background = element_blank()) +
  scale_y_continuous(limits = c(-2, 2), breaks = -2:2) +
  scale_x_reordered() +
  facet_grid(~ marker, scales = "free_x")


pdf("./figures/figure9.pdf", height = 3, width = 7.5)
p1
dev.off()


