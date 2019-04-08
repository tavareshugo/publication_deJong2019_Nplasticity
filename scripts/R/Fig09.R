##
# Fig 10
# Get BLUP founder effects for a few of the QTL
##

#### setup ####

library(tidyverse)
library(patchwork)

# Change ggplot2 default aesthetics
theme_set(theme_bw() + 
            theme(panel.grid = element_blank(), 
                  text = element_text(size = 8),
                  plot.title = element_text(hjust = -0.07, size = 10, face = "bold", 
                                            margin = margin(t = -5, b = 1)),
                  plot.subtitle = element_text(size = 10, hjust = 0.5)))

# Source custom functions
source("./functions/qtl2_tidiers.R")


#### Read data ####

# Genotype objects
magic_gen <- read_rds("../../data/qtl_magic/magic_qtl2_cross2.rds")
magic_prob <- read_rds("../../data/qtl_magic/magic_qtl2_genoprob.rds")

# QTL scan results
qtl_results <- read_csv("../../data/qtl_magic/qtl2_scans_lm.csv")
qtl_thresh <- read_csv("../../data/qtl_magic/qtl2_scans_lm_perm.csv") %>% 
  filter(level == 0.05)

# Linear mixed model results
qtl_lmm <- read_csv("../../data/qtl_magic/qtl_scans_lmm.csv") %>% 
  select(-LOD_full) %>% 
  gather("test", "LOD", matches("LOD")) %>% 
  mutate(test = str_remove(test, "LOD_"))
qtl_lmm_thresh <- read_csv("../../data/qtl_magic/qtl2_scans_lmm_perm.csv") %>% 
  select(-LOD_full, -seed) %>% 
  summarise_all(funs(threshold = quantile(., probs = 0.95))) %>% 
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

p1 <- hgt_ln_2_blup %>% 
  mutate(coef = fct_reorder(coef, estimate)) %>% 
  ggplot(aes(coef, estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(ymax = hi, ymin = lo)) +
  labs(x = "Accession", y = "Allele effect", title = hgt_ln_2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_line(colour = "grey")) +
  scale_y_continuous(limits = c(-2, 2), breaks = -2:2)

p2 <- ft_hn_5_blup %>% 
  mutate(coef = fct_reorder(coef, estimate)) %>% 
  ggplot(aes(coef, estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(ymax = hi, ymin = lo)) +
  labs(x = "Accession", y = "Allele effect", title = ft_hn_5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_line(colour = "grey")) +
  scale_y_continuous(limits = c(-2, 2), breaks = -2:2)

p3 <- totalbr_d_2_blup %>% 
  mutate(coef = fct_reorder(coef, estimate)) %>% 
  ggplot(aes(coef, estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(ymax = hi, ymin = lo)) +
  labs(x = "Accession", y = "Allele effect", title = totalbr_d_2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_line(colour = "grey")) +
  scale_y_continuous(limits = c(-2, 2), breaks = -2:2)

pdf("./figures/figure10.pdf", height = 3, width = 7.5)
p1 + p2 + p3
dev.off()


