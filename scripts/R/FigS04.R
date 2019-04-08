#-----------------------------#
#     de Jong et al 2019      #
# Code to reproduce Figure 07 #
#-----------------------------#


#### setup ####

library(tidyverse)
library(patchwork)

# Change ggplot2 default aesthetics
theme_set(theme_bw() + 
            theme(panel.grid = element_blank(), text = element_text(size = 8)))


#### Read data ####

# QTL scan results
qtl_results <- read_csv("./data_processed/qtl_magic/qtl2_scans_lm.csv")
qtl_thresh <- read_csv("./data_processed/qtl_magic/qtl2_scans_lm_perm.csv") %>% 
  filter(level == 0.05)

# Linear mixed model results
qtl_lmm <- read_csv("./data_processed/qtl_magic/qtl_scans_lmm.csv") %>% 
  select(-LOD_full) %>% 
  gather("test", "LOD", matches("LOD")) %>% 
  mutate(test = str_remove(test, "LOD_"))

qtl_lmm_thresh <- read_csv("./data_processed/qtl_magic/qtl2_scans_lmm_perm.csv") %>% 
  select(-LOD_full, -seed) %>% 
  summarise_all(list(threshold = ~ quantile(., probs = 0.95))) %>% 
  gather("test", "threshold") %>% 
  mutate(test = str_remove(test, "LOD_")) %>% 
  separate(test, c("test", "level"), sep = "_") 



#### make plot ####


qtl_results %>% 
  filter(grepl("height", trait)) %>% 
  ggplot(aes(bp/1e6, LOD)) +
  geom_hline(data = qtl_thresh %>% filter(grepl("height", trait)),
             aes(yintercept = threshold), colour = "grey") +
  geom_line() +
  facet_grid(toupper(nitrate) ~ chromosome, scales = "free", space = "free") +
  labs(x = "Mb") +
  theme(panel.grid = element_blank(), legend.position = "none") +
  scale_y_continuous(breaks = seq(0, 20, 4)) +
  scale_x_continuous(breaks = seq(0, 30, 10))


