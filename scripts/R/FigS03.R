#-----------------------------#
#     de Jong et al 2019      #
# Code to reproduce Figure S3 #
#-----------------------------#

#### setup ####

library(tidyverse)
library(patchwork)

# Change ggplot2 default aesthetics
theme_set(theme_bw() + 
            theme(panel.grid = element_blank(), text = element_text(size = 8)))


#### Read data ####

# heritabilities and bootstrap intervals calculated from other scripts
# ./scripts/R/02a_fit_lmm_magics.R
her_magic <- read_csv("./data_processed/lmm_traits/magic_trait_heritabilities.csv")

# ./scripts/R/02a_fit_lmm_accessions.R
her_accessions <- read_csv("./data_processed/lmm_traits/accessions_trait_heritabilities.csv")

# join
hers <- bind_rows(her_magic, her_accessions)


#### make plots ####

pdf("./figures/figureS3.pdf", width = 4.5, height = 2.5)
hers %>%
  mutate(set = factor(set, levels = c("magic", "acc")),
         trait = ifelse(stage == "sen", paste0(trait, " (sen)"), trait)) %>% 
  ggplot(aes(trait, her_i, fill = nitrate)) + 
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = her_i_lo, ymax = her_i_hi), width = 0, position = position_dodge(0.9)) +
  facet_grid(flowering ~ set, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c("grey48", "grey")) +
  labs(x = "Trait", y = "Broad-sense Heritability") +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank())
dev.off()
