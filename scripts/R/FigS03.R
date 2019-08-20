#-----------------------------#
#     de Jong et al 2019      #
# Code to reproduce Figure S3 #
#-----------------------------#

#### setup ####

library(tidyverse)
library(patchwork)

# Change ggplot2 default aesthetics
theme_set(theme_bw() + 
            theme(panel.grid = element_blank(), text = element_text(size = 10)))


#### Read data ####

# heritabilities and bootstrap intervals calculated from other scripts
# ./scripts/R/02a_fit_lmm_magics.R
her_magic <- read_csv("./data_processed/lmm_traits/magic_trait_heritabilities.csv")

# ./scripts/R/02a_fit_lmm_accessions.R
her_accessions <- read_csv("./data_processed/lmm_traits/accessions_trait_heritabilities.csv")

# join and tidy
hers <- bind_rows(her_magic, her_accessions) %>% 
  # rename set
  mutate(set = ifelse(set == "acc", "Accessions", "MAGIC")) %>% 
  # factorise for right order
  mutate(set = factor(set, levels = c("MAGIC", "Accessions"))) %>% 
  # rename traits
  mutate(trait = case_when(trait == "totalbr" ~ "branches",
                           trait == "height" ~ "height", 
                           trait == "totalsil" ~ "siliques", 
                           trait == "bolt" ~ "flowering")) %>% 
  # add senescence information for correct axis label
  mutate(trait = ifelse(stage == "sen", paste0(trait, "\n(sen)"), trait)) %>% 
  # factorise for right order
  mutate(trait = factor(trait, 
                        levels = c("flowering", "height", "branches", "branches\n(sen)", "siliques\n(sen)")))

#### make plots ####

pdf("./figures/figureS3.pdf", width = 4.5, height = 2.5)
hers %>%
  ggplot(aes(trait, her_i, fill = nitrate)) + 
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = her_i_lo, ymax = her_i_hi), width = 0, position = position_dodge(0.9)) +
  facet_grid(flowering ~ set, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c("grey48", "grey")) +
  labs(x = "Trait", y = "Broad-sense Heritability") +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
