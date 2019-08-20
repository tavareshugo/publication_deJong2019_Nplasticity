#-----------------------------#
#     de Jong et al 2019      #
# Code to reproduce Figure S7 #
#-----------------------------#

#### setup ####

# Load packages
library(tidyverse)
library(patchwork)

# Change ggplot2 default aesthetics
theme_set(theme_bw() + 
            theme(panel.grid = element_blank(), text = element_text(size = 10)))


#### read data ####

# These were obtained from GCTA, and compiled from two separate scripts
# ./scripts/shell/04_gwas_gcta_250k.sh 
# ./scripts/shell/04_gwas_gcta_imputed.sh 
her <- read_csv("./data_processed/gwas_accessions/compiled_heritabilities.csv")


# Tidy
her <- her %>% 
  filter(trait %in% c("bolt", "height", "totalbr", "totalsil")) %>% 
  mutate(nitrate = ifelse(nitrate == "D", "Plasticity", toupper(nitrate)),
         err_up = ifelse(her_var + her_se < 1, her_var + her_se, 1),
         err_lo = ifelse(her_var - her_se > 0, her_var - her_se, 0)) %>% 
  mutate(trait = case_when(trait == "totalbr" ~ "branches",
                         trait == "height" ~ "height", 
                         trait == "totalsil" ~ "siliques", 
                         trait == "bolt" ~ "flowering")) %>% 
  # add senescence information for correct axis label
  mutate(trait = ifelse(stage == "sen", paste0(trait, "\n(sen)"), trait)) %>% 
  # factorise for right order
  mutate(trait = factor(trait, 
                        levels = c("flowering", "height", "branches", "branches\n(sen)", "siliques\n(sen)")))


#### Make plot ####


pdf("./figures/figureS7.pdf", width = 4.5, height = 3)
her %>% 
  ggplot(aes(trait, her_var, fill = nitrate)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), aes(colour = nitrate)) +
  geom_errorbar(aes(ymax = err_up, ymin = err_lo), width = 0, position = position_dodge(width = 0.9)) +
  facet_grid(set ~ .) +
  scale_fill_manual(values = c("grey48", "grey", "black")) +
  scale_colour_manual(values = c("grey48", "grey", "black")) +
  labs(fill = "Nitrate", colour = "Nitrate", x = "Trait", y = expression(h[GWAS]^2)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
# 
# her %>% 
#   filter(trait %in% c("bolt", "height", "totalbr", "totalsil")) %>% 
#   mutate(trait = ifelse(stage == "sen", paste0(trait, "\n(sen)"), trait),
#          nitrate = ifelse(nitrate == "D", "Plasticity", toupper(nitrate)),
#          err_up = ifelse(her_var + her_se < 1, her_var + her_se, 1),
#          err_lo = ifelse(her_var - her_se > 0, her_var - her_se, 0)) %>% 
#   ggplot(aes(trait, her_var, fill = nitrate)) +
#   geom_bar(stat = "identity", position = position_dodge(0.9), aes(colour = nitrate)) +
#   #geom_point(aes(colour = nitrate), position = position_dodge(0.9)) +
#   geom_errorbar(aes(ymax = err_up, ymin = err_lo), width = 0, position = position_dodge(width = 0.9)) +
#   facet_grid(set ~ .) +
#   scale_fill_manual(values = c("grey48", "grey", "black")) +
#   scale_colour_manual(values = c("grey48", "grey", "black")) +
#   theme(panel.border = element_rect(colour = "black", fill=NA)) +
#   labs(fill = "Nitrate")
# 
# 
# 
# 
