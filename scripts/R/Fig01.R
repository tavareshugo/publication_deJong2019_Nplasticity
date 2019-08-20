#-----------------------------#
#     de Jong et al 2019      #
# Code to reproduce Figure 01 #
#-----------------------------#


#### setup ####

# Load packages
library(tidyverse)
library(patchwork)

# Change ggplot2 default aesthetics
theme_set(theme_classic() + 
            theme(text = element_text(size = 10)))


#### Read data ####

# Read variance partitioning data
magic_trait_vars <- read_csv("./data_processed/lmm_traits/magic_trait_variances.csv")
acc_trait_vars <- read_csv("./data_processed/lmm_traits/accessions_trait_variances.csv")
all_trait_vars <- bind_rows(magic_trait_vars, acc_trait_vars)

# Make table of correlations
corr_labels <- all_trait_vars %>% 
  distinct(trait, set, cor_ln_hn) %>% 
  mutate(cor_ln_hn = signif(cor_ln_hn, 2), set = factor(set, levels = c("MAGIC lines", "Accessions")))


#
# plot ----
#
# Panel A
p1 <- all_trait_vars %>% 
  mutate(trait = factor(trait, levels = c("bolt", "height", "totalbr", "totalbr_sen", "totalsil_sen")),
         set = factor(set, levels = c("MAGIC lines", "Accessions"))) %>% 
  group_by(set, trait) %>% 
  mutate(variance = variance/sum(variance)*100,
         component = factor(component, levels = c("var_res", "var_plas", "var_nitrate", "var_ln", "var_hn"))) %>% 
  ungroup() %>% 
  ggplot(aes(trait, variance, fill = component, label = round(variance))) +
  geom_col(colour = "black") +
  facet_grid(~ set, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c("#999999", "#56B4E9", "#009E73", "#E69F00", "#F0E442"),
                    labels = c("Residual", "Genotype\nx\nNitrate", "Nitrate", "Genotype\n(HN)", "Genotype\n(LN)")) +
  #scale_fill_manual(values = c("white", "gray48", "grey", "black")) +
  #geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  geom_text(data = corr_labels, 
            aes(x = trait, y = 104, label = sprintf("%0.2f", cor_ln_hn)),
            inherit.aes = FALSE, size = 3) +
  scale_x_discrete(labels = c("flowering", "height", "branches", 
                              "branches\n(sen)", "siliques\n(sen)")) + 
  labs(y = "% variance", tag = "A") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_blank(),
        legend.text.align = 0.5,
        legend.title = element_blank())

# Panel B
p2 <- all_trait_vars %>% 
  mutate(trait = factor(trait, levels = c("bolt", "height", "totalbr", "totalbr_sen", "totalsil_sen")),
         set = factor(set, levels = c("MAGIC lines", "Accessions"))) %>% 
  group_by(set, trait) %>% 
  mutate(gcv = sqrt(variance)/intercept,
         component = factor(component, levels = c("var_res", "var_plas", "var_nitrate", "var_ln", "var_hn"))) %>% 
  ungroup() %>% 
  ggplot(aes(trait, gcv, fill = component, label = round(gcv))) +
  geom_col(colour = "black", position = "dodge") +
  facet_grid( ~ set, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c("#999999", "#56B4E9", "#009E73", "#E69F00", "#F0E442"),
                    labels = c("Residual", "Genotype\nx\nNitrate", "Nitrate", "Genotype\n(HN)", "Genotype\n(LN)")) +
  # scale_fill_manual(values = c("white", "gray48", "grey", "black"),
  #                   labels = c("Res", "GxE", "E", "G"), name = NULL) +
  labs(y = "Coefficient of variation") +
  scale_x_discrete(labels = c("flowering", "height", "branches", "branches\n(sen)", "siliques\n(sen)")) + 
  labs(y = "CV", tag = "B") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_blank(), legend.position = "none")

# Arrange (and save) figure
pdf("./figures/figure1.pdf", width = 4.5, height = 5)
#tiff("./figures/figure1.tiff", width = 4.5, height = 5, units = "in", compression = "lzw", res = 600)
p1 + p2 + plot_layout(ncol = 1)
dev.off()
