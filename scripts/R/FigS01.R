#-----------------------------#
#     de Jong et al 2019      #
# Code to reproduce Figure S1 #
#-----------------------------#


#### setup ####

# Load packages
library(tidyverse)
library(patchwork)

# Change ggplot2 default aesthetics
theme_set(theme_classic() + 
            theme(text = element_text(size = 8)))

# Custom function 
source("./scripts/R/functions/corrLabel.R")


#
# Read data ----
#
# Read Accession data - silique stage
acc_sum_sil <- read_rds("./data_processed/phenotypes/accessions_summarised_silique.rds")

# Read Accession data - senescence stage
acc_sum_sen <- read_rds("./data_processed/phenotypes/accessions_summarised_senescence.rds")

# Join the silique and senescence tables
acc_sum <- inner_join(acc_sum_sil, acc_sum_sen, by = c("id_leyser", "nitrate")) %>% 
  select(id_leyser, nitrate, matches("totalbr_mean")) %>% 
  rename(totalbr_sil = totalbr_mean.x,
         totalbr_sen = totalbr_mean.y)


#
# Correlation ----
# 
# Correlation between traits at silique and senescence stages
cor_sil_sen <- data.frame(nitrate = c("HN", "LN"),
                          cor = c(corLabel(acc_sum$totalbr_sil[acc_sum$nitrate == "HN"], 
                                           acc_sum$totalbr_sen[acc_sum$nitrate == "HN"], 
                                           TRUE,
                                           use = "complete.obs"),
                                  corLabel(acc_sum$totalbr_sil[acc_sum$nitrate == "LN"], 
                                           acc_sum$totalbr_sen[acc_sum$nitrate == "LN"], 
                                           TRUE, 
                                           use = "complete.obs")))


#
# plot ----
#

# Figure
#pdf("./figures/figureS1.pdf", width = 5.2, height = 5)
acc_sum %>% 
  ggplot(aes(totalbr_sil, totalbr_sen)) +
  geom_point() +
  geom_text(data = cor_sil_sen, x = 0, y = 7, aes(label = cor), size = 2.5, hjust = 0) +
  geom_abline(linetype = "dashed") +
  facet_wrap(~ nitrate) +
  labs(x = "Total Branches (silique)", y = "Total Branches (senescence)")
#dev.off()
