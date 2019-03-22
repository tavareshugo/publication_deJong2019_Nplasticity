#-------------------------------#
#### Fit linear mixed models ####
#-------------------------------#

# This script fits linear mixed models to the trait data to extract variance components
# and get estimates of broad-sense heritability for each trait


#
# Setup ----
#
data_sil <- "./data_processed/phenotypes/accessions_individual_silique.rds"
data_sen <- "./data_processed/phenotypes/accessions_individual_senescence.rds"

if(!file.exists(data_sil) | !file.exists(data_sen)){
  stop("Could not find data files.")
}

# Create output directory
dir.create("./data_processed/lmm_traits/", recursive = TRUE)


#
# Load packages ----
#
library(tidyverse)
library(lme4)

# Custom functions
source("scripts/R/functions/extractVarsLmm.R")
source("scripts/R/functions/repeatHeritability.R")


#
# Read data ----
#
# Read accessions data
acc_ind_sil <- read_rds(data_sil)
acc_ind_sen <- read_rds(data_sen)

# For modeling purposes, set the variable nitrate as a factor, with LN as the reference value
acc_ind_sil$nitrate <- factor(acc_ind_sil$nitrate, levels = c("LN", "HN"))
acc_ind_sen$nitrate <- factor(acc_ind_sen$nitrate, levels = c("LN", "HN"))

# Get list of accessions that flower < 25 days on LN
early_acc <- acc_ind_sil %>% 
  filter(nitrate == "LN") %>% 
  group_by(id_leyser) %>% 
  filter(mean(bolt, na.rm = TRUE) < 25) %>% 
  pull(id_leyser) %>% 
  unique()


#
# Fit LMM ----
#
acc_lmm <- list(
  height = lmer(height ~ nitrate + (nitrate|id_leyser), 
                data = acc_ind_sil, 
                contrasts = list(nitrate = contr.helmert(2))),
  bolt = lmer(boltLog ~ nitrate + (nitrate|id_leyser), 
              data = acc_ind_sil, 
              contrasts = list(nitrate = contr.helmert(2))),
  totalbr = lmer(totalbr ~ nitrate + (nitrate|id_leyser), 
                 data = acc_ind_sil, 
                 contrasts = list(nitrate = contr.helmert(2))),
  totalbr_sen = lmer(totalbr ~ nitrate + (nitrate|id_leyser), 
                     data = acc_ind_sen, 
                     contrasts = list(nitrate = contr.helmert(2))),
  totalsil_sen = lmer(totalsilSqrt ~ nitrate + (nitrate|id_leyser), 
                      data = acc_ind_sen, 
                      contrasts = list(nitrate = contr.helmert(2)))
)


#
# Variance Components ----
#
# Extract variance components using custom function
acc_trait_vars <- map_df(acc_lmm, extractVarsLmm, .id = "trait") %>% 
  gather("component", "variance", var_ln, var_hn, var_plas, var_res, var_nitrate) %>% 
  mutate(set = "Accessions")


#
# Heritability ----
#
# Make list of phenotypes and respective genotypes to estimate heritabilities
## naming convention is "population_trait_stage_flowering", where:
## population: either MAGIC lines or Accessions
## trait: bolt, height, total branches, total siliques
## stage: 2-silique or senescence
## flowering: all lines or only early flowering lines (<25 days)
traits <- list(
  acc_bolt_sil_all = acc_ind_sil %>% select(nitrate, id_leyser, boltLog),
  acc_bolt_sil_early = acc_ind_sil %>% filter(id_leyser %in% early_acc) %>% select(nitrate, id_leyser, boltLog),
  acc_height_sil_all = acc_ind_sil %>% select(nitrate, id_leyser, height),
  acc_height_sil_early = acc_ind_sil %>% filter(id_leyser %in% early_acc) %>% select(nitrate, id_leyser, height),
  acc_totalbr_sil_all = acc_ind_sil %>% select(nitrate, id_leyser, totalbr),
  acc_totalbr_sil_early = acc_ind_sil %>% filter(id_leyser %in% early_acc) %>% select(nitrate, id_leyser, totalbr),
  acc_totalbr_sen_all = acc_ind_sen %>% select(nitrate, id_leyser, totalbr),
  acc_totalbr_sen_early = acc_ind_sen %>% filter(id_leyser %in% early_acc) %>% select(nitrate, id_leyser, totalbr),
  acc_totalsil_sen_all = acc_ind_sen %>% select(nitrate, id_leyser, totalsilSqrt),
  acc_totalsil_sen_early = acc_ind_sen %>% filter(id_leyser %in% early_acc) %>% select(nitrate, id_leyser, totalsilSqrt)
)

# Loop through list to calculate heritabilities per nitrate
hers <- map_df(traits, 
               function(x){
                 names(x) <- c("nitrate", "genotype", "phenotype")
                 x %>% 
                   group_by(nitrate) %>% 
                   do(repeatHeritability(.$phenotype, .$genotype, plot_diag = FALSE, 
                                         ci = 0.95, nboot = 1000, ncores = 2))
               }, .id = "group") %>% 
  separate(group, c("set", "trait", "stage", "flowering"), sep = "_")


#
# Save data ----
#
# Save variance tables
acc_trait_vars %>% 
  select(trait, component, cor_ln_hn, intercept, variance, set) %>% 
  write_csv("./data_processed/lmm_traits/accessions_trait_variances.csv")

# Save heritability table
hers %>% 
  write_csv("./data_processed/lmm_traits/accessions_trait_heritabilities.csv")



