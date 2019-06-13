#-------------------------------#
#### Fit linear mixed models ####
#-------------------------------#

# This script fits linear mixed models to the trait data to extract variance components
# and get estimates of broad-sense heritability for each trait


#
# Setup ----
#
data_file <- "./data_processed/phenotypes/magic_individual.rds"

if(!file.exists(data_file)){
  stop("Could not find", data_file, ". Please run 01a_process_magic_data.R script first.")
}

# Create output directory
dir.create("./data_processed/lmm_traits/", recursive = TRUE)


#
# Load packages ----
#
library(tidyverse)
library(lme4)

# Custom functions
source("./scripts/R/functions/extractVarsLmm.R")
source("./scripts/R/functions/repeatHeritability.R")


#
# Read data ----
#
# Read MAGIC line data
magic_ind <- read_rds(data_file)

# For modeling purposes, set the variable nitrate as a factor, with LN as the reference value
magic_ind$nitrate <- factor(magic_ind$nitrate, levels = c("LN", "HN"))

# Get list of MAGIC lines that flower < 25 days on LN
early_magic <- magic_ind %>% 
  filter(nitrate == "LN") %>% 
  group_by(id_leyser) %>% 
  filter(mean(bolt, na.rm = TRUE) < 25) %>% 
  pull(id_leyser) %>% 
  unique()


#
# Fit LMM ----
#
magic_lmm <- list(
  height = lmer(height ~ nitrate + (nitrate|id_leyser), 
                data = magic_ind, 
                contrasts = list(nitrate = contr.helmert(2))),
  bolt = lmer(boltLog ~ nitrate + (nitrate|id_leyser), 
              data = magic_ind, 
              contrasts = list(nitrate = contr.helmert(2))),
  totalbr = lmer(totalbr ~ nitrate + (nitrate|id_leyser), 
                 data = magic_ind, 
                 contrasts = list(nitrate = contr.helmert(2)))
)


#
# Variance Components ----
#
# Extract variance components using custom function
magic_trait_vars <- map_df(magic_lmm, extractVarsLmm, .id = "trait") %>% 
  gather("component", "variance", var_ln, var_hn, var_plas, var_res, var_nitrate) %>% 
  mutate(set = "MAGIC lines")


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
  magic_bolt_sil_all = magic_ind %>% select(nitrate, id_leyser, boltLog),
  magic_bolt_sil_early = magic_ind %>% filter(id_leyser %in% early_magic) %>% select(nitrate, id_leyser, boltLog),
  magic_height_sil_all = magic_ind %>% select(nitrate, id_leyser, height),
  magic_height_sil_early = magic_ind %>% filter(id_leyser %in% early_magic) %>% select(nitrate, id_leyser, height),
  magic_totalbr_sil_all = magic_ind %>% select(nitrate, id_leyser, totalbr),
  magic_totalbr_sil_early = magic_ind %>% filter(id_leyser %in% early_magic) %>% select(nitrate, id_leyser, totalbr)
)

# Loop through list to calculate heritabilities per nitrate
# takes a while - grab a cup of tea!
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
magic_trait_vars %>% 
  select(trait, component, cor_ln_hn, intercept, variance, set) %>% 
  write_csv("./data_processed/lmm_traits/magic_trait_variances.csv")

# Save heritability table
hers %>% 
  write_csv("./data_processed/lmm_traits/magic_trait_heritabilities.csv")


