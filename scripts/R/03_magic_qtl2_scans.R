#------------------------#
#### Run R/qtl2 scans ####
#------------------------#

# This script does the QTL scans in MAGIC lines

#
# Setup ----
#
data_file <- "./data_processed/phenotypes/magic_plasticity.rds"

if(!file.exists(data_file)){
  stop("Could not find", data_file, ". Please run 01a_process_magic_data.R script first.")
}

# Create output directory
dir.create("./data_processed/qtl_magic/", recursive = TRUE)


#
# Load packages ----
#

library(qtl2)
library(dplyr)
library(readr)
library(purrr)
library(tidyr)
library(stringr)

# Load custom function to add phenotypes to R/qtl2 object
source("scripts/R/functions/add_pheno.R")


#-----------------#
#### Read data ####
#-----------------#

# Read phenotypes
pheno <- read_rds(data_file)

# Retain only phenotypes of interest for early-flowering individuals
pheno <- pheno %>% 
  select(ID = id_kover, 
         height_mean_hn, height_mean_ln, 
         bolt_mean_hn, bolt_mean_ln,
         totalbr_mean_hn, totalbr_mean_ln, totalbr_mean_D) %>% 
  filter(bolt_mean_ln < 25)

# Get cross2 object from our package
magic_qtl2 <- atMAGIC::kover2009

# Alternative cross2 object available from:
#magic_qtl2 <- read_cross2("https://raw.githubusercontent.com/tavareshugo/qtl2data/ArabMAGIC/ArabMAGIC/arabmagic.zip")

# add phenotypes to it (retain only phenotyped individuals)
magic_qtl2 <- add_pheno(magic_qtl2, pheno, retain_all = FALSE)

# Calculate genotype probabilities assuming a 1% probability of error
magic_qtl2_gen <- calc_genoprob(magic_qtl2, error_prob = 0.01, cores = 2)



#
# QTL scans ----
#
# Run scan for all traits
qtl2_all_scans <- scan1(magic_qtl2_gen, magic_qtl2$pheno)

# Run permutation scans to get thresholds
qtl2_thresholds <- scan1perm(magic_qtl2_gen, magic_qtl2$pheno, n_perm = 1000, cores = 2) %>% 
  summary(alpha = c(0.01, 0.05, 0.1))

# Scan with covariates
qtl2_scan_covar <- list(
  totalbr_hn = scan1(magic_qtl2_gen, magic_qtl2$pheno[, "totalbr_mean_hn"], 
                     addcovar = magic_qtl2$pheno[, "bolt_mean_hn"]),
  totalbr_ln = scan1(magic_qtl2_gen, magic_qtl2$pheno[, "totalbr_mean_ln"], 
                     addcovar = magic_qtl2$pheno[, "bolt_mean_ln"]),
  totalbr_D = scan1(magic_qtl2_gen, magic_qtl2$pheno[, "totalbr_mean_D"], 
                    addcovar = magic_qtl2$pheno[, "bolt_mean_ln"])
)

# Run permutation scans on models with covariate
qtl2_thresholds_covar <- list(
  totalbr_hn = scan1perm(magic_qtl2_gen, magic_qtl2$pheno[, "totalbr_mean_hn"], 
                         addcovar = magic_qtl2$pheno[, "bolt_mean_hn"], n_perm = 1000, cores = 2),
  totalbr_ln = scan1perm(magic_qtl2_gen, magic_qtl2$pheno[, "totalbr_mean_ln"], 
                         addcovar = magic_qtl2$pheno[, "bolt_mean_ln"], n_perm = 1000, cores = 2),
  totalbr_D = scan1perm(magic_qtl2_gen, magic_qtl2$pheno[, "totalbr_mean_D"], 
                        addcovar = magic_qtl2$pheno[, "bolt_mean_ln"], n_perm = 1000, cores = 2)
) %>% 
  map(summary, alpha = c(0.01, 0.05, 0.1))



#
# Tidy scans ----
#
# Extract marker information into data.frame
markers <- map_dfr(magic_qtl2$pmap, function(i){
  tibble(marker = names(i), bp = as.numeric(i))
}, .id = "chromosome")

# Tidy threshold results
qtl2_thresholds <- qtl2_thresholds %>% 
  as.data.frame() %>% 
  as_tibble(rownames = "level") %>% 
  gather("trait", "threshold", -level) %>% 
  mutate(trait = str_remove(trait, "_mean")) %>% 
  separate(trait, c("trait", "nitrate"), sep = "_")

# Tidy threshold results with covariate models
qtl2_thresholds_covar <- qtl2_thresholds_covar %>% 
  map(function(i) as.data.frame(i) %>% as_tibble(rownames = "level")) %>% 
  bind_rows(.id = "trait") %>% 
  rename(threshold = pheno1) %>% 
  separate(trait, c("trait", "nitrate"), sep = "_")

# Tidy scan results
qtl2_all_scans <- qtl2_all_scans %>% 
  as.data.frame() %>% 
  as_tibble(rownames = "marker") %>% 
  left_join(markers, by = "marker") %>% 
  gather("trait", "LOD", -marker, -chromosome, -bp) %>% 
  mutate(trait = str_remove(trait, "_mean")) %>% 
  separate(trait, c("trait", "nitrate"), sep = "_")

# Tidy scans with covariate
qtl2_scan_covar <- qtl2_scan_covar %>% 
  map_dfr(function(i){
    i %>% 
      as.data.frame() %>% 
      as_tibble(rownames = "marker") %>% 
      left_join(markers, by = "marker")
  }, .id = "trait") %>% 
  separate(trait, c("trait", "nitrate"), sep = "_") %>% 
  rename(LOD = pheno1)


#
# Save output ----
#
# Bind both covariate and simple models
bind_rows(no = qtl2_thresholds, yes = qtl2_thresholds_covar, .id = "covariate") %>% 
  write_csv("./data_processed/qtl_magic/qtl2_scans_lm_perm.csv")

bind_rows(no = qtl2_all_scans, yes = qtl2_scan_covar, .id = "covariate") %>% 
  write_csv("./data_processed/qtl_magic/qtl2_scans_lm.csv")

# Save R/qtl2 cross2 objects also
write_rds(magic_qtl2, "./data_processed/qtl_magic/magic_qtl2_cross2.rds")
write_rds(magic_qtl2_gen, "./data_processed/qtl_magic/magic_qtl2_genoprob.rds")
