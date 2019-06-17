#-------------------------------#
#### Fit linear mixed models ####
#-------------------------------#

# This script fits linear mixed models to the trait data to extract variance components
# and get estimates of broad-sense heritability for each trait


#
# Setup ----
#
# MAGIC data
data_file <- "./data_processed/phenotypes/magic_individual.rds"

if(!file.exists(data_file)){
  stop("Could not find", data_file)
}

# Create output directory
dir.create("./data_processed/lmm_traits/", recursive = TRUE)


# Accession data
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
source("./scripts/R/functions/extractVarsLmm.R")
source("./scripts/R/functions/repeatHeritability.R")


#
# Read data ----
#
# Read MAGIC line data
magic_ind <- read_rds(data_file)

# For modeling purposes, set the variable nitrate as a factor, with LN as the reference value
magic_ind$nitrate <- factor(magic_ind$nitrate, levels = c("LN", "HN"))


# Read accessions data
acc_ind_sil <- read_rds(data_sil)
acc_ind_sen <- read_rds(data_sen)

# For modeling purposes, set the variable nitrate as a factor, with LN as the reference value
acc_ind_sil$nitrate <- factor(acc_ind_sil$nitrate, levels = c("LN", "HN"))
acc_ind_sen$nitrate <- factor(acc_ind_sen$nitrate, levels = c("LN", "HN"))


# Scale variables for easier comparison
magic_ind <- magic_ind %>% 
  mutate_at(vars(totalbr, boltLog, height), ~ (. - mean(., na.rm = TRUE))/sd(., na.rm = TRUE))
acc_ind_sil <- acc_ind_sil %>% 
  mutate_at(vars(totalbr, boltLog, height), ~ (. - mean(., na.rm = TRUE))/sd(., na.rm = TRUE))
acc_ind_sen <- acc_ind_sen %>% 
  mutate_at(vars(totalbr, totalsilSqrt), ~ (. - mean(., na.rm = TRUE))/sd(., na.rm = TRUE))



#
# Fit LMM ----
#
full_lmm <- list(
  magic_height = lmer(height ~ nitrate + (nitrate|id_leyser), 
                      data = magic_ind),
  magic_bolt = lmer(boltLog ~ nitrate + (nitrate|id_leyser), 
                    data = magic_ind),
  magic_totalbr = lmer(totalbr ~ nitrate + (nitrate|id_leyser), 
                       data = magic_ind),
  acc_height = lmer(height ~ nitrate + (nitrate|id_leyser), 
                    data = acc_ind_sil),
  acc_bolt = lmer(boltLog ~ nitrate + (nitrate|id_leyser), 
                  data = acc_ind_sil),
  acc_totalbr = lmer(totalbr ~ nitrate + (nitrate|id_leyser), 
                     data = acc_ind_sil),
  acc_totalbr_sen = lmer(totalbr ~ nitrate + (nitrate|id_leyser), 
                         data = acc_ind_sen),
  acc_totalsil_sen = lmer(totalsilSqrt ~ nitrate + (nitrate|id_leyser), 
                          data = acc_ind_sen)
)

# Models with no GxE interaction
reduced_lmm <- list(
  height = lmer(height ~ nitrate + (1|id_leyser), 
                data = magic_ind),
  bolt = lmer(boltLog ~ nitrate + (1|id_leyser), 
              data = magic_ind),
  totalbr = lmer(totalbr ~ nitrate + (1|id_leyser), 
                 data = magic_ind),
  acc_height = lmer(height ~ nitrate + (1|id_leyser), 
                    data = acc_ind_sil),
  acc_bolt = lmer(boltLog ~ nitrate + (1|id_leyser), 
                  data = acc_ind_sil),
  acc_totalbr = lmer(totalbr ~ nitrate + (1|id_leyser), 
                     data = acc_ind_sil),
  acc_totalbr_sen = lmer(totalbr ~ nitrate + (1|id_leyser), 
                         data = acc_ind_sen),
  acc_totalsil_sen = lmer(totalsilSqrt ~ nitrate + (1|id_leyser), 
                          data = acc_ind_sen)
)


#
# Tidy output ----
#
tidyLmm <- function(model, model2){
  # Extract the random terms (using our custom function)
  random_terms <- model %>% 
    extractVarsLmm() %>% 
    mutate(pct_ln = var_ln/(var_hn + var_ln + var_plas + var_res)*100,
           pct_hn = var_hn/(var_hn + var_ln + var_plas + var_res)*100,
           pct_plas = var_plas/(var_hn + var_ln + var_plas + var_res)*100,
           pct_res = var_res/(var_hn + var_ln + var_plas + var_res)*100) %>% 
    select(var_ln, var_hn, var_plas, var_res, cor_ln_plas, cor_ln_hn, pct_ln, pct_hn, pct_plas, pct_res) %>% 
    gather("term", "estimate") %>% 
    mutate(estimate = round(estimate, 2))
  
  # Extract fixed terms
  fixed_terms <- model %>% 
    broom::tidy(effects = "fixed") %>% 
    select(term, estimate) %>% 
    mutate(estimate = round(estimate, 2))
  
  # Compare with reduced model
  model_comp <- tibble(
    term = c("delta_AIC", "p-value"),
    estimate = c(round(AIC(model) - AIC(model2), 2), 
                 (anova(model, model2)$`Pr(>Chisq)`[2]))
  )
  
  # output
  bind_rows(fixed_terms, random_terms, model_comp) %>% 
    mutate(term = fct_inorder(term))
}


map2(full_lmm, reduced_lmm, tidyLmm) %>% 
  bind_rows(.id = "trait") %>% 
  mutate(estimate = prettyNum(estimate)) %>% 
  spread(trait, estimate) %>% 
  write_csv("./figures/TableS2_tidy_lmm.csv")
