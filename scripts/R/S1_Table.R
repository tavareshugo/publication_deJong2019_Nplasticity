####
# S1 Table - variance explained by QTL
###

library(qtl2)
library(tidyverse)
library(patchwork)

source("./functions/repeatHeritability.R")

# Change ggplot2 default aesthetics
theme_set(theme_bw() + 
            theme(panel.grid = element_blank(), 
                  text = element_text(size = 8)))


#### Read data ####

# QTL scan results
qtl_results <- read_csv("../../data/qtl_magic/qtl2_scans_lm.csv")
qtl_thresh <- read_csv("../../data/qtl_magic/qtl2_scans_lm_perm.csv") %>% 
  filter(level == 0.05)

# Linear mixed model results
qtl_lmm <- read_csv("../../data/qtl_magic/qtl_scans_lmm.csv") %>% 
  select(-LOD_full) %>% 
  gather("test", "LOD", matches("LOD")) %>% 
  mutate(test = str_remove(test, "LOD_"))
qtl_lmm_thresh <- read_csv("../../data/qtl_magic/qtl2_scans_lmm_perm.csv") %>% 
  select(-LOD_full, -seed) %>% 
  summarise_all(funs(threshold = quantile(., probs = 0.95))) %>% 
  gather("test", "threshold") %>% 
  mutate(test = str_remove(test, "LOD_")) %>% 
  separate(test, c("test", "level"), sep = "_") 

# Genotype objects
magic_gen <- read_rds("../../data/qtl_magic/magic_qtl2_cross2.rds")
magic_prob <- read_rds("../../data/qtl_magic/magic_qtl2_genoprob.rds")

# Phenotype data
pheno <- readr::read_csv("../../data/phenotypes/magic_individual.csv")

# Subset to retain only those that flower early (contained in the genotype object)
pheno <- subset(pheno, id_kover %in% ind_ids(magic_gen))


#### Get peaks above thresholds ####

# For simple linear model
sig_lm <- qtl_results %>% 
  # Join with thresholds table 
  full_join(qtl_thresh, by = c("covariate", "trait", "nitrate")) %>% 
  # Keep those above threshold
  filter(LOD >= threshold) %>%
  # Take top LOD score for each sign level, test and chromosome 
  group_by(chromosome, trait, nitrate, covariate) %>% 
  arrange(desc(LOD)) %>% 
  slice(1) %>% 
  ungroup()

# For mixed model
sig_lmm <- qtl_lmm %>% 
  # Join with thresholds table 
  full_join(qtl_lmm_thresh, by = "test") %>% 
  # Keep those above threshold
  filter(LOD >= threshold) %>%
  # Take top LOD score for each sign level, test and chromosome 
  group_by(level, test, chrom) %>% 
  arrange(desc(LOD)) %>% 
  slice(1) %>% 
  ungroup()

# Combine all markers
peak_markers <- bind_rows(sig_lm, sig_lmm)


#### Variance explained by individual QTL ####

# Add column of variance explained by QTL
sig_lm$var_explained <- NA

# First look at non-plasticity traits
for(i in 1:nrow(sig_lm)){
  if(sig_lm$nitrate[i] == "D") next
  
  # Get only relevant lines
  pheno_subset <- pheno %>% filter(nitrate == toupper(sig_lm$nitrate[i]))
  
  # Variables
  PHEN <- pheno_subset[[sig_lm$trait[i]]]  # phenotype
  ID <- pheno_subset$id_kover   # magic line ID
  
  # Get genotypes for all the MAGIC lines (including replicates)
  all_gen <- magic_prob[ID, ]
  
  # Fit null model and get between-line variance (this is the same for every marker)
  fit_null <- lme4::lmer(PHEN ~ 1 + (1|ID))
  null_var <- lme4::VarCorr(fit_null) %>% as.data.frame %>% pull(vcov)
  null_var <- null_var[1]/sum(null_var)
  
  # Fit QTL model and get between-line variance
  GEN <- all_gen[[sig_lm$chromosome[i]]][, , sig_lm$marker[i]]
  GEN <- (GEN > 0.5)*1
  fit_qtl <- lme4::lmer(PHEN ~ GEN + (1|ID))
  
  qtl_var <- lme4::VarCorr(fit_qtl) %>% as.data.frame %>% pull(vcov)
  qtl_var <- qtl_var[1]/sum(qtl_var)
  
  # Calculate difference to obtained variance explained by QTL
  sig_lm$var_explained[i] <- round((null_var - qtl_var)*100, 2)
  
}


# Now look at GxE trait
for(i in 1:nrow(sig_lm)){
  if(sig_lm$nitrate[i] != "D") next
  
  # All lines are used for this model
  pheno_subset <- pheno
  
  # Variables
  PHEN <- pheno_subset[[sig_lm$trait[i]]]  # phenotype
  COV <- pheno_subset$nitrate
  ID <- pheno_subset$id_kover   # magic line ID
  
  # Get genotypes for all the MAGIC lines (including replicates)
  all_gen <- magic_prob[ID, ]
  
  # Fit null model and get GxE variance (this is the same for every marker)
  fit_null <- lme4::lmer(PHEN ~ 1 + COV + (COV|ID))
  null_var <- lme4::VarCorr(fit_null) %>% as.data.frame %>% filter(is.na(var2)) %>% pull(vcov)
  null_var <- null_var[2]/sum(null_var)
  
  # Fit QTL model and get between-line variance
  GEN <- all_gen[[sig_lm$chromosome[i]]][, , sig_lm$marker[i]]
  GEN <- (GEN > 0.5)*1
  fit_qtl <- lme4::lmer(PHEN ~ COV + GEN + GEN:COV + (COV|ID))
  
  qtl_var <- lme4::VarCorr(fit_qtl) %>% as.data.frame %>% filter(is.na(var2)) %>% pull(vcov)
  qtl_var <- qtl_var[2]/sum(qtl_var)
  
  # Calculate difference to obtained variance explained by QTL
  sig_lm$var_explained[i] <- round((null_var - qtl_var)*100, 2)
  
}

# Results included in S1 table
sig_lm %>% arrange(covariate, trait, nitrate, chromosome, bp)


#### Joint SNP variance explained ####

# Flowering time markers
ft_markers <- c("MN1_23474588", "MN5_4327715", "MN1_24322296", "MN5_4327715")

# 

# Shoot branching markers
sb_markers <- c("MASC00557", "MN5_26708459", 
                "MN1_26278413", "MN3_5910420", 
                "SOC1_461", "MN5_5452600")




peak_markers <- tibble(
  chr = c(1, 5, 1, 1, 5, 1, 3, 2, 5),
  peak = c("FT.HN.1", "FT.HN.5/FT.LN.5", "FT.LN.1", 
           "SB.HN.1", "SB.HN.5", 
           "SB~FT.LN.1", "SB~FT.LN.3", 
           "SB.PL.2", "SB.PL.5"),
  marker = c("MN1_23474588", "MN5_4327715", "MN1_24322296", 
             "MASC00557", "MN5_26708459", 
             "MN1_26278413", "MN3_5910420", 
             "SOC1_461", "MN5_5452600")
)

# Get the genotype for each marker:
peak_geno <- vector("list", length = nrow(peak_markers))
names(peak_geno) <- peak_markers$marker

for(i in seq_len(nrow(peak_markers))){
  # Get the founder probabilities
  GEN <- all_gen[[peak_markers$chr[i]]][, , peak_markers$marker[i]]
  
  # Convert probabilities to most likely founder genotype
  GEN <- (GEN > 0.5)*1
  
  peak_geno[[peak_markers$marker[i]]] <- GEN
}


# Fit flowering time model #

# Null model (no markers)
PHEN <- pheno$bolt
COV <- pheno$nitrate
ID <- pheno$id_kover
ft_null <- lme4::lmer(PHEN ~ 1 + COV + (COV|ID))
ft_null_var <- lme4::VarCorr(ft_null) %>% as.data.frame %>% filter(is.na(var2)) %>% pull(vcov)
ft_null_var <- sum(ft_null_var[1:2])/sum(ft_null_var)

# Fit QTL model and get between-line variance
ft_qtl <- lme4::lmer(PHEN ~ COV + 
                       COV:peak_geno$MN1_23474588 + COV:peak_geno$MN5_4327715 + COV:peak_geno$MN1_24322296 + 
                       (COV|ID))

ft_qtl_var <- lme4::VarCorr(ft_qtl) %>% as.data.frame %>% filter(is.na(var2)) %>% pull(vcov)
ft_qtl_var <- sum(ft_qtl_var[1:2])/sum(ft_qtl_var)

# Calculate difference to obtained variance explained by QTL
round((ft_null_var - ft_qtl_var)*100, 2)


# Fit shoot branching model #

# Null model (no markers)
PHEN <- pheno$totalbr
COV <- pheno$nitrate
ID <- pheno$id_kover
sb_null <- lme4::lmer(PHEN ~ 1 + COV + (COV|ID))
sb_null_var <- lme4::VarCorr(sb_null) %>% as.data.frame %>% filter(is.na(var2)) %>% pull(vcov)
sb_null_var <- sum(sb_null_var[1:2])/sum(sb_null_var)

# Fit QTL model and get between-line variance
sb_qtl <- lme4::lmer(PHEN ~ COV + 
                       COV:peak_geno$MASC00557 + COV:peak_geno$MN5_26708459 + COV:peak_geno$MN1_26278413 + 
                       COV:peak_geno$MN3_5910420 + COV:peak_geno$SOC1_461 + COV:peak_geno$MN5_5452600 +
                       (COV|ID))

sb_qtl_var <- lme4::VarCorr(sb_qtl) %>% as.data.frame %>% filter(is.na(var2)) %>% pull(vcov)
sb_qtl_var <- sum(sb_qtl_var[1:2])/sum(sb_qtl_var)

# Calculate difference to obtained variance explained by QTL
round((sb_null_var - sb_qtl_var)*100, 2)

