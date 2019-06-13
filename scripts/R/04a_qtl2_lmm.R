########## LMM model with R/qtl2 imputations ############

# This scripts runs a linear mixed model on each marker using individual data

library(qtl2)
library(magrittr)


#-----------------#
#### Read data ####
#-----------------#

data_cross <- "./data_processed/qtl_magic/magic_qtl2_cross2.rds"
data_genoprob <- "./data_processed/qtl_magic/magic_qtl2_genoprob.rds"
data_pheno <- "./data_processed/phenotypes/magic_individual.rds"

if(!file.exists(data_cross) | !file.exists(data_genoprob) | !file.exists(data_pheno)){
  stop("Could not find all data files. Please run '01a_process_magic_data.R' and '03_magic_qtl2_scans.R' scripts first.")
}

# QTL2 data
magic_qtl2 <- readRDS(data_cross)
magic_qtl2_gen <- readRDS(data_genoprob)

# Phenotype data
pheno <- read_rds(data_pheno)

# Subset to retain only those that flower early (contained in the genotype object)
#early_magic <- rownames(magic_qtl2$pheno)[magic_qtl2$pheno[, "bolt_mean_ln"] < 25]
pheno <- subset(pheno, id_kover %in% ind_ids(magic_qtl2))



# #-------------------------#
# #### Prepare genotypes ####
# #-------------------------#
# 
# # Calculate genotype probabilities
# magic_qtl2_gen <- calc_genoprob(magic_qtl2, error_prob = 0.01) %>% 
#   genoprob_to_alleleprob()



#-------------------------------#
#### Prepare model variables ####
#-------------------------------#

# Get genotypes for all the MAGIC lines (including replicates)
all_gen <- magic_qtl2_gen[pheno$id_kover, ]

# Variables
PHEN <- pheno$totalbr  # phenotype
COV <- pheno$nitrate   # nitrate covariate
ID <- pheno$id_kover   # magic line ID



##### fitting models ####
# Fit null model (this is the same for every marker)
fit_null <- lme4::lmer(PHEN ~ COV + (COV|ID))

# fit genetic models for each marker
lod_scores <- lapply(all_gen, function(chr){ 
  lods <- lapply(dimnames(chr)[[3]], function(marker){
    
    # Get genotype matrix
    GEN <- chr[, , marker]
    
    # convert it to binary values
    GEN <- (GEN > 0.5)*1
    
    tryCatch(
      {
        # Fit each of the models
        # Full model including interaction term
        fit_full <- lme4::lmer(PHEN ~ COV + GEN + GEN:COV + (COV|ID))
        
        # Genetic model including genetic but no interaction term
        fit_gen <- lme4::lmer(PHEN ~ COV + GEN + (COV|ID))
        
        # Output a table with LOD scores for each
        data.frame(
          LOD_full = (logLik(fit_full) - logLik(fit_null))/log(10),
          LOD_GxE = (logLik(fit_full) - logLik(fit_gen))/log(10),
          LOD_G = (logLik(fit_gen) - logLik(fit_null))/log(10)
        )
        
      }, error=function(err) data.frame(LOD_full = NA, LOD_GxE = NA, LOD_G = NA))    
  })
  
  names(lods) <- dimnames(chr)[[3]]
  
  return(lods)
})


# fit genetic model but return result of F-test
f_tests <- lapply(all_gen, function(chr){ 
  result <- lapply(dimnames(chr)[[3]], function(marker){
    
    # Get genotype matrix
    GEN <- chr[, , marker]
    
    # convert it to binary values
    GEN <- (GEN > 0.5)*1
    
    tryCatch(
      {
        fit1 <- lmerTest::lmer(PHEN ~ COV + GEN + GEN:COV + (COV|ID))
        anova(fit1, type = 2) %>% broom::tidy()
      }, error=function(err) tibble(term = NA, sumsq = NA, meansq = NA, NumDF = NA, DenDF = NA, statistic = NA, p.value = NA))
    
  })
  
  names(result) <- dimnames(chr)[[3]]
  
  return(result)
})



#----------------------#
#### Prepare output ####
#----------------------#

# Bind the results of the scan to tidy tables
names(lod_scores) <- NULL
lod_scores <- lod_scores %>% 
  unlist(recursive = FALSE) %>% 
  dplyr::bind_rows(.id = "marker")

names(f_tests) <- NULL
f_tests <- f_tests %>% 
  unlist(recursive = FALSE) %>% 
  dplyr::bind_rows(.id = "marker")

# Get marker information from qtl2 object
markers <- lapply(magic_qtl2$pmap, function(i){
  data.frame(marker = names(i), bp = as.numeric(i), stringsAsFactors = FALSE)
}) %>% 
  dplyr::bind_rows(.id = "chrom")

# Join the marker and test tables
lod_scores <- merge(markers, lod_scores, by = "marker", all = TRUE)

f_tests <- merge(markers, f_tests, by = "marker", all = TRUE)

# Save output
readr::write_csv(lod_scores,
                 "./data_processed/qtl_magic/qtl_scans_lmm.csv")

readr::write_csv(f_tests,
                 "./data_processed/qtl_magic/qtl_scans_lmm_ftest.csv")

