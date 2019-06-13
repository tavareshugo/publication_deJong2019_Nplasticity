########## Permute LMM model with R/qtl2 imputations ############

# This scripts runs a linear mixed model on each marker using permuted genotypes
# The script is meant to run non-interactively, so that it can operate in parallel


#----------------------------------#
#### Read user provided options ####
#----------------------------------#

# Load packages silently
suppressMessages({
  if(!require(optparse) | !require(qtl2) | !require(magrittr) | !require(lme4)){
    stop("please install the `optparse`, `qtl2`, `lme4` and `magrittr` packages.")
  }
})

# Make list of options
option_list = list(
  make_option(c("--workdir"), type = "character", default = NULL, 
              help = "Working directory for the project", metavar = "dir"),
  make_option(c("--seed"), type = "character", default = NULL, 
              help = "An integer to set the seed for reproducible pseudo-random \
              number generator", metavar = "int")
)

# Option parser
opt_parser = OptionParser(option_list=option_list,
                          description = "
                          Script description...")
opt = parse_args(opt_parser)

# Check that working directory exists and set it
if(!dir.exists(opt$workdir)) stop("Working directory does not exist:", opt$workdir)
setwd(opt$workdir)



#-----------------#
#### Read data ####
#-----------------#

data_genoprob <- "./data_processed/qtl_magic/magic_qtl2_genoprob.rds"
data_pheno <- "./data_processed/phenotypes/magic_individual.rds"

if(!file.exists(data_genoprob) | !file.exists(data_pheno)){
  stop("Could not find all data files. Please run '01a_process_magic_data.R' and '03_magic_qtl2_scans.R' scripts first.")
}

# QTL2 genotype data
magic_qtl2_gen <- readRDS(data_genoprob)

# Phenotype data
pheno <- read_rds(data_pheno)

# Subset to retain only those that flower early (contained in the genotype object)
early_magic <- rownames(magic_qtl2_gen[[1]])
pheno <- subset(pheno, id_kover %in% early_magic)



#-------------------------#
#### Prepare genotypes ####
#-------------------------#

# # Calculate genotype probabilities
# magic_qtl2_gen <- calc_genoprob(magic_qtl2, error_prob = 0.01) %>% 
#   genoprob_to_alleleprob()

# Shuffle genotypes 
## doing it this way ensures that replicates of the same MAGIC line have the same genotype
set.seed(opt$seed)
pheno$permuted_id <- plyr::mapvalues(pheno$id_kover, 
                                     from = early_magic, 
                                     to = sample(early_magic))


#-------------------------------#
#### Prepare model variables ####
#-------------------------------#

# Get genotypes using the permuted IDs
permuted_gen <- magic_qtl2_gen[pheno$permuted_id, ]

# Variables
PHEN <- pheno$totalbr  # phenotype
COV <- pheno$nitrate   # nitrate covariate
ID <- pheno$id_kover   # magic line ID



#-----------------------#
##### fitting models ####
#-----------------------#

# Fit null model (this is the same for every marker)
fit_null <- lme4::lmer(PHEN ~ COV + (COV|ID))

# fit genetic models for each marker
lod_scores <- lapply(permuted_gen, function(chr){ 
  lapply(dimnames(chr)[[3]], function(marker){
    
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
})



#----------------------#
#### Prepare output ####
#----------------------#

# Bind the results to a table
lod_scores <- lod_scores %>% 
  unlist(recursive = FALSE) %>% 
  dplyr::bind_rows()

# Get the maximum LOD observed across the permutations
# (the seed number is added for reproducibility)
max_lod <- data.frame(seed = opt$seed,
                      LOD_full = max(lod_scores$LOD_full),
                      LOD_GxE = max(lod_scores$LOD_GxE),
                      LOD_G = max(lod_scores$LOD_G))

# Save into a file 
outfile <- "./data_processed/qtl_magic/qtl2_scans_lmm_perm.csv"

# If file exists append, otherwise write new file (to ensure we get column names)
if(file.exists(outfile)){
  readr::write_csv(max_lod,
                   path = outfile,
                   append = TRUE)
} else {
  readr::write_csv(max_lod,
                   path = outfile)
}

## Deprecated - writing individual files
# readr::write_csv(max_lod,
#                  path = paste0("./data/qtl_magic/permuted_lod_seed", opt$seed, ".csv"))
