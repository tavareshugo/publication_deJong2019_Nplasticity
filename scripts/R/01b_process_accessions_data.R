#--------------------------#
#### Process trait data ####
#--------------------------#

# Take individual data and
# summarise trait per line
# calculate trait plasticity

#
# Setup ----
#
# Check raw files exist
if(!file.exists("./data/phenotypes/accessions_silique.csv")){
  stop("Could not find ./data/phenotypes/accessions_silique.csv")
}

if(!file.exists("./data/phenotypes/accessions_senescence.csv")){
  stop("Could not find ./data/phenotypes/accessions_senescence.csv")
}

# Create output directory
dir.create("./data_processed/phenotypes/", recursive = TRUE)


#
# Load packages ----
#
library(tidyverse)

# Custom functions
source("scripts/R/functions/rdpiCalculator.R")


#
# Read data ----
#
acc_sil <- read_csv("./data/phenotypes/accessions_silique.csv", 
                    col_types = cols(id_leyser = col_character(),
                                     id_gwapp = col_character(),
                                     name = col_character(),
                                     nitrate = col_character(),
                                     height = col_double(),
                                     rosette = col_integer(),
                                     cauline = col_integer(),
                                     totalbr = col_integer(),
                                     nodes = col_integer(),
                                     bolt = col_integer()
                    ))

acc_sen <- read_csv("./data/phenotypes/accessions_senescence.csv", 
                    col_types = cols(id_leyser = col_character(),
                                     id_gwapp = col_character(),
                                     name = col_character(),
                                     totalbr = col_integer(),
                                     totalsil = col_integer()
                    ))


#
# Clean data ----
#
# log transform flowering time
acc_sil <- acc_sil %>%
  mutate(boltLog = log2(bolt))

# square-root transform silique number
acc_sen <- acc_sen %>%
  mutate(totalsilSqrt = sqrt(totalsil))


#
# Summary stats ----
#
# Calculate summary stats: mean, sd, median
acc_sum_sil <- acc_sil %>% 
  group_by(id_leyser, id_gwapp, name, nitrate) %>%
  summarise_each(list(n = ~sum(!is.na(.)),
                      mean = ~mean(., na.rm=T),
                      sd = ~sd(., na.rm=T),
                      cv = ~sd(., na.rm=T)/mean(., na.rm=T),
                      median = ~median(., na.rm=T))) %>% 
  ungroup()

acc_sum_sen <- acc_sen %>% 
  group_by(id_leyser, id_gwapp, name, nitrate) %>%
  summarise_each(list(n = ~sum(!is.na(.)),
                      mean = ~mean(., na.rm=T),
                      sd = ~sd(., na.rm=T),
                      cv = ~sd(., na.rm=T)/mean(., na.rm=T),
                      median = ~median(., na.rm=T))) %>% 
  ungroup()



#
# Plasticity indexes ----
#
# Calculate plasticity as difference of means
acc_plas_sil <- acc_sum_sil %>%
  group_by(id_leyser, id_gwapp, name) %>%
  summarise_each(list(hn = ~ .[nitrate == "HN"],
                      ln = ~ .[nitrate == "LN"],
                      D = ~ .[nitrate == "HN"] - .[nitrate == "LN"]),
                 matches("_mean")) %>%
  ungroup()

acc_plas_sen <- acc_sum_sen %>%
  group_by(id_leyser, id_gwapp, name) %>%
  summarise_each(list(hn = ~ .[nitrate == "HN"],
                      ln = ~ .[nitrate == "LN"],
                      D = ~ .[nitrate == "HN"] - .[nitrate == "LN"]),
                 matches("_mean")) %>%
  ungroup()


# Calculate relative distance plasticity index (RDPI) and add to table
rdpi_sil <- acc_sil %>% 
  group_by(id_leyser) %>% 
  summarise_each(list(rdpi = ~ rdpiCalculator(.[nitrate == "HN"], .[nitrate == "LN"])),
                 height:boltLog) %>% 
  ungroup()

rdpi_sen <- acc_sen %>% 
  group_by(id_leyser) %>% 
  summarise_each(list(rdpi = ~ rdpiCalculator(.[nitrate == "HN"], .[nitrate == "LN"])),
                 totalbr:totalsilSqrt) %>% 
  ungroup()

# Merge both tables
acc_plas_sil <- left_join(acc_plas_sil, rdpi_sil, by= "id_leyser")
acc_plas_sen <- left_join(acc_plas_sen, rdpi_sen, by= "id_leyser")


rm(rdpi_sil, rdpi_sen)


#
# Write data ----
#
# Writing files in RDS format so that column types are retained when read later on
acc_sil %>% 
  write_rds("data_processed/phenotypes/accessions_individual_silique.rds")

acc_sum_sil %>% 
  write_rds("data_processed/phenotypes/accessions_summarised_silique.rds")

acc_plas_sil %>% 
  write_rds("data_processed/phenotypes/accessions_plasticity_silique.rds")

acc_sen %>% 
  write_rds("data_processed/phenotypes/accessions_individual_senescence.rds")

acc_sum_sen %>% 
  write_rds("data_processed/phenotypes/accessions_summarised_senescence.rds")

acc_plas_sen %>% 
  write_rds("data_processed/phenotypes/accessions_plasticity_senescence.rds")



