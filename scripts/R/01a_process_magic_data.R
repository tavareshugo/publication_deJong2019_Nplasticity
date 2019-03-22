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
if(!file.exists("./data/phenotypes/magic_silique.csv")){
  stop("Could not find ./data/phenotypes/magic_silique.csv")
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
magic <- read_csv("./data/phenotypes/magic_silique.csv", 
                  col_types = cols(id_leyser = col_character(),
                                   id_kover = col_character(),
                                   nitrate = col_character(),
                                   height = col_double(),
                                   rosette = col_integer(),
                                   cauline = col_integer(),
                                   totalbr = col_integer(),
                                   nodes = col_integer(),
                                   bolt = col_integer()
                  ))


#
# Clean data ----
#
# Store this version of the data in a separate variable, as it contains all the data
magic_raw <- magic

# Retain only lines with at least 4 replicates in each nitrate treatment
magic <- magic %>% 
  group_by(id_kover, nitrate) %>% 
  filter(n() >= 4)

# Retain only lines with data for both nitrate treatments
magic <- magic %>% 
  group_by(id_kover) %>% 
  filter(n_distinct(nitrate) == 2) %>% 
  ungroup()

# log-transform flowering time
magic <- magic %>%
  mutate(boltLog = log2(bolt))


#
# Summary stats ----
#
# Calculate summary stats: mean, sd, median
magic_sum <- magic %>% 
  group_by(id_leyser, id_kover, nitrate) %>%
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
magic_plas <- magic_sum %>% 
  group_by(id_leyser, id_kover) %>% 
  summarise_each(list(hn = ~ .[nitrate == "HN"],
                      ln = ~ .[nitrate == "LN"],
                      D = ~ .[nitrate == "HN"] - .[nitrate == "LN"]), 
                 matches("_mean")) %>% 
  ungroup()

# Calculate RDPI (from individual data)
rdpi <- magic %>% 
  group_by(id_leyser, id_kover) %>% 
  summarise_each(list(rdpi = ~ rdpiCalculator(.[nitrate == "HN"], .[nitrate == "LN"])),
                 height:boltLog) %>% 
  ungroup()

# Merge both tables
magic_plas <- left_join(magic_plas, rdpi, by= c("id_leyser", "id_kover"))

rm(rdpi)


#
# Write data ----
#
# Writing files in RDS format so that column types are retained when read later on
magic %>% 
  write_rds("data_processed/phenotypes/magic_individual.rds")

magic_sum %>% 
  write_rds("data_processed/phenotypes/magic_summarised.rds")

magic_plas %>% 
  write_rds("data_processed/phenotypes/magic_plasticity.rds")



