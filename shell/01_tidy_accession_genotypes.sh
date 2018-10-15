#!/bin/bash
#$-N tidy_gen
#$-o $HOME/projects/20150501_accessions/scripts/shell/01_tidy_accession_genotypes.stout
#$-e $HOME/projects/20150501_accessions/scripts/shell/01_tidy_accession_genotypes.sterr
#$-pe smp 6
#$-l mem_free=30G
#$-S /bin/bash

# qsub ~/projects/20150501_accessions/scripts/shell/01_tidy_accession_genotypes.sh

cd ~/projects/20150501_accessions

# Genotypes downloaded from the link below, kindly provided by Umit Seren (GMI)
# Link will expire, so just here for reference
#  http://www.scidrive.org/vospace-2.0/data/3b7572cd-412a-42c6-b357-e35c4f373284



#### Convert genotypes to ped format ####

# 250k
/home/hugot/software/anaconda2/2.4.0/bin/python ./scripts/python/hdf5_to_ped_converter.py \
-i ./genotypes/PYGWAS_GENOTYPES/1/all_chromosomes_binary.hdf5 \
-o ./genotypes/snp250k/pygwas_genotypes

# Imputed
/home/hugot/software/anaconda2/2.4.0/bin/python ./scripts/python/hdf5_to_ped_converter.py \
-i ./genotypes/PYGWAS_GENOTYPES/5/all_chromosomes_binary.hdf5 \
-o ./genotypes/imputed/pygwas_genotypes



#### Convert ped files to binary .bed format ####

# 250k
plink2 --file ./genotypes/snp250k/pygwas_genotypes \
--out ./genotypes/snp250k/pygwas_genotypes \
--make-bed --allow-no-sex

# Imputed
plink2 --file ./genotypes/imputed/pygwas_genotypes \
--out ./genotypes/imputed/pygwas_genotypes \
--make-bed --allow-no-sex

# Remove .ped files, to save space
rm ./genotypes/imputed/pygwas_genotypes.ped 
rm ./genotypes/snp250k/pygwas_genotypes.ped


#### Make prunned .bed files #####

# 250k
plink2 --bfile ./genotypes/snp250k/pygwas_genotypes \
--out ./genotypes/snp250k/pygwas_genotypes_maf5 \
--maf 0.05 --make-bed --allow-no-sex

# Imputed
plink2 --bfile ./genotypes/imputed/pygwas_genotypes \
--out ./genotypes/imputed/pygwas_genotypes_maf5 \
--maf 0.05 --make-bed --allow-no-sex


#### Make relatedness matrix ####

# 250k
gcta64 --bfile ./genotypes/snp250k/pygwas_genotypes_maf5 \
--make-grm-bin \
--autosome-num 5 --autosome \
--out ./genotypes/snp250k/pygwas_genotypes_maf5 \
--thread-num 6

# Imputed
gcta64 --bfile ./genotypes/imputed/pygwas_genotypes_maf5 \
--make-grm-bin \
--autosome-num 5 --autosome \
--out ./genotypes/imputed/pygwas_genotypes_maf5 \
--thread-num 6


#### Convert files back to HDF5 format for LIMIX ####

# This seems ridiculous, but I think LIMIX is in a different format 
# then the original HDF5 files from Umem
/home/hugot/software/anaconda2/2.4.0/bin/limix_converter/limix_converter \
--plink=./genotypes/snp250k/pygwas_genotypes_maf5 \
--outfile=./genotypes/snp250k/pygwas_genotypes_maf5_limix.hdf5

# Also convert phenotype (this should not be on this script, but just 
# ran it straight on the command line)
/home/hugot/software/anaconda2/2.4.0/bin/limix_converter/limix_converter \
--csv=./phenotypes/limix/accession_phenotypes_silique.csv \
--outfile=./phenotypes/limix/accession_phenotypes_silique.hdf5

/home/hugot/software/anaconda2/2.4.0/bin/limix_converter/limix_converter \
--csv=./phenotypes/limix/accession_phenotypes_silique_early.csv \
--outfile=./phenotypes/limix/accession_phenotypes_silique_early.hdf5

/home/hugot/software/anaconda2/2.4.0/bin/limix_converter/limix_converter \
--csv=./phenotypes/limix/accession_phenotypes_senescence.csv \
--outfile=./phenotypes/limix/accession_phenotypes_senescence.hdf5
