#!/bin/bash

#$-N basic_imputed
#$-o $HOME/projects/20150501_accessions/scripts/shell/02_basic_stats_imputed.stout
#$-e $HOME/projects/20150501_accessions/scripts/shell/02_basic_stats_imputed.sterr
#$-pe smp 6
#$-S /bin/bash
#$-hold_jid tidy_gen

## qsub ~/projects/20150501_accessions/scripts/shell/02_basic_stats_imputed.sh

cd ~/projects/20150501_accessions/

bed=./genotypes/imputed/pygwas_genotypes_maf5  #bed files
pheno=./phenotypes/plink/accession_phenotypes  #phenotype files



#### Make directories ####
# Plink outputs
mkdir -p ./gwas/imputed/basic_stats/all_accessions
mkdir -p ./gwas/imputed/basic_stats/phenotyped_silique_stage
mkdir -p ./gwas/imputed/basic_stats/phenotyped_senescence_stage



#### Basic Plink summary stats #####
# Allele frequency; missing data reports; pruned list SNPs based on LD;

outdir=./gwas/imputed/basic_stats  # Variable with output directory

# Get the allele frequency across all individuals
plink2 --bfile $bed \
--out $outdir/all_accessions/all_accessions \
--allow-no-sex \
--freq \
--missing

# Make a prunned list of SNPs based on their LD properties
plink2 --bfile $bed \
--out $outdir/all_accessions/all_accessions \
--allow-no-sex \
--maf 0.05 \
--indep-pairwise 50 5 0.5

# Phenotyped individuals
for set in silique senescence
do
	plink2 --bfile $bed \
	--out $outdir/phenotyped_${set}_stage/phenotyped_${set} \
	--allow-no-sex \
	--freq \
	--missing

	plink2 --bfile $bed \
	--out $outdir/phenotyped_${set}_stage/phenotyped_${set} \
	--maf 0.05 \
	--allow-no-sex \
	--indep-pairwise 50 5 0.5 \
	--blocks \
	--blocks-max-kb 500 \
	--blocks-min-maf 0.05 \
	--prune \
	--pheno ${pheno}_${set}.tsv \
	--pheno-name tss_mean_D \
	--missing-phenotype -100
done



#### Structure ####
# PCA; Clustering by IBS; MDS based on IBS
plink2 --bfile $bed \
--out $outdir/all_accessions/all_accessions \
--allow-no-sex \
--cluster \
--pca 10 header \
--mds-plot 10 \
--extract $outdir/all_accessions/all_accessions.prune.in





#~ #### Structure analysis ####
#~ mkdir -p ./gwas/imputed/admixture
#~ 
#~ cd ./gwas/imputed/admixture
#~ 
#~ # Run for 1 to 20 clusters
#~ for K in {1..20}
#~ do
  #~ ~/software/admixture/1.23/admixture -j6 \
  #~ --cv ~/projects/20150501_accessions/genotypes/imputed/pygwas_genotypes_maf5.bed $K > all.admix.${K}.log
#~ done
#~ 
#~ for K in 25 30 35
#~ do
  #~ ~/software/admixture/1.23/admixture -j6 \
  #~ --cv ~/projects/20150501_accessions/genotypes/imputed/pygwas_genotypes_maf5.bed $K > all.admix.${K}.log
#~ done
#~ 
#~ for K in 40 45
#~ do
  #~ ~/software/admixture/1.23/admixture -j6 \
  #~ --cv ~/projects/20150501_accessions/genotypes/imputed/pygwas_genotypes_maf5.bed $K > all.admix.${K}.log
#~ done
#~ 
#~ 
#~ # Extract cross-validation values to a nicely formated file
#~ printf "k cv_error\n" > all.admix.cross-validation_k1-20.txt
#~ 
#~ grep -h CV all.admix.*.log  | \
#~ cut -d " " -f3,4 | \
#~ sed 's/(K=//g' | \
#~ sed 's/)://g' | \
#~ sort -n >> all.admix.cross-validation_k1-20.txt



