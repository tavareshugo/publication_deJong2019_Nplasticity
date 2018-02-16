#!/bin/bash

#$-N plink_test
#$-o $HOME/projects/20150501_accessions/scripts/shell/03_gwas_tests.stout
#$-e $HOME/projects/20150501_accessions/scripts/shell/03_gwas_tests.sterr
#$-pe smp 6

## qsub ~/projects/20150501_accessions/scripts/shell/03_gwas_tests.sh

cd ~/projects/20150501_accessions/

bed='./genotypes/snp250k/pygwas_genotypes'
pheno='./phenotypes/plink/accession_phenotypes'
pca='./gwas/snp250k/basic_stats/all_accessions/all_accessions.eigenvec'



#### Create directories ####

mkdir -p ./gwas/snp250k/tests

outdir='./gwas/snp250k/tests'



#### Run several approaches ####

traits=`head -1 ${pheno}_silique.tsv | sed 's/FID\tIID\t//'`

ntraits=`echo $traits | wc -w`

for i in $(seq 1 $ntraits)
do
	trait=`echo $traits | cut -d " " -f $i`
	
	# MLMA
	gcta64 --bfile $bed \
	--out $outdir/${trait}_silique \
	--mlma \
	--grm-bin $bed \
	--pheno ${pheno}_silique.tsv \
	--mpheno $i \
	--autosome-num  5 \
	--maf 0.05 \
	--thread-num 6
	
	# MLMA-LOCO
	gcta64 --bfile $bed \
	--out $outdir/${trait}_silique \
	--mlma-loco \
	--grm-bin $bed \
	--pheno ${pheno}_silique.tsv \
	--mpheno $i \
	--autosome-num  5 \
	--maf 0.05 \
	--thread-num 6
	
	# MLMA + PCA
	gcta64 --bfile $bed \
	--out $outdir/${trait}_silique_pca \
	--mlma \
	--qcovar ./gwas/snp250k/basic_stats/all_accessions/all_accessions.eigenvec \
	--grm-bin $bed \
	--pheno ${pheno}_silique.tsv \
	--mpheno $i \
	--autosome-num  5 \
	--maf 0.05 \
	--thread-num 6
	
	# MLMA-LOCO + PCA
	gcta64 --bfile $bed \
	--out $outdir/${trait}_silique_pca \
	--mlma-loco \
	--qcovar ./gwas/snp250k/basic_stats/all_accessions/all_accessions.eigenvec \
	--grm-bin $bed \
	--pheno ${pheno}_silique.tsv \
	--mpheno $i \
	--autosome-num  5 \
	--maf 0.05 \
	--thread-num 6
	
	# PCA (plink)
	plink2 --bfile $bed \
	--out $outdir/${trait}_silique_pca \
	--allow-no-sex \
	--prune \
	--pheno ${pheno}_silique.tsv \
	--pheno-name $trait \
	--missing-phenotype -100 \
	--maf 0.05 \
	--linear \
	--covar ./gwas/snp250k/basic_stats/all_accessions/all_accessions.eigenvec

	# Simple linear model without correction (plink)
	plink2 --bfile $bed \
	--out $outdir/${trait}_silique \
	--allow-no-sex \
	--prune \
	--pheno ${pheno}_silique.tsv \
	--pheno-name $trait \
	--missing-phenotype -100 \
	--maf 0.05 \
	--linear
	
done


# With FT as co-variate
plink2 --bfile $bed \
--out $outdir/tss_mean_ln_covar \
--allow-no-sex \
--prune \
--pheno $pheno \
--covar $pheno \
--covar-name ft_mean_ln \
--pheno-name tss_mean_ln \
--missing-phenotype -100 \
--maf 0.05 \
--linear
