#!/bin/bash

#$-N gcta_250k
#$-o $HOME/projects/20150501_accessions/scripts/shell/04_gwas_gcta_250k.stout
#$-e $HOME/projects/20150501_accessions/scripts/shell/04_gwas_gcta_250k.sterr
#$-pe smp 6
#$-hold_jid basic_250k

## qsub ~/projects/20150501_accessions/scripts/shell/04_gwas_gcta_250k.sh

cd ~/projects/20150501_accessions/

#### Create output directory ####

mkdir -p ./gwas/snp250k/

cd ./gwas/snp250k


#### Important files ####
# Genotype file prefix
bed='../../genotypes/snp250k/pygwas_genotypes_maf5'

# Phenotype file prefix
pheno='../../phenotypes/plink/accession_phenotypes'

# Gene annotation file
genes='/home/hugot/reference/arabidopsis/tair10/annotation/Arabidopsis_thaliana.TAIR10.27.tsv'


#### Run association tests ####

## Two different methods both with relatedness correction
## For each I apply a set-based method (centered on genes)

for set in silique silique_early senescence senescence_early
do
	# Make output directories
	mkdir -p ./gwas_heritability/${set}
	mkdir -p ./mlm_gcta/${set}
	mkdir -p ./mlm-loco_gcta/${set}
	
	# List of traits
	traits=$(head -1 ${pheno}_${set}.tsv | sed 's/FID\tIID\t//')
	
	# Number of traits
	ntraits=$(echo $traits | wc -w)

	# For each trait
	for i in $(seq 1 $ntraits)
	do
		# Get the trait
		trait=$(echo $traits | cut -d " " -f $i)
		
		# Estimate GWAS heritability
		gcta64 --reml \
		--out ./gwas_heritability/${set}/${trait} \
		--grm-bin $bed \
		--pheno ${pheno}_${set}.tsv \
		--mpheno $i \
		--autosome-num  5 \
		--maf 0.05 \
		--thread-num 6
		
		
		### MLM ###
		# Association test
		gcta64 --bfile $bed \
		--out ./mlm_gcta/${set}/${trait} \
		--mlma \
		--grm-bin $bed \
		--pheno ${pheno}_${set}.tsv \
		--mpheno $i \
		--autosome-num  5 \
		--maf 0.05 \
		--thread-num 6
		
		# Produce a file of p-values for set-based GWAS
		cat ./mlm_gcta/${set}/${trait}.mlma | \
		awk '{print $2 " " $9}' > ./mlm_gcta/${set}/${trait}.pvals
		
		# Run set-based GWAS
		gcta64 \
			--bfile $bed \
			--autosome-num  5 \
			--maf 0.05 \
			--thread-num 6 \
			--fastBAT ./mlm_gcta/${set}/${trait}.pvals \
			--fastBAT-gene-list ${genes} \
			--fastBAT-wind 5 \
			--fastBAT-ld-cutoff 0.9 \
			--out ./mlm_gcta/${set}/${trait}
		
		rm ./mlm_gcta/${set}/${trait}.pvals
		
		
		#### MLMA LOCO ####
		# Run association test
		gcta64 --bfile $bed \
		--out ./mlm-loco_gcta/${set}/${trait} \
		--mlma-loco \
		--grm-bin $bed \
		--pheno ${pheno}_${set}.tsv \
		--mpheno $i \
		--autosome-num  5 \
		--maf 0.05 \
		--thread-num 6
		
		# Produce a file of p-values for set-based GWAS
		cat ./mlm-loco_gcta/${set}/${trait}.loco.mlma | \
		awk '{print $2 " " $9}' > ./mlm-loco_gcta/${set}/${trait}.pvals
		
		# Run set-based GWAS
		gcta64 \
			--bfile $bed \
			--autosome-num  5 \
			--maf 0.05 \
			--thread-num 6 \
			--fastBAT ./mlm-loco_gcta/${set}/${trait}.pvals \
			--fastBAT-gene-list $genes \
			--fastBAT-wind 5 \
			--fastBAT-ld-cutoff 0.9 \
			--out ./mlm-loco_gcta/${set}/${trait}
		
		rm ./mlm-loco_gcta/${set}/${trait}.pvals
		
				
		#~ # Clump significant SNPs from MLM approach
		#~ plink2 --bfile $bed \
		#~ --clump ${outdir}/${set}/${trait}.mlma \
		#~ --out ${outdir}/${set}/${trait}.mlma \
		#~ --clump-snp-field SNP \
		#~ --clump-field p \
		#~ --clump-p1 0.0001 \
		#~ --clump-kb 500 \
		#~ --clump-r2 0.2 \
		#~ --clump-range ~/reference/arabidopsis/tair10/annotation/Arabidopsis_thaliana.TAIR10.27.tsv \
		#~ --clump-range-border 10
		
	done
	
done


#### Extract heritability estimates ####
for set in silique silique_early senescence senescence_early
do
	# Go to respective directory
	cd ./gwas_heritability/${set}

	# Print file header
	printf "trait\tgenetic_var\tgenetic_se\tresidual_var\tresidual_se\tvp_var\tvp_se\ther_var\ther_se\tpval\n" > snp_heritabilities.tsv

	files=$(ls *.hsq)

	for file in $files
	do
		genvar=$(cat $file | grep "V(G)" | head -1 | cut -f 2,3)
		resvar=$(cat $file | grep "V(e)" | head -1 | cut -f 2,3)
		vp=$(cat $file | grep "Vp" | head -1 | cut -f 2,3)
		genvp=$(cat $file | grep "V(G)/Vp" | cut -f 2,3)
		pval=$(cat $file | grep "Pval" | cut -f 2)
		
		printf "$file\t$genvar\t$resvar\t$vp\t$genvp\t$pval\n" >> snp_heritabilities.tsv
	done
	
	# Get back to starting directory
	cd ../../

done


