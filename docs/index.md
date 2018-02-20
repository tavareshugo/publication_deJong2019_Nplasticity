
[This repository](https://github.com/tavareshugo/publication_deJong2018_Nplasticity) 
contains several scripts with analysis code for de Jong et al (2018).

A citable (static) version is also available from the Cambridge University data 
repository Appolo (DOI: [10.17863/CAM.20379](https://doi.org/10.17863/CAM.20379)).

There are three main directories:

### R

This directory contains all the R code to do the analysis and plots shown in the 
paper. 

All the data used by the scripts is available for download from the Cambridge 
University data repository Appolo (DOI: [10.17863/CAM.20379](https://doi.org/10.17863/CAM.20379)). 
A description of the data is provided with the data in a README file.

The directory structure necessary for the scripts to work is:

```
project_folder
 |
 |_ data - unzip the data file from the Appolo repository
 |
 |_ scripts - download/clone this repository to a separate directory
```

These analysis scripts were produced using Rmarkdown and are rendered as individual 
pages, with detailed information and code:

* [Partitioning of trait variances](01_variance_partitioning.html)
* [Analysis of trait variation](02_analysis_of_variation.html)
* [Analysis of N-pulse experiments](03_pulse_analysis.html)
* [Analysis of qPCR and grafting experiments](04_qpcr_graft_analysis.html)
* [Estimating broad-sense heritabilities](05_estimating_heritabilities.html)
* [QTL mapping in MAGIC lines](06_qtl_mapping.html)
* [QTL mapping in accessions](07_accession_gwas.html)


### python


This directory contains a single python script that was writted to convert 
genotypes from HDF5 format to plink format. It was run from the shell script 
`01_tidy_accession_genotypes.sh` (see below).

### shell

This directory contains several shell scripts used to perform GWAS in the 
natural accessions. 

The scripts are numbered according to the order in which they were run. 

These scripts use specific paths of the filesystem used to run them, so would 
need to be adjusted to reproduce our analysis. 

Summarily, these scripts do the following:

* `01_tidy_accession_genotypes` - convert genotypes to plink format; calculate 
genetic relatedness matrix
* `02_basic_stats_250k` - calculate some basic descriptive statistics such as 
allele frequencies and PCA based on the 250K SNP data
* `02_basic_stats_imputed` - same as above but for the imputed dataset
* `03_gwas_tests` - run the GWAS tests using different models. This was used to 
test different approaches. In the paper we report the results for the 
`--mlma-loco` method in GCTA.
* `04_gwas_gcta_250k` - GWAS using the 250K SNP set. It does the GWAS scan but also 
calculates GWAS-heritability and does the gene set-based test implemented in GCTA.
* `04_gwas_gcta_imputed` - same as above but using the imputed SNP set.
