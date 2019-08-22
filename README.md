# Analysis code for _de Jong et al. 2019_

This repository contains code used in the analysis of data from [de Jong et al. 2019](). 

All of the code (and data) is deposited as a static snapshot in the Cambridge Apollo repository 
(https://doi.org/10.17863/CAM.20379), corresponding to release 0.1 of the 
[GitHub repository](https://github.com/tavareshugo/publication_deJong2019_Nplasticity). 

The code in the GitHub repository is non-static and may be updated if we find any 
bugs (see "NEWS" section at the bottom).


#### Setup

1. Clone the repository (or [download as a zip file](https://github.com/tavareshugo/publication_deJong2019_Nplasticity/archive/master.zip)
and extract the files) - you should then have a project directory called _"publication_deJong2019_Nplasticity"_
2. Download and extract the zip files with the [raw data]() and the [processed data]() 
from the Appolo repository. Place them inside the previous project directory. 

In the end you should have this directory structure: 

```
publication_deJong2019_Nplasticity
          ├── data
          ├── data_processed
          └── scripts
```


#### Running analysis

Most of the analysis was done using `R`. 

To recreate the analysis and figures for the paper, open the `2019_deJong_Nplas.Rproj` 
RStudio project file (which will ensure you have the correct working directory 
set - all R scripts use relative paths from the project directory).

There is also a document detailing [how the mixed models were fit to trait data](./docs/01_variance_partitioning.md) 

Some details about the files on `scripts/R`:

- All scripts named with prefix `Fig`/`Table` recreate each figure/table in the paper. 
Some of the graphs were slightly edited in inkscape for the published version.
- The scripts with numbered prefix `01`-`04` do different data processing steps and produce files already provided in the `data_processed` folder in our repository. However, they can be used to recreate those analysis. 
  - `01` - tidy the phenotype data
  - `02` - fit linear mixed models on trait data (used to estimate variance components and heritabilities)
  - `03` - run standard QTL scans using `R/qtl2` package
  - `04` - run custom QTL scans using a mixed model approach, as well as a permutation for significance thresholds.

The GWAS was done using standalone command-line programs. 
The code is in the `scripts/shell` directory. 
Those scripts are numbered according to the order in which analysis was run. 
This part of the analysis is not completely reproducible in its current form, as we used a pre-release version of the Arabidopsis genotype data. However, a similar analysis should be possible to run using the [publicly available SNP data](https://github.com/Gregor-Mendel-Institute/atpolydb/wiki).

Finally, we used the python package `LIMIX` to fit a multi-SNP QTL model (`scripts/python/limix_multitrait.py`) as published by [Sasaki et al 2015](https://doi.org/10.1371/journal.pgen.1005597).


#### Software versions

```
# R and packages
R 3.5.3
tidyverse 1.2.1
patchwork 0.0.1
broom 0.5.2
qtl2 0.18
atMAGIC 0.1.0

# python and LIMIX
python 2.7.12
limix 0.7.12
```

