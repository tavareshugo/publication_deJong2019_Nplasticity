# Analysis code for publication de Jong et al. 2019

# Analysis code for _de Jong et al. 2019_

**TODO**

* add link to appolo
* add link to paper
* add more concrete explanation for downloading data and having right directory structure

This repository contains code used in the analysis of data from [de Jong et al. 2019](). 

All of the code (and data) is deposited as a static snapshot in the Appolo repository 
(**DOI link**). This repository is non-static and may be updated if we find any 
bugs (a "NEWS" section will be created in this README if that happens.)


#### Setup

1. Clone this repository (or [download as a zip file](https://github.com/tavareshugo/publication_deJong2019_Nplasticity/archive/master.zip)
and extract the files)
2. Download and extract the zip files with the [raw data]() and the [processed data]() 
from the Appolo repository

In the end you should have this directory structure: 

```
publication_deJong2019_Nplasticity
          ├── data
          ├── data_processed
          └── scripts
```

To recreate analysis and figures for the paper, open the `2019_deJong_Nplas.Rproj` 
RStudio project file (which will ensure you have the correct working directory 
set - all R scripts use relative paths from the project directory).


There is also a document detailing [how the mixed models were fit to trait data](./docs/01_variance_partitioning.md) 
