# Import the (many) libraries
# All of the below could have been achieved with
# from setup import *
## I've just put it here to be more self-contained

import sys
sys.path.append('./..')
#import data as tutorial_data


import scipy as sp
import pylab as pl
from matplotlib import cm
import scipy.stats as st
import h5py
import pdb
import pandas as pd
sp.random.seed(0)
# import LIMIX
import sys
import limix.modules.varianceDecomposition as var
import limix.modules.qtl as qtl
import limix.io.data as data
import limix.io.genotype_reader as gr
import limix.io.phenotype_reader as phr
import limix.io.data_util as data_util
import limix.utils.preprocess as preprocess
# plotting and visualization utilties
from limix.utils.plot import *
# genotype summary stats
from limix.deprecated.stats.geno_summary import *
import os
import cPickle
import sys
import numpy as np
import pandas as pd


#--------------------#
#### Prepare data ####
#--------------------#

# Reader instance for genotypes
geno_reader  = gr.genotype_reader_tables('/home/hugot/projects/20150501_accessions/genotypes/snp250k/pygwas_genotypes_limix.hdf5')

# Reader instance for phenotypes
pheno_reader = phr.pheno_reader_tables('/home/hugot/projects/20150501_accessions/phenotypes/limix/accession_phenotypes_silique_early.hdf5')

# Combine genotypes and phenotypes into limix-specific object
dataset = data.QTLData(geno_reader = geno_reader, pheno_reader = pheno_reader)

# Get SNPs, phenotypes and positions in respective variables
snps = dataset.getGenotypes()

phenotypes = dataset.getPhenotypes(intersection = True)[0]

pos = dataset.getPos()
pos, chromBounds = data_util.estCumPos(position = pos, offset = 0)

# Subset only TSS trait for multi-trait LMM
phenotypes_tss = phenotypes[['totalbr_mean_ln', 'totalbr_mean_hn']]


# Estimate relatedness matrix
sample_relatedness = dataset.getCovariance(normalize = True, center = True)

# This extra step is in the tutorial. 
## However, the denominator is 1, so this effectively has no effect here
sample_relatedness = sample_relatedness / sample_relatedness.diagonal().mean()


#-----------------------#
#### Any effect test ####
#-----------------------#

# Number of samples (N) and phenotypes (P)
N, P = phenotypes_tss.shape

# Set all parameters of model
covs = None                 #covariates
Acovs = None                #the design matrix for the covariates   
Asnps = sp.eye(P)           #the design matrix for the SNPs
K1r = sample_relatedness    #the first sample-sample covariance matrix (non-noise)
K2r = sp.eye(N)             #the second sample-sample covariance matrix (noise)
K1c = None                  #the first phenotype-phenotype covariance matrix (non-noise)
K2c = None                  #the second phenotype-phenotype covariance matrix (noise)
covar_type = 'freeform'     #the type of the trait/trait covariance to be estimated 
searchDelta = False         #specify if delta should be optimized for each SNP
test="lrt"                  #specify type of statistical test

# Running the analysis
mtlmm_any, pvalues = qtl.test_lmm_kronecker(snps=snps,phenos=phenotypes_tss.values,covs=covs,Acovs=Acovs, Asnps=Asnps,K1r=K1r,trait_covar_type=covar_type)

# Format data in nice DataFrame
pvalues_any = pd.DataFrame(data = pvalues.T, index = dataset.geno_ID, columns = ['mtlmm_any'])


#--------------------------#
#### Common effect test ####
#--------------------------#

# Set new design matrix for the SNPs
Asnps = sp.ones((1,P))

# Run the analysis
mtlmm_common, pvalues = qtl.test_lmm_kronecker(snps = snps, phenos = phenotypes_tss.values, covs = covs, Acovs = Acovs, Asnps = Asnps, K1r = K1r, trait_covar_type = covar_type)

# Format the data nicely
pvalues_common = pd.DataFrame(data = pvalues.T, index=dataset.geno_ID, columns = ['mtlmm_common'])


#-----------------------#
#### GxE effect test ####
#-----------------------#

# Set new parameters for the model
#the null model design matrix for the SNPs
Asnps0 = sp.ones((1,P))

#the alternative model design matrix for the SNPs
Asnps1 = sp.zeros((2,P))
Asnps1[0,:] = 1.0
Asnps1[1,0] = 1.0

# Run the analysis
mtlmm_inter = qtl.test_interaction_lmm_kronecker(snps=snps, phenos=phenotypes_tss.values, covs=covs, Acovs=Acovs, Asnps0=Asnps0, Asnps1=Asnps1, K1r=K1r, trait_covar_type=covar_type)

# Tidy into dataframe. Also concatenate with previous tests
pvalues_inter = pd.DataFrame(data = sp.concatenate(mtlmm_inter).T, index = dataset.geno_ID, columns = ['specific', 'common', 'any'])


#-------------------#
#### Tidy output ####
#-------------------#

# Function to calculate minor allele frequencies
def afCalc(M):
    hom_minor = (M==0).sum(axis = 0)
    het = (M==1).sum(axis = 0)
    #hom_major = (snps==2).sum(axis = 0)
    
    maf = (2*hom_minor + het)/sp.double(2*M.shape[0])
    
    return(maf) 

# Get MAF to a dataframe
af_df = pd.DataFrame(data = afCalc(geno_reader.getGenotypes()), columns = ['maf'])

# Concatenate with test results and SNP positions
all_tests = pd.concat([pos, pvalues_inter, af_df], axis = 1)

all_tests.to_csv('/home/hugot/projects/20150501_accessions/gwas/snp250k/limix/multi_trait_silique_early.csv')
