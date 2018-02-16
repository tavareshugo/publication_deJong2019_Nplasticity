#!/bin/python

import h5py
import numpy as np
from argparse import ArgumentParser
import os


#### Argument definition ####
def getArguments():
    # Usage message
    usage = """
      hdf5_to_ped_converter.py -i genotype/directory/genotypes.hdf5
                               -o plink/output/directory

    Details:
      This script converts SNP genotypes stored in hdf5 format to plink
      format. Output files are named the same as the input file (excluding 
      file extension).
      Note: currently this loads entire file to memory. For extremely large
      files this might be a problem. This might be solved in the future if
      data is read by chunks.

      Output plink files: 
          Two files are produced: .ped and .map
          The family ID and individual ID are the same.
          Sex and phenotype are set to missing.
          Alleles coded as 1 and 2.            
    """

    parser = ArgumentParser(description = "Convert .hdf5 genotypes to plink format", usage=usage)

    parser.add_argument("-i", "--input", required = True, help = "input hdf5 with genotypes.")
    parser.add_argument("-o", "--output", required = True, help = "output file name prefix (two files are created with .map and .ped extensions). If output directory is non-existent it will be created.")

    args = parser.parse_args()
    
    return(args)


#### Make plink .map ####
def makeMap(h5in, outfile):
	
	# Get SNP positions
	pos = h5py.File(h5in, 'r')['positions']
	
	# Print some verbose information
	print( "Creating .map file.\nThere are {} SNPs.\n".format(pos.size) )
	
	# Produce .map file
	with open(outfile + ".map", "w") as mapout:
		for chromosome in pos.attrs['chrs']:
			start = pos.attrs['chr_regions'][int(chromosome)-1][0]
			end = pos.attrs['chr_regions'][int(chromosome)-1][1]
			positions = pos[start:end]
			
			for position in positions:
				mapout.write( "{} {}_{} 0 {}\n".format(chromosome, chromosome, position, position))


#### Make plink .ped ####
def makePed(h5in, outfile):
	
	# Get genotypes and accession IDs
	snp = h5py.File(h5in, 'r')['snps']
	acc = h5py.File(h5in, 'r')['accessions']
	
	# Number of individuals and SNPs
	nind = snp.attrs['num_accessions']
	nsnps = snp.attrs['num_snps']
	
	# Print some verbose information
	print( "Creating .ped file.\nThere are {} SNPs and {} accessions.\n".format(nsnps, nind) )
	
	# Transpose SNP array so individuals are rows and genotypes are columns
	snp_tr = snp[0:nsnps].T
	
	del(snp) # to save some memory
	
	# Repeat each element of array twice (because plink expects two alleles for every SNP)
	snp_tr = np.repeat(snp_tr, 2, axis = 1)
	
	with open(outfile + ".ped", "w") as pedout:
		for i in range(0, nind):
				pedout.write( "{} {} 0 0 0 -9 {}\n".format(acc[i], acc[i], " ".join(map(str, snp_tr[i,:]+1))) )
		

#### Main ####
if __name__ == "__main__":
	# Get user arguments
	args = getArguments()
	
	# Create output directory
	if not os.path.exists(os.path.dirname(args.output)):
		os.makedirs(os.path.dirname(args.output))
	
	# Map file:
	makeMap(args.input, args.output)
	
	# Ped file:
	makePed(args.input, args.output)
	
