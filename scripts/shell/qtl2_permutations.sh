#!/bin/bash
#SBATCH --job-name=qtl2
#SBATCH -A LEYSER-SL3-CPU
#SBATCH --workdir=/home/hm533/scratch/projects/2018_deJong_NPlas/supplementary_data
#SBATCH --array=1-1000
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=00:20:00
#SBATCH -p skylake  #for 6G per CPU use "-p skylake"
#SBATCH -o ./scripts/shell/qtl2_permutations.stdout

# sbatch /home/hm533/scratch/projects/2018_deJong_NPlas/supplementary_data/scripts/shell/qtl2_permutations.sh

# Launch the R script setting the seed based on job array number
Rscript --vanilla ./scripts/R/run_qtl2_lmm_perm.R --workdir "$(pwd)" --seed "1000$SLURM_ARRAY_TASK_ID"

wait
