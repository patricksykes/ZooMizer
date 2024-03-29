#!/bin/bash
#
#SBATCH --array=0-651
#SBATCH --job-name=20220611zoomizergrid_Rlevel0
#SBATCH --output=slurm_%a.out
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=5G
#SBATCH --time=0-05:00

module load R/4.0.3
R CMD BATCH zoomizerslurm.R