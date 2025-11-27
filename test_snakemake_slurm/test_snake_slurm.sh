#!/bin/bash

#SBATCH --partition=shortterm
#SBATCH --nodes=1
#SBATCH -c 16
#SBATCH --mem=128GB
#SBATCH --tmp=100G
#SBATCH --job-name=test_snake_slurm
#SBATCH --output=logs/%j_%u_%N_test_snake_slurm.out
#SBATCH --error=logs/%j_%u_%N_test_snake_slurm.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alihassan1697@gmail.com



# # ## Load Snakemake module
# module load snakemake/9.13.7


# set up path for conda env
PATH=/work/hassan/hassan/miniforge:$PATH
source /work/hassan/hassan/miniforge/etc/profile.d/conda.sh

# Activate conda environment
conda activate snakemake

echo "Snakemake version:"
snakemake --version


# Run Snakemake workflow with slurm profile
snakemake --profile slurm_profile --cores 5