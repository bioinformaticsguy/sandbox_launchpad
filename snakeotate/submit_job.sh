#!/bin/bash

#SBATCH --partition=shortterm
#SBATCH --nodes=1
#SBATCH -c 16
#SBATCH --mem=128GB
#SBATCH --tmp=100G
#SBATCH --job-name=snakotate
#SBATCH --output=logs/%j_%u_%N_snakotate.out
#SBATCH --error=logs/%j_%u_%N_snakotate.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alihassan1697@gmail.com

## Usage check


# ## Load Snakemake module
# module load snakemake/9.13.7

# set up path for conda env
PATH=/work/hassan/hassan/miniforge:$PATH
source /work/hassan/hassan/miniforge/etc/profile.d/conda.sh

# Activate conda environment
conda activate snakemake


echo "Snakemake version:"
snakemake --version

# Print working directory
echo "Working directory: $(pwd)"

## Run Snakemake with small test data config (Works)
# snakemake -c 16 --use-conda --snakefile snakeotate/workflow/Snakefile --configfile snakeotate/resources/snv_annotation/config.yaml

# # Run Snakemake workflow with real data config
# snakemake -c 32 --use-conda --snakefile snakeotate/workflow/Snakefile --configfile snakeotate/resources/snv_annotation/real_config.yaml

# # # Run Snakemake workflow with slurm profile
snakemake -c 16 --use-conda --snakefile snakeotate/workflow/Snakefile --configfile snakeotate/resources/snv_annotation/config.yaml --profile snakeotate/myprofile

# snakemake --profile snakeotate/myprofile --configfile snakeotate/resources/snv_annotation/config.yaml

# # # Run Snakemake workflow
# snakemake -c 16 --use-conda --snakefile snakeotate/workflow/Snakefile --configfile snakeotate/resources/snv_annotation/config.yaml --executor slurm  
#  slurm_partition=shortterm --jobs 2