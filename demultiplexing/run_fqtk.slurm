#!/bin/bash

#SBATCH --partition=shortterm
#SBATCH --nodes=1
#SBATCH -c 8
#SBATCH --mem=128GB
#SBATCH --tmp=100G
#SBATCH --output=logs/%j_%u_%N_running_slurm_job.out
#SBATCH --error=logs/%j_%u_%N_running_slurm_job.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alihassan1697@gmail.com


# set up path for conda env
PATH=/work/hassan/hassan/miniforge:$PATH
source /work/hassan/hassan/miniforge/etc/profile.d/conda.sh

# Activate conda environment
conda activate fqtk

# Run FQTK
fqtk demux \
    --inputs combined_R1.fq.gz combined_R2.fq.gz \
    --read-structures 150T10B 150T10B \
    --sample-metadata metadata_GS369.tsv \
    --output output_folder