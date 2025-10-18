#!/bin/bash

#SBATCH --partition=shortterm
#SBATCH --nodes=1
#SBATCH -c 8
#SBATCH --mem=128GB
#SBATCH --tmp=100G
#SBATCH --job-name=running_cadd
#SBATCH --output=slurm_logs/%j_%u_%N_running_cadd.out
#SBATCH --error=slurm_logs/%j_%u_%N_running_cadd.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alihassan1697@gmail.com

## Usage check
## Real Test Command for just running CADD without slurm:
## ./CADD.sh -g GRCh38 -o /data/humangen_kircherlab/Users/hassan/repos/sandbox_launchpad/final_outputs/cadd_annotations test/input.vcf
## sbatch --job-name=anno_cadd_svs cadd_test.bash  /data/humangen_sfb1665_seqdata/short_read/processed_data/batch_001_size10/A4842_DNA_01/call_sv/genome/case_A4842_DNA_01_sv.vcf.gz  /data/humangen_kircherlab/Users/hassan/repos/sandbox_launchpad/final_outputs/cadd_annotations

# Check if input file and output directory are provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 <input_vcf_file> <output_directory>"
    echo "Example: $0 input.vcf.gz /path/to/output"
    exit 1
fi

INPUT_VCF="$1"
OUTPUT_DIR="$2"

# Check if input file exists
if [ ! -f "${INPUT_VCF}" ]; then
    echo "Error: Input file ${INPUT_VCF} does not exist"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

# Extract base filename without path and extension
BASENAME=$(basename "${INPUT_VCF}" .vcf.gz)

# Define intermediate and final output file paths
TEMP_VCF="${OUTPUT_DIR}/${BASENAME}_nochr.vcf.gz"
OUTPUT_CADD="${OUTPUT_DIR}/${BASENAME}_CADD_GRCh38.tsv.gz"

PATH=/work/hassan/hassan/miniforge:$PATH
source /work/hassan/hassan/miniforge/etc/profile.d/conda.sh

module load singularity/v4.1.3

# Remove chr prefix and recompress
echo "Removing 'chr' prefix from VCF..."
zcat ${INPUT_VCF} | sed 's/^chr//' | bgzip > ${TEMP_VCF}

# Index the new VCF
echo "Indexing the new VCF..."
tabix -p vcf ${TEMP_VCF}

echo "activating conda environment..."
conda activate cadd

cd /data/humangen_kircherlab/Users/hassan/chaos_lab/cadd

# Run CADD with specific input and output files
echo "Running CADD..."
./CADD.sh -g GRCh38 -o ${OUTPUT_CADD} ${TEMP_VCF}

rm ${TEMP_VCF} ${TEMP_VCF}.tbi

echo "CADD analysis complete. Output file: ${OUTPUT_CADD}"