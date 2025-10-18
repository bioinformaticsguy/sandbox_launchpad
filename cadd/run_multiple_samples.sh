#!/bin/bash

# Check if CSV file and output directory are provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 <csv_file> <output_directory>"
    echo "Example: $0 input_files.csv /path/to/output"
    echo ""
    echo "CSV format: one VCF file path per line"
    exit 1
fi

CSV_FILE="$1"
OUTPUT_DIR="$2"

# Check if CSV file exists
if [ ! -f "${CSV_FILE}" ]; then
    echo "Error: CSV file ${CSV_FILE} does not exist"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

# Create slurm_logs directory if it doesn't exist
mkdir -p slurm_logs

# Counter for jobs
JOB_COUNT=0

# Read CSV file line by line
while IFS= read -r VCF_PATH || [ -n "$VCF_PATH" ]; do
    # Skip empty lines
    if [ -z "$VCF_PATH" ]; then
        continue
    fi
    
    # Remove leading/trailing whitespace
    VCF_PATH=$(echo "$VCF_PATH" | xargs)
    
    # Check if VCF file exists
    if [ ! -f "$VCF_PATH" ]; then
        echo "Warning: File $VCF_PATH does not exist. Skipping..."
        continue
    fi
    
    # Extract sample name from VCF path
    SAMPLE_NAME=$(basename "$VCF_PATH" .vcf.gz)
    
    # Submit job
    echo "Submitting job for: $SAMPLE_NAME"
    sbatch --job-name="cadd_${SAMPLE_NAME}" cadd_test.bash "$VCF_PATH" "$OUTPUT_DIR"
    
    JOB_COUNT=$((JOB_COUNT + 1))
    
    # Optional: Add a small delay to avoid overwhelming the scheduler
    sleep 0.5
    
done < "$CSV_FILE"

echo ""
echo "Submitted $JOB_COUNT jobs to the queue"
echo "Monitor with: squeue -u $USER"