#!/bin/bash

# Configuration
OUTPUT_BASE_DIR="/data/humangen_kircherlab/Users/hassan/chaos_lab/final_outputs/combind_annotations"
SCRIPT_PATH="/data/humangen_kircherlab/Users/hassan/chaos_lab/snv_annotations/updated.bash"

# Default annotation databases (gnomAD and ClinVar are same for all)
GNOMAD_FILE="/data/humangen_kircherlab/Users/hassan/data/databases/genomAD/gnomad.all.AC_AN_AF.vcf.bgz"
CLINVAR_FILE="/data/humangen_kircherlab/Users/hassan/data/databases/clinvar/clinvar.vcf.gz"
VEP_DIR="ensamble_vep/vep_data/"
VEP_SIF="ensamble_vep/vep.sif"

# Function to display usage
usage() {
    cat << EOF
Usage: $0 -c <csv_file> [OPTIONS]

Required arguments:
  -c, --csv         CSV file with sample information

Optional arguments:
  -o, --output      Base output directory (default: $OUTPUT_BASE_DIR)
  -s, --script      Path to annotation script (default: $SCRIPT_PATH)
  -h, --help        Display this help message

CSV Format:
  sample_id,sample_varients,cadd_annotation_path

Example:
  $0 -c sample_list.csv
  $0 -c sample_list.csv -o /path/to/results
EOF
    exit 1
}

# Initialize variables
CSV_FILE=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -c|--csv)
            CSV_FILE="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_BASE_DIR="$2"
            shift 2
            ;;
        -s|--script)
            SCRIPT_PATH="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Check required arguments
if [[ -z "$CSV_FILE" ]]; then
    echo "Error: CSV file is required."
    usage
fi

# Check if CSV file exists
if [[ ! -f "$CSV_FILE" ]]; then
    echo "Error: CSV file '$CSV_FILE' not found."
    exit 1
fi

# Check if script exists
if [[ ! -f "$SCRIPT_PATH" ]]; then
    echo "Error: Annotation script '$SCRIPT_PATH' not found."
    exit 1
fi

# Create base output directory
mkdir -p "$OUTPUT_BASE_DIR"

# Create logs directory
mkdir -p logs

# Counter for submitted jobs
JOB_COUNT=0

echo "=========================================="
echo "Starting batch job submission"
echo "=========================================="
echo "CSV file: $CSV_FILE"
echo "Output directory: $OUTPUT_BASE_DIR"
echo "Annotation script: $SCRIPT_PATH"
echo "=========================================="

# Read CSV file and submit jobs (skip header)
tail -n +2 "$CSV_FILE" | while IFS=',' read -r sample_id sample_vcf cadd_file; do
    # Skip empty lines
    if [[ -z "$sample_id" ]]; then
        continue
    fi
    
    # Create sample-specific output directory
    SAMPLE_OUTPUT_DIR="$OUTPUT_BASE_DIR/$sample_id"
    mkdir -p "$SAMPLE_OUTPUT_DIR"
    
    # Check if input files exist
    if [[ ! -f "$sample_vcf" ]]; then
        echo "Warning: VCF file '$sample_vcf' not found for sample '$sample_id'. Skipping..."
        continue
    fi
    
    if [[ ! -f "$cadd_file" ]]; then
        echo "Warning: CADD file '$cadd_file' not found for sample '$sample_id'. Skipping..."
        continue
    fi
    
    echo "Submitting job for sample: $sample_id"
    echo "  VCF: $sample_vcf"
    echo "  CADD: $cadd_file"
    echo "  Output: $SAMPLE_OUTPUT_DIR"
    
    # Submit the job
    sbatch --job-name="${sample_id}_annotation" \
           "$SCRIPT_PATH" \
           -i "$sample_vcf" \
           -o "$SAMPLE_OUTPUT_DIR" \
           -c "$cadd_file" \
           -g "$GNOMAD_FILE" \
           -v "$CLINVAR_FILE" \
           -d "$VEP_DIR" \
           -s "$VEP_SIF"
    
    # Increment counter
    ((JOB_COUNT++))
    
    echo "  Job submitted successfully!"
    echo "---"
done

echo "=========================================="
echo "Batch submission completed!"
echo "Total jobs submitted: $JOB_COUNT"
echo "=========================================="
echo ""
echo "Monitor jobs with: squeue -u $USER"
echo "Check logs in: logs/"
echo "Results will be in: $OUTPUT_BASE_DIR/<sample_id>/"