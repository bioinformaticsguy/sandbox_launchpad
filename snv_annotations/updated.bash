#!/bin/bash

#SBATCH --partition=shortterm
#SBATCH --nodes=1
#SBATCH -c 16
#SBATCH --mem=64GB
#SBATCH --tmp=100G
#SBATCH --output=logs/%j_%u_%N_testing_slurmJob.out
#SBATCH --error=logs/%j_%u_%N_testing_slurmJob.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alihassan1697@gmail.com

# Default values for annotation databases
DEFAULT_CADD="/data/humangen_kircherlab/Users/hassan/chaos_lab/snv_annotations/data/case_GS608_snv_nochr.tsv.gz"
DEFAULT_GNOMAD="/data/humangen_kircherlab/Users/hassan/data/databases/genomAD/gnomad.all.AC_AN_AF.vcf.bgz"
DEFAULT_CLINVAR="/data/humangen_kircherlab/Users/hassan/data/databases/clinvar/clinvar.vcf.gz"
DEFAULT_VEP_DIR="ensamble_vep/vep_data/"
DEFAULT_VEP_SIF="ensamble_vep/vep.sif"

# Function to display usage
usage() {
    cat << EOF
Usage: $0 -i <input_vcf> -o <output_dir> [OPTIONS]

Required arguments:
  -i, --input       Input VCF file (can be .vcf or .vcf.gz)
  -o, --output      Output directory for results

Optional arguments:
  -c, --cadd        CADD annotation file (default: $DEFAULT_CADD)
  -g, --gnomad      gnomAD annotation file (default: $DEFAULT_GNOMAD)
  -v, --clinvar     ClinVar annotation file (default: $DEFAULT_CLINVAR)
  -d, --vep-dir     VEP data directory (default: $DEFAULT_VEP_DIR)
  -s, --vep-sif     VEP singularity image (default: $DEFAULT_VEP_SIF)
  -k, --keep-temp   Keep intermediate files (default: remove them)
  -h, --help        Display this help message

Example:
  $0 -i data/variants.vcf.gz -o results/my_analysis
  $0 -i data/variants.vcf.gz -o results/my_analysis -c custom_cadd.vcf.gz -k
  
Real Example:
  sbatch --job-name=testing_annotations updated.bash -i /data/humangen_kircherlab/Users/hassan/chaos_lab/snv_annotations/data/subset_1000.vcf.gz -o /data/humangen_kircherlab/Users/hassan/chaos_lab/temp_files
EOF
    exit 1
}

# Initialize variables
INPUT_VCF=""
OUTPUT_DIR=""
CADD_FILE="$DEFAULT_CADD"
GNOMAD_FILE="$DEFAULT_GNOMAD"
CLINVAR_FILE="$DEFAULT_CLINVAR"
VEP_DIR="$DEFAULT_VEP_DIR"
VEP_SIF="$DEFAULT_VEP_SIF"
KEEP_TEMP=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            INPUT_VCF="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -c|--cadd)
            CADD_FILE="$2"
            shift 2
            ;;
        -g|--gnomad)
            GNOMAD_FILE="$2"
            shift 2
            ;;
        -v|--clinvar)
            CLINVAR_FILE="$2"
            shift 2
            ;;
        -d|--vep-dir)
            VEP_DIR="$2"
            shift 2
            ;;
        -s|--vep-sif)
            VEP_SIF="$2"
            shift 2
            ;;
        -k|--keep-temp)
            KEEP_TEMP=true
            shift
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
if [[ -z "$INPUT_VCF" ]] || [[ -z "$OUTPUT_DIR" ]]; then
    echo "Error: Input VCF and output directory are required."
    usage
fi

# Check if input file exists
if [[ ! -f "$INPUT_VCF" ]]; then
    echo "Error: Input file '$INPUT_VCF' not found."
    exit 1
fi

# Load required modules
module load singularity/v4.1.3

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Extract base name without extension for output files
BASENAME=$(basename "$INPUT_VCF" .vcf.gz)
BASENAME=$(basename "$BASENAME" .vcf)

echo "=========================================="
echo "Starting variant annotation pipeline"
echo "=========================================="
echo "Input VCF: $INPUT_VCF"
echo "Output directory: $OUTPUT_DIR"
echo "Base name: $BASENAME"
echo "CADD file: $CADD_FILE"
echo "gnomAD file: $GNOMAD_FILE"
echo "ClinVar file: $CLINVAR_FILE"
echo "=========================================="

# Define output files
CADD_VCF="$OUTPUT_DIR/${BASENAME}.cadd_prepared.vcf.gz"
FINAL_VCF="$OUTPUT_DIR/${BASENAME}.annotated.vcf"

# Step 1: Prepare CADD annotations for VEP
echo "Step 1: Preparing CADD annotations for VEP..."
if [[ "$CADD_FILE" == *.tsv.gz ]]; then
    echo "Converting CADD TSV to VCF format..."
    
    # Convert TSV to VCF with proper header
    {
        echo "##fileformat=VCFv4.2"
        echo "##INFO=<ID=CADD_RAW,Number=1,Type=Float,Description=\"CADD raw score\">"
        echo "##INFO=<ID=CADD_PHRED,Number=1,Type=Float,Description=\"CADD PHRED score\">"
        echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO"
        zcat "$CADD_FILE" | awk 'BEGIN {OFS="\t"}
        /^##/ {next}
        /^#/ {next}
        {
            chr = ($1 ~ /^chr/) ? $1 : "chr" $1
            print chr, $2, ".", $3, $4, ".", ".", "CADD_RAW=" $5 ";CADD_PHRED=" $6
        }' | sort -k1,1V -k2,2n
    } | bgzip -c > "$CADD_VCF"
    
    echo "Indexing CADD VCF..."
    tabix -p vcf "$CADD_VCF"
    CADD_PREPARED="$CADD_VCF"
else
    # Assume CADD_FILE is already in VCF format
    echo "CADD file is already in VCF format, using directly..."
    CADD_PREPARED="$CADD_FILE"
fi

# Step 2: Run VEP with all custom annotations
echo "Step 2: Running VEP with CADD, gnomAD, and ClinVar annotations..."
singularity exec \
  "$VEP_SIF" vep \
  --dir "$VEP_DIR" \
  --fork 4 \
  -i "$INPUT_VCF" \
  --cache --offline --vcf \
  -o "$FINAL_VCF" \
  --custom "file=$CADD_PREPARED,short_name=CADD,format=vcf,type=exact,coords=0,fields=CADD_RAW%CADD_PHRED" \
  --custom "file=$GNOMAD_FILE,short_name=gnomAD,format=vcf,type=exact,coords=0,fields=AC%AN%AF" \
  --custom "file=$CLINVAR_FILE,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN%CLNDISDB"

# Step 3: Clean up intermediate files
if [[ "$KEEP_TEMP" == false ]]; then
    echo "Cleaning up intermediate files..."
    if [[ "$CADD_FILE" == *.tsv.gz ]]; then
        rm -f "$CADD_VCF" "$CADD_VCF.tbi"
    fi
fi

echo "=========================================="
echo "Annotation pipeline completed successfully!"
echo "Final output: $FINAL_VCF"
echo "=========================================="