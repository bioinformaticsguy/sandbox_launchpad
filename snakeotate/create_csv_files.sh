#!/bin/bash

# Function to display usage
usage() {
    echo "Usage: $0 <sample_path> <output_path> <batch_size>"
    echo "  sample_path: Directory containing all sample directories"
    echo "  output_path: Directory where CSV files will be created"
    echo "  batch_size:  Number of samples per CSV file"
    echo ""
    echo "Example: $0 /data/samples /data/output 5"
    exit 1
}

# Check if correct number of arguments provided
if [ "$#" -ne 3 ]; then
    echo "Error: Incorrect number of arguments"
    usage
fi

SAMPLE_PATH="$1"
OUTPUT_PATH="$2"
BATCH_SIZE="$3"

# Validate inputs
if [ ! -d "$SAMPLE_PATH" ]; then
    echo "Error: Sample path '$SAMPLE_PATH' does not exist or is not a directory"
    exit 1
fi

if [ ! -d "$OUTPUT_PATH" ]; then
    echo "Creating output directory: $OUTPUT_PATH"
    mkdir -p "$OUTPUT_PATH"
fi

if ! [[ "$BATCH_SIZE" =~ ^[0-9]+$ ]] || [ "$BATCH_SIZE" -le 0 ]; then
    echo "Error: Batch size must be a positive integer"
    exit 1
fi

echo "Sample path: $SAMPLE_PATH"
echo "Output path: $OUTPUT_PATH"
echo "Batch size: $BATCH_SIZE"
echo ""

# Get list of sample directories (only directories, not files)
SAMPLES=($(find "$SAMPLE_PATH" -maxdepth 1 -type d -not -path "$SAMPLE_PATH" -printf "%f\n" | sort))

if [ ${#SAMPLES[@]} -eq 0 ]; then
    echo "Error: No sample directories found in $SAMPLE_PATH"
    exit 1
fi

echo "Found ${#SAMPLES[@]} samples:"
printf '%s\n' "${SAMPLES[@]}" | head -10
if [ ${#SAMPLES[@]} -gt 10 ]; then
    echo "... and $((${#SAMPLES[@]} - 10)) more"
fi
echo ""

# Calculate number of batches needed
TOTAL_SAMPLES=${#SAMPLES[@]}
NUM_BATCHES=$(( (TOTAL_SAMPLES + BATCH_SIZE - 1) / BATCH_SIZE ))

echo "Creating $NUM_BATCHES batch files with $BATCH_SIZE samples each (last batch may have fewer)"
echo ""

# Create batch CSV files
for ((batch=1; batch<=NUM_BATCHES; batch++)); do
    # Calculate start and end indices for this batch
    START_IDX=$(( (batch - 1) * BATCH_SIZE ))
    END_IDX=$(( START_IDX + BATCH_SIZE - 1 ))
    
    # Adjust end index if it exceeds array bounds
    if [ $END_IDX -ge $TOTAL_SAMPLES ]; then
        END_IDX=$(( TOTAL_SAMPLES - 1 ))
    fi
    
    # Create CSV filename
    CSV_FILE="$OUTPUT_PATH/batch_${batch}_sample_list.csv"
    
    # Write CSV header
    echo "sample_name,file_path" > "$CSV_FILE"
    
    # Write samples for this batch
    BATCH_COUNT=0
    for ((i=START_IDX; i<=END_IDX; i++)); do
        SAMPLE_NAME="${SAMPLES[i]}"
        SAMPLE_FULL_PATH="$SAMPLE_PATH/$SAMPLE_NAME"
        echo "$SAMPLE_NAME,$SAMPLE_FULL_PATH" >> "$CSV_FILE"
        ((BATCH_COUNT++))
    done
    
    echo "Created $CSV_FILE with $BATCH_COUNT samples"
done

echo ""
echo "Batch creation completed!"
echo "Generated files:"
ls -la "$OUTPUT_PATH"/*.csv