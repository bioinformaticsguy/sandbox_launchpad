#!/bin/bash

# Function to get reverse complement
reverse_complement() {
    echo "$1" | rev | tr "ACGTacgt" "TGCAtgca"
}

# Read the input file and create output
input_file="metadata.tsv"
output_file="metadata_combined.tsv"

# Write header
echo -e "sample_id\tbarcode" > "$output_file"

# Process each line (skip header)
tail -n +2 "$input_file" | while IFS=$'\t' read -r sample_id barcode_r1 barcode_r2; do
    # Get reverse complement of R2 barcode
    barcode_r2_rc=$(reverse_complement "$barcode_r2")
    
    # Combine R1 + R2_reverse_complement
    combined_barcode="${barcode_r1}${barcode_r2_rc}"
    
    # Write to output file
    echo -e "${sample_id}\t${combined_barcode}" >> "$output_file"
done

echo "Processing complete. Output saved to $output_file"
