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



# Command used to create the updated fastq files with combined barcodes:

# zcat */*_1.fq.gz | awk 'BEGIN{ OFS="\t"; count=0}{ if (count == 0) { split($2,a,":"); split(a[4],b,"+"); print $1 }; if (count == 1) { print $1""b[1] } if (count == 2) print; if (count == 3) { print $1"IIIIIIIIII"; } count+=1; if (count > 3) count = 0;}' | gzip -c > combined_R1.fq.gz

# zcat */*_1.fq.gz  | awk 'BEGIN{ OFS="\t"; count=0}{ if (count == 0) { split($2,a,":"); split(a[4],b,"+"); print $1 }; if (count == 1) { print $1""b[2] } if (count == 2) print; if (count == 3) { print $1"IIIIIIIIII"; } count+=1; if (count > 3) count = 0;}' | gzip -c > combined_R2.fq.gz

# Hi Martin, I was getting a bit confused because of line 3 it should like this right:

# zcat */*_2.fq.gz  | awk 'BEGIN{ OFS="\t"; count=0}{ if (count == 0) { split($2,a,":"); split(a[4],b,"+"); print $1 }; if (count == 1) { print $1""b[2] } if (count == 2) print; if (count == 3) { print $1"IIIIIIIIII"; } count+=1; if (count > 3) count = 0;}' | gzip -c > combined_R2.fq.gz

# for the other part of barcode we need to apply that only to R2 files... but in the commands that you shared it is apparently applied to R1 in both cases.

