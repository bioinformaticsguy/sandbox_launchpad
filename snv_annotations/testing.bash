#!/bin/bash

#SBATCH --partition=shortterm
#SBATCH --nodes=1
#SBATCH -c 4
#SBATCH --mem=16GB
#SBATCH --tmp=100G
#SBATCH --output=logs/%j_%u_%N_testing_slurmJob.out
#SBATCH --error=logs/%j_%u_%N_testing_slurmJob.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alihassan1697@gmail.com


module load singularity/v4.1.3

# # remove old results
# rm results/case_GS608_snv.all.genomad.tsv_summary.* 
# rm results/case_GS608_snv.all.genomad.*

# # Generate VCF output
# singularity exec \
#   ensamble_vep/vep.sif vep \
#   --dir ensamble_vep/vep_data/ \
#   --fork 4 \
#   -i /data/humangen_kircherlab/Users/hassan/chaos_lab/snv_annotations/data/case_GS608_snv.vcf.gz \
#   --cache --offline --vcf \
#   -o /data/humangen_kircherlab/Users/hassan/chaos_lab/snv_annotations/results/case608/case_GS608_snv.all.genomad.vcf \
#   --custom /data/humangen_kircherlab/Users/hassan/data/databases/genomAD/gnomad.all.AC_AN_AF.vcf.bgz,gnomAD,vcf,exact,0,AC,AN,AF 


# combine cadd and gnomad annotations
echo "Combining CADD and gnomAD annotations..."
echo "changing directory to results/case608"
cd /data/humangen_kircherlab/Users/hassan/chaos_lab/snv_annotations/results/case608

echo "printing working directory"
pwd


echo "Converting CADD TSV to VCF format..."
# Convert CADD TSV to VCF format for annotation
zcat case_GS608_snv_nochr.tsv.gz | awk 'BEGIN {OFS="\t"}
/^##/ {next}  # Skip CADD header lines
/^#/ {next}   # Skip column header
{
  # Add "chr" prefix if needed to match your gnomAD VCF
  chr = "chr" $1
  # Format: CHROM POS ID REF ALT QUAL FILTER INFO
  print chr, $2, ".", $3, $4, ".", ".", "CADD_RAW=" $5 ";CADD_PHRED=" $6
}' | \
# Sort the output properly before compressing
sort -k1,1V -k2,2n | \
bgzip -c > case_GS608_snv.cadd.vcf.gz

echo "Adding proper VCF header to CADD file..."
# Create a proper VCF with header
{
  echo "##fileformat=VCFv4.2"
  cat cadd_header.txt
  echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO"
  zcat case_GS608_snv.cadd.vcf.gz
} | bgzip -c > case_GS608_snv.cadd_with_header.vcf.gz

# Use the file with proper header
mv case_GS608_snv.cadd_with_header.vcf.gz case_GS608_snv.cadd.vcf.gz



echo "Indexing CADD VCF..."
# Index it
tabix -p vcf case_GS608_snv.cadd.vcf.gz

echo "Annotating gnomAD VCF with CADD scores..."
# Annotate gnomAD VCF with CADD scores
bcftools annotate \
  -a case_GS608_snv.cadd.vcf.gz \
  -c CHROM,POS,REF,ALT,INFO/CADD_RAW,INFO/CADD_PHRED \
  -h cadd_header.txt \
  case_GS608_snv.all.genomad.vcf \
  -Oz -o case_GS608_snv.all.genomad_cadd.vcf.gz


echo "Indexing combined VCF..."
# Index the final file
tabix -p vcf case_GS608_snv.all.genomad_cadd.vcf.gz

cd /data/humangen_kircherlab/Users/hassan/chaos_lab/snv_annotations/

# Generate VCF output with both gnomAD and CADD annotations
singularity exec \
  ensamble_vep/vep.sif vep \
  --dir ensamble_vep/vep_data/ \
  --fork 4 \
  -i /data/humangen_kircherlab/Users/hassan/chaos_lab/snv_annotations/results/case608/case_GS608_snv.all.genomad.vcf \
  --cache --offline --vcf \
  -o results/case608/case_GS608_snv.all.cadd_and_genomad_together.vcf \
  --custom results/case608/case_GS608_snv.cadd.vcf.gz,CADD,vcf,exact,0,CADD_RAW,CADD_PHRED \
  --custom /data/humangen_kircherlab/Users/hassan/data/databases/genomAD/gnomad.all.AC_AN_AF.vcf.bgz,gnomAD,vcf,exact,0,AC,AN,AF \
  --custom file=/data/humangen_kircherlab/Users/hassan/data/databases/clinvar/clinvar.vcf.gz,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN

echo "Combined file created: case_GS608_snv.all.genomad_cadd.vcf.gz"



