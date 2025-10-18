# singularity exec \
#   ensamble_vep/vep.sif vep \
#   --dir ensamble_vep/vep_data/ \
#   -i data/debug_sample.vcf \
#   --cache --offline --vcf \
#   -o results/debug_sample.vcf \
#   --plugin CADD,snv=/data/humangen_kircherlab/Users/hassan/run_rare/rare-disease-pipeline/data/grch38/vcfanno/cadd/whole_genome_SNVs.tsv.gz,indels=/data/humangen_kircherlab/Users/hassan/chaos_lab/snv_annotations/databases/gnomad.genomes.r4.0.indel.tsv.gz

singularity exec \
  ensamble_vep/vep.sif vep \
  --dir ensamble_vep/vep_data/ \
  -i data/case_GS608_snv.chrY.vcf.gz \
  --cache --offline --vcf \
  -o results/case_GS608_snv.chrY.annotated.vcf \
  --custom /data/humangen_kircherlab/Users/hassan/data/databases/genomAD/gnomad.all.AC_AN_AF.vcf.bgz,gnomAD,vcf,exact,0,AC,AN,AF \
  --plugin CADD,snv=/data/humangen_kircherlab/Users/hassan/run_rare/rare-disease-pipeline/data/grch38/vcfanno/cadd/whole_genome_SNVs.tsv.gz,indels=/data/humangen_kircherlab/Users/hassan/chaos_lab/snv_annotations/databases/gnomad.genomes.r4.0.indel.tsv.gz

# Generate TSV output
singularity exec \
  ensamble_vep/vep.sif vep \
  --dir ensamble_vep/vep_data/ \
  -i data/case_GS608_snv.chrY.vcf.gz \
  --cache --offline --tab \
  -o results/case_GS608_snv.chrY.annotated.tsv \
  --custom databases/gnomad.genomes.v4.1.sites.chrY.vcf.bgz,gnomADg,vcf,exact,0,AF_grpmax \
  --plugin CADD,snv=/data/humangen_kircherlab/Users/hassan/run_rare/rare-disease-pipeline/data/grch38/vcfanno/cadd/whole_genome_SNVs.tsv.gz,indels=/data/humangen_kircherlab/Users/hassan/chaos_lab/snv_annotations/databases/gnomad.genomes.r4.0.indel.tsv.gz


# Extract clean table with relevant information
# Extract clean table with relevant information
echo "Extracting clean table..."
bcftools +split-vep -f '%CHROM\t%POS\t%REF\t%ALT\t%gnomADg\t%Consequence\t%SYMBOL\t%gnomADg_AF_grpmax\t%CADD_PHRED\t%CADD_RAW[\t%AD]\n' \
  -A tab \
  results/case_GS608_snv.chrY.annotated.vcf \
  | awk 'BEGIN {FS="\t"; OFS="\t"} 
  
  # Function to remove duplicate values from comma-separated string
  function unique_values(str) {
    if (str == "" || str == ".") return str
    n = split(str, arr, ",")
    delete seen
    result = ""
    for (i=1; i<=n; i++) {
      if (!(arr[i] in seen)) {
        seen[arr[i]] = 1
        result = result (result == "" ? "" : ",") arr[i]
      }
    }
    return result
  }
  
  NR==1 {
    print "Chromosome", "Position", "Ref", "Alt", "rsID", "Effect", "Gene", "gnomAD_AF_grpmax", "CADD_Phred", "CADD_Raw", "VAF"
    next
  }
  {
    # Calculate VAF from AD field
    if ($11 != ".") {
      split($11, ad, ",")
      total = ad[1] + ad[2]
      if (total > 0) {
        vaf = sprintf("%.3f", ad[2]/total)
      } else {
        vaf = "."
      }
    } else {
      vaf = "."
    }
    
    # Remove duplicates from comma-separated fields
    rsID = unique_values($5)
    effect = unique_values($6)
    gene = unique_values($7)
    gnomad_af = unique_values($8)
    cadd_phred = unique_values($9)
    cadd_raw = unique_values($10)
    
    print $1, $2, $3, $4, rsID, effect, gene, gnomad_af, cadd_phred, cadd_raw, vaf
  }' > results/case_GS608_snv.chrY.clean_table.tsv

echo "Clean table created: results/case_GS608_snv.chrY.clean_table.tsv"


## Filtering the clean table based on gnomAD_AF_grpmax < 0.005 and CADD_Phred > 20
echo "Filtering clean table..."
awk 'BEGIN {FS="\t"; OFS="\t"}
NR==1 {print; next}
{
  # Check gnomAD_AF_grpmax (column 8): missing or < 0.005
  gnomad_pass = ($8 == "." || $8 == "" || $8 < 0.005)
  
  # Check CADD_Phred (column 9): > 20 (allow missing values)
  cadd_pass = ($9 == "." || $9 == "" || $9 > 20)
  
  # Keep variant if both conditions are met
  if (gnomad_pass && cadd_pass) {
    print
  }
}' results/case_GS608_snv.chrY.clean_table.tsv > results/case_GS608_snv.chrY.clean_table.filtered.tsv

echo "Filtered clean table created: results/case_GS608_snv.chrY.clean_table.filtered.tsv"

echo "Clean table created: results/case_GS608_snv.chrY.clean_table.tsv"