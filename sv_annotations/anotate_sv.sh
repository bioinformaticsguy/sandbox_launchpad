
# singularity exec \
#   ensamble_vep/vep.sif vep \
#   --dir ensamble_vep/vep_data/ \
#   -i data/case_GS608_snv.chrY.vcf.gz \
#   --cache --offline --vcf \
#   -o results/case_GS608_snv.chrY.annotated.vcf \
#   --custom /data/humangen_kircherlab/Users/hassan/data/databases/genomAD/gnomad.all.AC_AN_AF.vcf.bgz,gnomAD,vcf,exact,0,AC,AN,AF \
#   --plugin CADD,snv=/data/humangen_kircherlab/Users/hassan/run_rare/rare-disease-pipeline/data/grch38/vcfanno/cadd/whole_genome_SNVs.tsv.gz,indels=/data/humangen_kircherlab/Users/hassan/chaos_lab/snv_annotations/databases/gnomad.genomes.r4.0.indel.tsv.gz


#   --custom "file=$CADD_PREPARED,short_name=CADD,format=vcf,type=exact,coords=0,fields=CADD_RAW%CADD_PHRED" \

# Step 2: Run VEP with all custom annotations
echo "Step 2: Running VEP with CADD and gnomAD annotations..."
singularity exec \
  /data/humangen_kircherlab/Users/hassan/repos/sandbox_launchpad/snv_annotations/ensamble_vep/vep.sif vep \
  --dir /data/humangen_kircherlab/Users/hassan/repos/sandbox_launchpad/snv_annotations/ensamble_vep/vep_data \
  --fork 4 \
  -i /data/humangen_kircherlab/Users/hassan/repos/sandbox_launchpad/sv_annotations/data/case_GS608_sv.vcf.gz \
  --cache --offline --vcf \
  -o results/case_GS608_sv.annotated.vcf \
  --custom "file=/data/humangen_kircherlab/Users/hassan/repos/sandbox_launchpad/sv_annotations/data/case_GS608_sv_score.bed.gz,short_name=CADD_SV,format=bed,type=exact,coords=0,fields=CADD-SV_PHRED-score%CADD-SV_Raw-score" \
  --custom "file=/data/humangen_kircherlab/Users/hassan/data/databases/genomAD/gnomad.v4.1.sv.non_neuro_controls.sites.vcf.gz,short_name=gnomAD,format=vcf,type=exact,coords=0,fields=AC%AN%AF" 
