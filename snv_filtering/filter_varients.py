from cyvcf2 import VCF, Writer
import sys
import os


def get_format_fields(variant):
    """Extract FORMAT fields for the first sample"""
    format_data = {}
    
    # Genotype (GT)
    if variant.genotypes:
        gt = variant.genotypes[0]  # [allele1, allele2, phased]
        format_data['GT'] = f"{gt[0]}/{gt[1]}" if not gt[2] else f"{gt[0]}|{gt[1]}"
    else:
        format_data['GT'] = 'NA'
    
    # Depth (DP)
    dp = variant.format('DP')
    format_data['DP'] = str(dp[0][0]) if dp is not None and len(dp) > 0 else 'NA'
    
    # Allelic depths (AD)
    ad = variant.format('AD')
    format_data['AD'] = ','.join(map(str, ad[0])) if ad is not None and len(ad) > 0 else 'NA'
    
    # Genotype Quality (GQ)
    gq = variant.format('GQ')
    format_data['GQ'] = str(gq[0][0]) if gq is not None and len(gq) > 0 else 'NA'
    
    # Phred-scaled likelihoods (PL)
    pl = variant.format('PL')
    format_data['PL'] = ','.join(map(str, pl[0])) if pl is not None and len(pl) > 0 else 'NA'
    
    # RNC (if exists)
    rnc = variant.format('RNC')
    format_data['RNC'] = ','.join(map(str, rnc[0])) if rnc is not None and len(rnc) > 0 else 'NA'
    
    return format_data

def get_csq_format(vcf):
    """Extract CSQ field format from VCF header"""
    for line in vcf.raw_header.split('\n'):
        if line.startswith('##INFO=<ID=CSQ'):
            # Extract format string
            format_start = line.find('Format: ') + 8
            format_end = line.find('">')
            return line[format_start:format_end]
    return None


def parse_csq_field(csq_string, csq_format):
    """Parse CSQ field and return list of dictionaries"""
    if not csq_string:
        return []
    
    # Get field names from CSQ format description
    field_names = csq_format.split('|')
    
    # Parse each CSQ entry (comma-separated)
    csq_entries = []
    for entry in csq_string.split(','):
        values = entry.split('|')
        csq_dict = dict(zip(field_names, values))
        csq_entries.append(csq_dict)
    
    return csq_entries


def filter_by_gnomad_af(input_vcf, output_prefix, af_threshold=0.005, cadd_threshold=None, additional_fields=None):
    """Filter variants by gnomAD_AF in CSQ field and optionally by CADD_PHRED"""
    
    if additional_fields is None:
        additional_fields = []
    
    vcf = VCF(input_vcf, strict_gt=True)
    csq_format = get_csq_format(vcf)
    
    if not csq_format:
        print("Error: CSQ field format not found in VCF header")
        return
    
    print(f"CSQ format: {csq_format}")
    field_names = csq_format.split('|')
    
    # Check if gnomAD_AF is in CSQ fields
    if 'gnomAD_AF' not in field_names:
        print("Error: gnomAD_AF not found in CSQ fields")
        return
    
    gnomad_af_idx = field_names.index('gnomAD_AF')
    print(f"gnomAD_AF index: {gnomad_af_idx}")
    
    # Get CADD_PHRED index if filtering on CADD
    cadd_phred_idx = None
    if 'CADD_CADD_PHRED' in field_names:
        cadd_phred_idx = field_names.index('CADD_CADD_PHRED')
        print(f"CADD_PHRED index: {cadd_phred_idx}")
        if cadd_threshold:
            print(f"Filtering on CADD_PHRED >= {cadd_threshold}")
    elif cadd_threshold:
        print("Warning: CADD_CADD_PHRED not found in CSQ fields, skipping CADD filtering")
        cadd_threshold = None
    
    # Verify additional fields exist
    if additional_fields:
        valid_fields = []
        for field in additional_fields:
            if field in field_names:
                valid_fields.append(field)
            else:
                print(f"Warning: Field '{field}' not found in CSQ")
        additional_fields = valid_fields
    
    # include AF and CADD thresholds in filename (replace '.' with 'p' to avoid extra dots)
    af_str = f"af{af_threshold}".replace('.', 'p')
    cadd_str = f"cadd{cadd_threshold}".replace('.', 'p') if cadd_threshold is not None else "noCADD"
    suffix = f".{af_str}.{cadd_str}"
    output_vcf = f"{output_prefix}{input_basename}{suffix}.rare_variants.vcf"
    print(f"Output VCF: {output_vcf}")
    writer = Writer(output_vcf, vcf)
    
    # Statistics
    total_variants = 0
    rare_variants = 0
    missing_af = 0
    failed_cadd = 0
    
    for variant in vcf:
        total_variants += 1
        
        # Get CSQ field
        csq_field = variant.INFO.get('CSQ')
        
        if not csq_field:
            continue
        
        # Parse CSQ entries
        is_rare = False
        passes_cadd = True
        
        for csq_entry in csq_field.split(','):
            values = csq_entry.split('|')
            
            if len(values) > gnomad_af_idx:
                gnomad_af_value = values[gnomad_af_idx]
                
                # Check if AF is missing or below threshold
                if gnomad_af_value == '' or gnomad_af_value == '.':
                    missing_af += 1
                    is_rare = True
                else:
                    try:
                        af = float(gnomad_af_value)
                        if af < af_threshold:
                            is_rare = True
                    except ValueError:
                        print(f"Warning: Could not parse AF value '{gnomad_af_value}' at {variant.CHROM}:{variant.POS}")
                        continue
                
                # Check CADD threshold if specified
                if is_rare and cadd_threshold and cadd_phred_idx is not None:
                    if len(values) > cadd_phred_idx:
                        cadd_value = values[cadd_phred_idx]
                        if cadd_value and cadd_value != '' and cadd_value != '.':
                            try:
                                cadd_score = float(cadd_value)
                                if cadd_score < cadd_threshold:
                                    passes_cadd = False
                                    failed_cadd += 1
                            except ValueError:
                                passes_cadd = False
                        else:
                            passes_cadd = False
                
                if is_rare and passes_cadd:
                    break
        
        if is_rare and passes_cadd:
            rare_variants += 1
            writer.write_record(variant)
    
    writer.close()
    vcf.close()
    
    # Print statistics
    print(f"\n=== Filtering Results ===")
    print(f"Total variants: {total_variants}")
    print(f"Rare variants (AF < {af_threshold}): {rare_variants}")
    print(f"Variants with missing AF: {missing_af}")
    if cadd_threshold:
        print(f"Variants filtered out by CADD < {cadd_threshold}: {failed_cadd}")
    print(f"Final variants passing all filters: {rare_variants}")
    print(f"Output vcf: {output_vcf}")
    
    # Also create a TSV summary
    create_summary_table(input_vcf, output_prefix, af_threshold, csq_format, cadd_threshold, additional_fields)


def create_summary_table(input_vcf, output_prefix, af_threshold, csq_format, cadd_threshold=None, additional_fields=None):
    """Create a TSV file with filtered variants"""
    
    if additional_fields is None:
        additional_fields = []
    
    vcf = VCF(input_vcf)
    field_names = csq_format.split('|')
    gnomad_af_idx = field_names.index('gnomAD_AF')
    
    # Try to get CADD indices
    cadd_raw_idx = field_names.index('CADD_CADD_RAW') if 'CADD_CADD_RAW' in field_names else None
    cadd_phred_idx = field_names.index('CADD_CADD_PHRED') if 'CADD_CADD_PHRED' in field_names else None
    
    # Get indices for additional fields
    additional_indices = {}
    for field in additional_fields:
        if field in field_names:
            additional_indices[field] = field_names.index(field)
    

    # Create output TSV    
    input_basename_local = os.path.splitext(os.path.basename(input_vcf))[0]
    if input_basename_local.endswith('.vcf'):
        input_basename_local = os.path.splitext(input_basename_local)[0]
    af_str = f"af{af_threshold}".replace('.', 'p')
    cadd_str = f"cadd{cadd_threshold}".replace('.', 'p') if cadd_threshold is not None else "noCADD"
    suffix = f".{af_str}.{cadd_str}"
    output_tsv = f"{output_prefix}{input_basename_local}{suffix}.rare_variants.tsv"
    print(f"Output TSV: {output_tsv}")

    with open(output_tsv, 'w') as f:
        # Write header
        # header = ['CHROM', 'POS', 'REF', 'ALT', 'gnomAD_AF']
        header = ['CHROM', 'POS', 'REF', 'ALT', 'GT', 'DP', 'AD', 'GQ', 'PL', 'RNC', 'gnomAD_AF']
        if cadd_raw_idx is not None:
            header.append('CADD_RAW')
        if cadd_phred_idx is not None:
            header.append('CADD_PHRED')
        # Add additional fields to header
        header.extend(additional_fields)
        f.write('\t'.join(header) + '\n')
        
        for variant in vcf:
            csq_field = variant.INFO.get('CSQ')
            
            if not csq_field:
                continue
            
            for csq_entry in csq_field.split(','):
                values = csq_entry.split('|')
                
                if len(values) > gnomad_af_idx:
                    gnomad_af_value = values[gnomad_af_idx]
                    
                    # Check if rare
                    is_rare = False
                    passes_cadd = True
                    
                    if gnomad_af_value == '' or gnomad_af_value == '.':
                        is_rare = True
                        gnomad_af_value = 'NA'
                    else:
                        try:
                            af = float(gnomad_af_value)
                            if af < af_threshold:
                                is_rare = True
                        except ValueError:
                            continue
                    
                    # Check CADD if threshold specified
                    if is_rare and cadd_threshold and cadd_phred_idx:
                        if len(values) > cadd_phred_idx:
                            cadd_value = values[cadd_phred_idx]
                            if cadd_value and cadd_value != '' and cadd_value != '.':
                                try:
                                    cadd_score = float(cadd_value)
                                    if cadd_score < cadd_threshold:
                                        passes_cadd = False
                                except ValueError:
                                    passes_cadd = False
                            else:
                                passes_cadd = False
                    
                    if is_rare and passes_cadd:
                        # Get FORMAT fields
                        format_data = get_format_fields(variant)

                        row = [
                            variant.CHROM,
                            str(variant.POS),
                            variant.REF,
                            ','.join(variant.ALT),
                            format_data['GT'],
                            format_data['DP'],
                            format_data['AD'],
                            format_data['GQ'],
                            format_data['PL'],
                            format_data['RNC'],
                            gnomad_af_value
                            ]
                        
                        if cadd_raw_idx is not None and len(values) > cadd_raw_idx:
                            row.append(values[cadd_raw_idx] or 'NA')
                        if cadd_phred_idx is not None and len(values) > cadd_phred_idx:
                            row.append(values[cadd_phred_idx] or 'NA')
                        
                        # Add additional fields
                        for field in additional_fields:
                            if field in additional_indices:
                                idx = additional_indices[field]
                                if len(values) > idx:
                                    row.append(values[idx] or 'NA')
                                else:
                                    row.append('NA')
                            else:
                                row.append('NA')
                        
                        f.write('\t'.join(row) + '\n')
                        break  # Only write once per variant
    
    vcf.close()
    print(f"Summary table: {output_tsv}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python filter_varients.py <input_vcf> <output_prefix> [af_threshold] [cadd_threshold] [additional_fields]")
        print("\nArguments:")
        print("  input_vcf          : Input VCF file (can be .vcf or .vcf.gz)")
        print("  output_prefix      : Prefix for output files")
        print("  af_threshold       : gnomAD allele frequency threshold (default: 0.005)")
        print("  cadd_threshold     : CADD_PHRED score threshold (optional, e.g., 15)")
        print("  additional_fields  : Comma-separated list of CSQ fields to extract (e.g., 'Gene,Consequence,SYMBOL')")
        print("\nExamples:")
        print("  python filter_varients.py input.vcf.gz output")
        print("  python filter_varients.py input.vcf.gz output 0.005")
        print("  python filter_varients.py input.vcf.gz output 0.005 15")
        print("  python filter_varients.py input.vcf.gz output 0.005 15 'Gene,Consequence,SYMBOL,gnomAD_AC,gnomAD_AN'")
        sys.exit(1)
    
    input_vcf = sys.argv[1]
    output_prefix = sys.argv[2]
    af_threshold = float(sys.argv[3]) if len(sys.argv) > 3 else 0.005
    cadd_threshold = float(sys.argv[4]) if len(sys.argv) > 4 else None

    input_basename = os.path.splitext(os.path.basename(input_vcf))[0]
    if input_basename.endswith('.vcf'):
        input_basename = os.path.splitext(input_basename)[0]

    
    # Parse additional CSQ fields to extract (comma-separated)
    additional_fields = []
    if len(sys.argv) > 5:
        additional_fields = [field.strip() for field in sys.argv[5].split(',')]
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_prefix)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    print(f"Input VCF: {input_vcf}")
    print(f"Output prefix: {output_prefix}")
    print(f"gnomAD AF threshold: < {af_threshold}")
    if cadd_threshold:
        print(f"CADD_PHRED threshold: >= {cadd_threshold}")
    else:
        print("CADD_PHRED threshold: Not applied")
    if additional_fields:
        print(f"Additional CSQ fields to extract: {', '.join(additional_fields)}")
    print()
    
    filter_by_gnomad_af(input_vcf, output_prefix, af_threshold, cadd_threshold, additional_fields)


    ## get max output

    ##  python3 filter_varients.py /data/humangen_kircherlab/Users/hassan/chaos_lab/temp_files/case_GS608_snv.chr21.annotated.vcf /data/humangen_kircherlab/Users/hassan/chaos_lab/temp_files/filtered 0.005 20 'Allele,Consequence,IMPACT,SYMBOL,Gene,Feature_type,Feature,BIOTYPE,EXON,INTRON,HGVSc,HGVSp,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,DISTANCE,STRAND,FLAGS,SYMBOL_SOURCE,HGNC_ID,SOURCE,CADD,CADD_CADD_RAW,CADD_CADD_PHRED,gnomAD,gnomAD_AC,gnomAD_AN,gnomAD_AF,ClinVar,ClinVar_CLNSIG,ClinVar_CLNREVSTAT,ClinVar_CLNDN'


    # Usefull
    # 'Consequence,IMPACT,SYMBOL,Gene,Feature_type,Feature,BIOTYPE,SYMBOL_SOURCE,HGNC_ID'


    ## GT Info
    ## GT:DP:AD:GQ:PL:RNC	./1:40:.,12:.:0,0,0:O.  GT will be -1/1 it is because there was a . in GT.