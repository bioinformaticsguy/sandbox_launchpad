import pandas as pd
import argparse
import os
from pathlib import Path


def annotate_variants_with_hpo(input_file, hpo_df, hpo_ids, output_directory):
    """
    Annotate a single variant file with HPO information.
    
    Args:
        input_file: Path to the input TSV file
        hpo_df: DataFrame containing HPO annotations
        hpo_ids: List of HPO IDs to filter by
        output_directory: Directory to save the annotated file
    """
    # Read the filtered TSV file to get the dataframe for annotation
    variants_df = pd.read_csv(input_file, sep='\t')
    
    # Filter the HPO dataframe for rows where 'hpo_id' is in hpo_ids
    filtered_df = hpo_df[hpo_df['hpo_id'].isin(hpo_ids)]
    
    # Add columns hpo_id, hpo_name and disease_id to variants_df
    variants_df = variants_df.merge(
        filtered_df[['gene_symbol', 'hpo_id', 'hpo_name', 'disease_id']], 
        left_on='SYMBOL', 
        right_on='gene_symbol', 
        how='left'
    )
    
    # Get base filename without extension
    base_filename = Path(input_file).stem
    
    # Export the annotated variants dataframe to a new TSV file
    output_file = os.path.join(output_directory, f'{base_filename}.hpo_annotated.tsv')
    variants_df.to_csv(output_file, sep='\t', index=False)
    
    print(f'Annotated variants saved to {output_file}')
    return output_file


def main():
    parser = argparse.ArgumentParser(description='Annotate variant TSV files with HPO information')
    parser.add_argument('-i', '--input_dir', required=True, 
                        help='Directory containing input TSV files')
    parser.add_argument('-o', '--output_dir', required=True,
                        help='Directory to save annotated output files')
    parser.add_argument('--hpo_file', 
                        default='/data/humangen_kircherlab/Users/hassan/chaos_lab/hpo_annotations/genes_to_phenotype.txt',
                        help='Path to genes_to_phenotype.txt file')
    parser.add_argument('--hpo_ids', nargs='+', default=['HP:0000078'],
                        help='List of HPO IDs to filter by (space-separated)')
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Read genes_to_phenotype.txt as a dataframe
    print(f'Loading HPO annotations from {args.hpo_file}')
    hpo_df = pd.read_csv(args.hpo_file, sep='\t')
    
    # Find all TSV files in the input directory
    input_dir = Path(args.input_dir)
    tsv_files = list(input_dir.glob('*.tsv'))
    
    if not tsv_files:
        print(f'No TSV files found in {args.input_dir}')
        return
    
    print(f'Found {len(tsv_files)} TSV file(s) to process')
    
    # Process each TSV file
    for tsv_file in tsv_files:
        print(f'\nProcessing {tsv_file.name}...')
        try:
            annotate_variants_with_hpo(
                str(tsv_file), 
                hpo_df, 
                args.hpo_ids, 
                args.output_dir
            )
        except Exception as e:
            print(f'Error processing {tsv_file.name}: {str(e)}')
            continue
    
    print(f'\nAll files processed. Output saved to {args.output_dir}')


if __name__ == '__main__':
    main()