#!/usr/bin/env python3
import pandas as pd
import sys

def main(input_file, output_file):
    # Define the column names for the GFF file
    gff_column_names = [
        "seqname", "source", "feature", "Start", "End", "score", "Strand", "frame", "Attributes"
    ]

    # Read the GFF file with the specified column names
    df = pd.read_csv(input_file, sep='\t', names=gff_column_names, comment='#', header=None)
    print(df.head())
    # Initialize variables
    ensembl_gene_ids = []
    current_ensg = "NA"

    # Loop through the GFF file to extract gene IDs and add them to the DataFrame
    for idx, row in df.iterrows():
        if row['feature'] == 'gene':
            current_ensg = row['Attributes'].split('ID=gene:')[1].split(';')[0]
            
        ensembl_gene_ids.append(current_ensg)
        if row['Attributes'].startswith('###'):
            current_ensg = "NA"  # Reset the gene ID at the end of a block

    # Add the ensembl_gene_id column to the DataFrame
    df['ensembl_gene_id'] = ensembl_gene_ids
    print(df)

    # Group by GeneID and apply a lambda function to select the most 3' transcript
    def select_most_3_transcript(group):
        non_gene_features = group[group['feature'] == 'exon']
        if not non_gene_features.empty:
            if group['Strand'].iloc[0] == '+':
                return non_gene_features.loc[non_gene_features['End'].idxmax()]
            else:
                return non_gene_features.loc[non_gene_features['Start'].idxmin()]
        else:
            return None

    result = df.groupby('ensembl_gene_id').apply(select_most_3_transcript)
    # Filter out None values
    result = result.dropna().reset_index(drop=True)

    # Save the result to a new file
    result.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python 3primeGrab.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    main(input_file, output_file)
