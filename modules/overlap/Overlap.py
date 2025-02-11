"""
This script processes a DataFrame of genomic features to identify transcripts that end 3' of their corresponding genes.

Steps:
1. Load the DataFrame from a CSV file.
2. Extract transcript IDs from the 'Attributes' column.
3. Use the BioMart service to fetch gene information for each transcript, including gene ID, start position, end position, and strand.
4. Merge the gene information with the original DataFrame.
5. Filter the DataFrame to include only rows where the transcript ends 3' of the gene:
   - For genes on the '+' strand, the transcript end position must be greater than the gene end position.
   - For genes on the '-' strand, the transcript start position must be less than the gene start position.
"""
import argparse
import io
import pandas as pd
from io import StringIO
import time
from bisect import bisect_left, bisect_right

def findIntersectingEnsgs(pipeline_res_df, gff_file):

    gff_file["start"] = pd.to_numeric(gff_file["start"], errors='coerce')
    gff_file["end"] = pd.to_numeric(gff_file["end"], errors='coerce')
    matching_genes = []
    chromosomes = pipeline_res_df["Chromosome"].unique()

    # Ensure chromosomes is always a list
    # Ensure chromosomes is always a list of strings
    chromosomes = [str(chromosome) for chromosome in chromosomes]
    for chromosome in chromosomes:
        chromosome_df = pipeline_res_df[pipeline_res_df["Chromosome"] == chromosome]
        ensembl_df = gff_file
        start_positions = ensembl_df["start"].astype(int).tolist()
        end_positions = ensembl_df["end"].astype(int).tolist()
        #print(ensembl_df)
        for index, row in chromosome_df.iterrows():
            strand = row["Strand"].strip()  # Remove any leading/trailing whitespace
            current_ensg = row["ensembl_gene_id"]
            polyA_start = int(row["PolyA_Start"])

            if strand == "+":
                start_idx = bisect_left(start_positions, polyA_start)
                
                
                matching_indices = []
                for i in range(start_idx):
                    if start_positions[i] <= polyA_start <= end_positions[i]:
                        #print(start_positions[i], polyA_start, end_positions[i], "all positions")
                        matching_indices.append(i)

                # Process all matching indices
                for matching_index in matching_indices:
                    result_row = ensembl_df.iloc[matching_index]
                    #print(result_row, "RESULTROW")
                    start_position = int(result_row["start"])
                    end_position = int(result_row["end"])

                    if (str(result_row["chromosome_name"]) == str(chromosome) and
                        str(result_row["strand"]) == str(strand) and
                        polyA_start >= start_position and polyA_start <= end_position):
                        #print("MATCH")
                        if result_row["ensembl_gene_id"] != current_ensg:
                            matching_genes.append((result_row["ensembl_gene_id"], current_ensg, chromosome, strand, start_position, end_position, polyA_start))
            elif strand == "-":
                start_idx = bisect_left(start_positions, polyA_start)
                # Iterate through potential matches to find the correct one
                matching_indices = []
                for i in range(start_idx):
                    if start_positions[i] <= polyA_start <= end_positions[i]:
                        matching_indices.append(i)

                # Process all matching indices
                for matching_index in matching_indices:
                    result_row = ensembl_df.iloc[matching_index]
                    start_position = int(result_row["start"])
                    end_position = int(result_row["end"])

                    if (str(result_row["chromosome_name"]) == str(chromosome) and
                        #str(result_row["strand"]) == str(strand) and
                        polyA_start >= start_position and polyA_start <= end_position):
                        #print("MATCH")
                        if result_row["ensembl_gene_id"] != current_ensg:
                            matching_genes.append((result_row["ensembl_gene_id"], current_ensg, chromosome, strand, start_position, end_position, polyA_start))
    return pd.DataFrame(matching_genes, columns=["Matching_ENSG", "Current_ENSG", "Chromosome", "Strand", "Start_Position", "End_Position", "polyApos"])
# Function to extract gene name from the attribute column
def extract_gene_name(attribute):
    for item in attribute.split(';'):
        if item.startswith('ID=gene:'):
            return item.split(':')[1]
    return None

def main(pipeline_input, gff_input, output_file):
    column_names = [
        "chromosome_name", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"
    ]

    # Read the GFF file with the specified column names
    gff_file = pd.read_csv(gff_input, sep="\t", comment='#', names=column_names)
    # Filter rows where feature is 'gene'
    gene_rows = gff_file[gff_file['feature'] == 'gene']
    gene_rows['ensembl_gene_id'] = gene_rows['attribute'].apply(extract_gene_name)
    gff_file = gene_rows
    filter_file = pd.read_csv(pipeline_input, sep='\t', low_memory=False, header=None)
    print(filter_file.head())
    # Convert the result to a DataFrame
    filter_file.columns=["Chromosome", "Source", "Feature",  "Start", "End", "Score", "Strand", "Phase", "Attributes", "ensembl_gene_id", "PolyA_Start", "PolyA_End", "Fantom_Start", "Fantom_End", "Fantom_Name"]
    print(filter_file.head())
    print(gff_file.head())
    # Save the DataFrame to a CSV file
    result = findIntersectingEnsgs(filter_file, gff_file)
    result.to_csv(output_file, index=False)
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Find intersecting ENSGs.')
    parser.add_argument('pipeline_input', help='Path to the input CSV file')
    parser.add_argument('gff_input', help='Path to the GFF input file')
    parser.add_argument('output_file', help='Path to the output CSV file')
    args = parser.parse_args()
    main(args.pipeline_input, args.gff_input, args.output_file)


