"""
This script filters a DataFrame of genomic features based on whether the transcript ends 3' of the gene it belongs to. 

It first loads the DataFrame from a CSV file, then extracts the transcript IDs from the 'Attributes' column. 

It uses the BioMart service to fetch gene information for each transcript, including the gene ID, start position, end position, and strand. 

The gene information is then merged with the original DataFrame. 

Finally, it filters the DataFrame to only include rows where the transcript ends 3' of the gene, i.e., for genes on the '+' strand, the transcript end position is greater than the gene end position, and for genes on the '-' strand, the transcript start position is less than the gene start position.
"""
import argparse
from bioservices import BioMart
import io
import pandas as pd
from io import StringIO
import time
from bisect import bisect_left, bisect_right

def findIntersectingEnsgs(pipeline_res_df, chromosome, gff_file):

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
        print(ensembl_df)
        
        start_positions = ensembl_df["start"].astype(int).tolist()
        end_positions = ensembl_df["end"].astype(int).tolist()

        for index, row in chromosome_df.iterrows():
            strand = "1" if row["Strand"] == "forward" else "-1"
            current_ensg = row["ensembl_gene_id"]
            polyA_start = int(row["PolyA_Start"])

            if strand == "1":
                start_idx = bisect_left(start_positions, polyA_start)
                
                matching_indices = []
                for i in range(start_idx):
                    if start_positions[i] <= polyA_start <= end_positions[i]:
                        matching_indices.append(i)

                # Process all matching indices
                for matching_index in matching_indices:
                    result_row = ensembl_df.iloc[matching_index]
                    start_position = int(result_row["start_position"])
                    end_position = int(result_row["end_position"])

                    if (str(result_row["chromosome_name"]) == str(chromosome) and
                        str(result_row["strand"]) == str(strand) and
                        polyA_start >= start_position and polyA_start <= end_position):
                        print("MATCH")
                        if result_row["ensembl_gene_id"] != current_ensg:
                            matching_genes.append((result_row["ensembl_gene_id"], current_ensg, chromosome, strand, start_position, end_position))
            elif strand == "-1":
                start_idx = bisect_left(start_positions, polyA_start)
                print("reverse positions", polyA_start)
                # Iterate through potential matches to find the correct one
                matching_indices = []
                for i in range(start_idx):
                    if start_positions[i] <= polyA_start <= end_positions[i]:
                        matching_indices.append(i)

                # Process all matching indices
                for matching_index in matching_indices:
                    result_row = ensembl_df.iloc[matching_index]
                    start_position = int(result_row["start_position"])
                    end_position = int(result_row["end_position"])

                    if (str(result_row["chromosome_name"]) == str(chromosome) and
                        #str(result_row["strand"]) == str(strand) and
                        polyA_start >= start_position and polyA_start <= end_position):
                        print("MATCH")
                        if result_row["ensembl_gene_id"] != current_ensg:
                            matching_genes.append((result_row["ensembl_gene_id"], current_ensg, chromosome, strand, start_position, end_position))
    return pd.DataFrame(matching_genes, columns=["Matching_ENSG", "Current_ENSG", "Chromosome", "Strand", "Start_Position", "End_Position"])


def main(pipeline_input, gff_input, output_file):
    column_names = [
        "seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"
    ]

    # Read the GFF file with the specified column names
    gff_file = pd.read_csv(gff_input, sep="\t", comment='#', names=column_names)
    filter_file = pd.read_csv(pipeline_input, sep='\t', low_memory=False, header=None)
    # Convert the result to a DataFrame
    filter_file = pd.DataFrame(filter_file, columns=["Chromosome", "Source", "Type", "Start", "End", "Score", "Strand", "Phase", "Attributes", "ENSG", "PolyA_Start", "PolyA_End", "Fantom_Start", "Fantom_End", "Fantom_Name"])
    # Save the DataFrame to a CSV file
    result = findIntersectingEnsgs(filter_file, 1, gff_file)
    result.to_csv(output_file, index=False)
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Find intersecting ENSGs.')
    parser.add_argument('pipeline_input', help='Path to the input CSV file')
    parser.add_argument('gff_input', help='Path to the GFF input file')
    parser.add_argument('output_file', help='Path to the output CSV file')
    args = parser.parse_args()
    main(args.pipeline_input, args.gff_input, args.output_file)


