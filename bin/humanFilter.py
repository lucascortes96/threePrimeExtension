#!/usr/bin/env python3
import re
import sys
import pandas as pd

def load_readthrough_list(readthrough_file):
    readthrough_df = pd.read_csv(readthrough_file, sep='\t')
    return set(readthrough_df['stable_id'])

def filter_protein_coding_genes(input_file, output_file, readthrough_file, single_exon):
    readthrough_ids = load_readthrough_list(readthrough_file)

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        gene_block = []
        keep_block = False
        if single_exon == True:
            exon_count = 0
        else:
            exon_count = 1 

        for line in infile:
            if line.startswith("###"):
                if keep_block and exon_count > 1:
                    outfile.write("".join(gene_block))
                    outfile.write(line)
                gene_block = []
                keep_block = False
                exon_count = 0
                
            else:
                gene_block.append(line)
                if "biotype=protein_coding" in line:
                    keep_block = True
                if any(str(readthrough_id) in line for readthrough_id in readthrough_ids):
                    keep_block = False
                if "\texon\t" in line:
                    exon_count += 1
                    

        # Write the last block if it should be kept
        if keep_block and exon_count > 1:
            outfile.write("".join(gene_block))

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python proteinCodingGeneFilterGff.py <input_file> <output_file> <readthrough_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    readthrough_file = sys.argv[3]
    single_exon = sys.argv[4].lower() == 'true'

    if single_exon:
        print("Single exon is True")
    else:
        print("Single exon is False")
    filter_protein_coding_genes(input_file, output_file, readthrough_file, single_exon)