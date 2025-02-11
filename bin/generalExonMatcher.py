#!/usr/bin/env python3

import pandas as pd
import re
import sys
import os
import subprocess
import argparse
import csv

def validate_gff(file_path):
    result = subprocess.run(['gffread', file_path, '-E'], capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Validation failed for {file_path}:\n{result.stderr}")
    else:
        print(f"{file_path} is valid.")

def preprocess_gff(file_path):
    processed_file = 'processed_' + os.path.basename(file_path)
    with open(file_path, 'r') as infile, open(processed_file, 'w') as outfile:
        for line in infile:
            if not line.startswith('#'):
                outfile.write(line)
    return processed_file

def extract_transcript_id_and_exon_number(df):
    if 'Name=' in df.iloc[0, 8]:
        df['transcript_id'] = df[8].str.extract('Name="(.*?)\..*?"')
        df['block_num'] = df[8].str.extract('Name=".*?_block(.*?)"')
    else:
        df['transcript_id'] = df[8].str.extract('transcript_id "(.*?)"')
        df['block_num'] = df[8].str.extract('exon_number "(.*?)"')
    return df

def match_exons_with_blocks(human_df, block_df, single_exon):
    # Check if 'chr' is present in any of the entries in the column
    if block_df[0].str.contains('chr').any():
        block_df[0] = block_df[0].str.replace('^chr', '', regex=True)
    print(block_df.head())
    human_df['transcript_id'] = human_df['Attributes'].str.extract('Parent=transcript:(.*?);')
    block_df = extract_transcript_id_and_exon_number(block_df)

    block_df['block_num'] = pd.to_numeric(block_df['block_num'], errors='coerce')
    block_df = block_df.dropna(subset=['block_num'])
    # FILTER OUT SINGLE EXON GENES
    if single_exon == True:
        max_block_num = block_df.groupby('transcript_id')['block_num'].transform('max')
        block_df = block_df[max_block_num > 1]
    blocks_forward = block_df[block_df[6] == '+']
    blocks_reverse = block_df[block_df[6] == '-']
    last_blocks_forward = blocks_forward.loc[blocks_forward.groupby('transcript_id')['block_num'].idxmax()]
    last_blocks_reverse = blocks_reverse.loc[blocks_reverse.groupby('transcript_id')['block_num'].idxmin()]

    last_blocks_forward[0] = last_blocks_forward[0].astype(str)
    last_blocks_reverse[0] = last_blocks_reverse[0].astype(str)

    matched_human_exons_forward = human_df[human_df['Strand'] == '+']
    matched_human_exons_forward = matched_human_exons_forward[matched_human_exons_forward['Start'].isin(last_blocks_forward[3])]
    matched_blocks_forward = last_blocks_forward[last_blocks_forward[3].isin(matched_human_exons_forward['Start'])]

    matched_human_exons_reverse = human_df[human_df['Strand'] == '-']
    matched_human_exons_reverse = matched_human_exons_reverse[matched_human_exons_reverse['End'].isin(last_blocks_reverse[4])]
    matched_blocks_reverse = last_blocks_reverse[last_blocks_reverse[4].isin(matched_human_exons_reverse['End'])]

    matched_human_exons_forward = matched_human_exons_forward.sort_values(by=['Start'])
    matched_human_exons_reverse = matched_human_exons_reverse.sort_values(by=['End'])
    matched_blocks_forward = matched_blocks_forward.sort_values(by=[3])
    matched_blocks_reverse = matched_blocks_reverse.sort_values(by=[4])

    matched_human_exons = pd.concat([matched_human_exons_forward, matched_human_exons_reverse])
    matched_blocks = pd.concat([matched_blocks_forward, matched_blocks_reverse])

    return matched_human_exons, matched_blocks

def filter_by_ensg(df1, df2):
    ensg_set = set(df2.iloc[:, -2])
    filtered_df = df1[df1.iloc[:, -2].isin(ensg_set)]
    return filtered_df

def write_file(file_path, data):
    data.to_csv(file_path, sep='\t', index=False, header=False)

def main():
    if len(sys.argv) < 4:
        print("Usage: python exonMatcher.py <human_file.gff> <fantom_file.gff> <additional_file.gff> <single_exon?> <output_directory>")
        sys.exit(1)

    human_file = sys.argv[1]
    fantom_file = sys.argv[2]
    additional_file = sys.argv[3]
    single_exon = sys.argv[4].lower() == 'true' if len(sys.argv) > 4 else True
    output_dir = sys.argv[5] if len(sys.argv) > 5 else os.getcwd()

    validate_gff(human_file)
    validate_gff(fantom_file)
    validate_gff(additional_file)

    processed_human_file = preprocess_gff(human_file)
    processed_fantom_file = preprocess_gff(fantom_file)
    processed_additional_file = preprocess_gff(additional_file)

    human_df = pd.read_csv(processed_human_file, sep='\t', on_bad_lines='skip')
    fantom_df = pd.read_csv(processed_fantom_file, sep='\t', header=None, comment='#', on_bad_lines='skip')
    additional_df = pd.read_csv(processed_additional_file, sep='\t', header=None, comment='#', on_bad_lines='skip')

    matched_human_exons_fantom, matched_fantom_blocks = match_exons_with_blocks(human_df, fantom_df, single_exon)
    matched_human_exons_additional, matched_additional_blocks = match_exons_with_blocks(human_df, additional_df, single_exon)

    output_dir = os.path.join(output_dir, '')
    matched_human_exons_fantom.to_csv(output_dir + 'matched_human_exons_fantom.gff', sep='\t', index=False, header=False)
    matched_fantom_blocks.to_csv(output_dir + 'matched_fantom_blocks.gff', sep='\t', index=False, header=False)
    matched_human_exons_additional.to_csv(output_dir + 'matched_human_exons_additional.gff', sep='\t', index=False, header=False)
    matched_additional_blocks.to_csv(output_dir + 'matched_encode_blocks.gff', sep='\t', index=False, header=False)

    # ENSG filtering
    
    filtered_human_exons_fantom = filter_by_ensg(matched_human_exons_fantom, matched_human_exons_additional)
    write_file(output_dir + 'filtered_matched_human_exons.gff' , filtered_human_exons_fantom)
    os.remove(processed_human_file)
    os.remove(processed_fantom_file)
    os.remove(processed_additional_file)

if __name__ == "__main__":
    main()
