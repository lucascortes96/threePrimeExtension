#!/usr/bin/env python3
import pandas as pd
import sys
import os

# Get the input file name from the command line arguments
input_file = sys.argv[1]
identity = sys.argv[2]

# Construct the output file name
base_name, ext = os.path.splitext(input_file)
output_file = f"{identity}_nextRun.gff"


df = pd.read_csv(input_file, sep = '\t', usecols = ['Chromosome', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Phase', 'Attributes', 'ENSG'] )
print(df.head())



# Save the modified DataFrame to a new CSV file
df.to_csv(output_file, index=False, sep='\t')