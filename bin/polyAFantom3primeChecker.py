#!/usr/bin/env python3
import pandas as pd
import sys

"""
Make sure you use the outputs taken from exonMatcher.py as input for this script. exonMatcher.py is a script
that compares the last exons of FANTOM transcripts and normal human transcripts to find matches. 

This script takes 5 arguments:
1. Path to the transcript file.
2. Path to the polyA file.
3. Path to the GFF file.
4. Path to the output file.
5. Chromosome value.


"""

def importGffs(human_file, fantom_file, encode_file, polyA_file):
    human_exons = pd.read_csv(human_file, sep='\t', header=None)
    fantom_transcripts = pd.read_csv(fantom_file, sep='\t', header=None).iloc[:, :9]
    polyA_sites = pd.read_csv(polyA_file, sep='\t', header=None, skiprows=1).iloc[:, :6]
    encode_transcripts = pd.read_csv(encode_file, sep='\t', header=None).iloc[:, :9]
    


    # Convert numeric chromosomes to int, leave 'X' and 'Y' as str
    def convert_chromosome(chromosome):
        try:
            return int(chromosome)
        except ValueError:
            return chromosome

    human_exons[0] = human_exons[0].apply(convert_chromosome)
    fantom_transcripts[0] = fantom_transcripts[0].apply(convert_chromosome)
    encode_transcripts[0] = encode_transcripts[0].apply(convert_chromosome)
    
    polyA_sites[0] = polyA_sites[0].apply(convert_chromosome)
    # Convert 'Start' and 'End' columns to integers
    human_exons[[3, 4]] = human_exons[[3, 4]].apply(pd.to_numeric)
    fantom_transcripts[[3, 4]] = fantom_transcripts[[3, 4]].apply(pd.to_numeric)
    encode_transcripts[[3, 4]] = encode_transcripts[[3, 4]].apply(pd.to_numeric)
    polyA_sites[[1, 2]] = polyA_sites[[1, 2]].apply(pd.to_numeric)

    # Add column names
    # Assign column names to the first few columns
    human_column_names = ['Chromosome', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Phase', 'Attributes', 'ENSG']
    human_exons.columns = human_column_names + list(human_exons.columns[len(human_column_names):])
    # Drop any extra columns if necessary
    human_exons = human_exons[human_column_names]

    fantom_column_names = ['Chromosome', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Phase', 'Attributes']
    fantom_transcripts.columns = fantom_column_names + list(fantom_transcripts.columns[len(fantom_column_names):])
    # Drop any extra columns if necessary
    fantom_transcripts = fantom_transcripts[fantom_column_names]
    encode_column_names = ['Chromosome', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Phase', 'Attributes']
    encode_transcripts.columns = encode_column_names + list(encode_transcripts.columns[len(encode_column_names):])
    # Drop any extra columns if necessary
    encode_transcripts = encode_transcripts[encode_column_names]
    
    polyA_sites.columns = ['Chromosome', 'chromStart', 'chromEnd', 'Name', 'Score', 'Strand']
    # Extract 'Name' from the 9th column
    fantom_transcripts['Name'] = fantom_transcripts['Attributes'].str.extract('Name="([^"]*)')
    encode_transcripts['Name'] = encode_transcripts['Attributes'].str.extract('gene_id\s+"([^"]+)"')
    print(human_exons)

    return human_exons, fantom_transcripts, encode_transcripts, polyA_sites

def findMatches(human_exons, fantom_transcripts, encode_transcripts, polyA_sites):
    # Initialize an empty dataframe to store the results
    results = pd.DataFrame()
    # Loop over human exons
    for i, exon in human_exons.iterrows():
        # Filter polyA sites based on the strand, chromosome, and position
        polyA_filtered = polyA_sites[
            (polyA_sites['Strand'] == exon['Strand']) &
            (polyA_sites['Chromosome'] == exon['Chromosome']) &
            ((polyA_sites['chromStart'] > exon['End']) if exon['Strand'] == '+' else (polyA_sites['chromStart'] < exon['Start'])) &
            ((polyA_sites['chromStart'] <= (exon['End'] + 10000)) if exon['Strand'] == '+' else (polyA_sites['chromStart'] >= (exon['Start'] - 10000)))
        ]

        # Loop over filtered polyA sites
        for j, polyA in polyA_filtered.iterrows():
            # Filter FANTOM transcripts based on the strand, chromosome, and position
            fantom_filtered = fantom_transcripts[
                (fantom_transcripts['Strand'] == exon['Strand']) &
                (fantom_transcripts['Chromosome'] == exon['Chromosome']) &
                ((fantom_transcripts['Start'] == exon['Start']) if exon['Strand'] == '+' else (fantom_transcripts['End'] == exon['End'])) &
                ((fantom_transcripts['End'] >= polyA['chromStart']) if exon['Strand'] == '+' else (fantom_transcripts['Start'] <= polyA['chromStart']))
            ]
            # Add the matching polyA chromStart to the fantom_filtered DataFrame
            fantom_filtered['polyA_chromStart'] = polyA['chromStart']

            # Filter ENCODE transcripts based on the same criteria as FANTOM transcripts
            encode_filtered = encode_transcripts[
                (encode_transcripts['Strand'] == exon['Strand']) &
                (encode_transcripts['Chromosome'] == exon['Chromosome']) &
                ((encode_transcripts['Start'] == exon['Start']) if exon['Strand'] == '+' else (encode_transcripts['End'] == exon['End'])) &
                ((encode_transcripts['End'] >= polyA['chromStart']) if exon['Strand'] == '+' else (encode_transcripts['Start'] <= polyA['chromStart']))
            ]
            # Add the matching polyA chromStart to the fantom_filtered DataFrame
            encode_filtered['polyA_chromStart'] = polyA['chromStart']
            # Append the filtered rows and drop duplicates
            combined_filtered = fantom_filtered[fantom_filtered['polyA_chromStart'].isin(encode_filtered['polyA_chromStart'])]


            print(f"PolyA site at index {j} has {len(combined_filtered)} matching transcripts")

            # If there are matching transcripts, add them to the results
            if not combined_filtered.empty:
                for k, transcript in combined_filtered.iterrows():
                    result = exon.copy()
                    result['PolyA_Start'] = polyA['chromStart']
                    result['PolyA_End'] = polyA['chromStart']
                    result['Transcript_Start'] = transcript['Start']
                    result['Transcript_End'] = transcript['End']
                    result['Transcript_Name'] = transcript.get('Name', 'Unknown')  # Assuming 'Name' column exists, otherwise default to 'Unknown'
                    results = pd.concat([results, pd.DataFrame([result])], ignore_index=True)
    return(results)

    # Save the results to a CSV file
def main(human_file, fantom_transcripts, encode_transcripts, polyA_sites, output_file, chromosome_value):
    imported = importGffs(human_file, fantom_transcripts, encode_transcripts, polyA_sites)
    human_exons = imported[0]
    fantom_transcripts = imported[1]
    encode_transcripts = imported[2]
    polyA_sites = imported[3]
    #print(polyA_sites['Chromosome'])
    #print(polyA_sites['Chromosome'].dtype)

    # Check if chromosome_value is 'X', 'Y', 'MT', or starts with 'G' or 'J'
    if chromosome_value in ['X', 'Y', 'MT'] or chromosome_value.startswith(('G', 'J')):
        polyA_sites['Chromosome'] = polyA_sites['Chromosome'].astype(str)
    else:
        try:
            chromosome_value = int(chromosome_value)
            polyA_sites['Chromosome'] = pd.to_numeric(polyA_sites['Chromosome'], errors='coerce')
        except ValueError:
            polyA_sites['Chromosome'] = polyA_sites['Chromosome'].astype(str)

    # Filter dataframes based on the chromosome value
    human_exons = human_exons[human_exons['Chromosome'] == chromosome_value]
    fantom_transcripts = fantom_transcripts[fantom_transcripts['Chromosome'] == chromosome_value]
    print(encode_transcripts)
    encode_transcripts = encode_transcripts[encode_transcripts['Chromosome'] == chromosome_value]
    polyA_sites = polyA_sites[polyA_sites['Chromosome'] == chromosome_value]
    print(polyA_sites)
    results = findMatches(human_exons, fantom_transcripts, encode_transcripts, polyA_sites)
    
    # Include the chromosome value in the output file name
    output_file = f"{output_file}_matched_chr{chromosome_value}.csv"
    results.to_csv(output_file, sep='\t', index=False)


if __name__ == "__main__":
    human_file = sys.argv[1]
    fantom_file = sys.argv[2]
    encode_file = sys.argv[3]
    polyA_file = sys.argv[4]
    output_file = sys.argv[5]
    chromosome_value = sys.argv[6]
    
    main(human_file, fantom_file, encode_file,  polyA_file, output_file, chromosome_value)