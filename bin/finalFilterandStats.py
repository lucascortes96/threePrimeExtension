#!/usr/bin/env python3
import pandas as pd
import argparse
import re
import csv
import statistics
import matplotlib.pyplot as plt
import os

# Function to filter the DataFrame
def filter_group(group):
    if group.iloc[0]['Strand'] == '+':
        return group.loc[group['PolyA_End'].idxmax()]
    else:
        return group.loc[group['PolyA_Start'].idxmin()]

def extract_transcript_name(attributes):
    match = re.search(r'Parent=transcript:(ENST\d+)', attributes)
    return match.group(1) if match else None

def main(input_file, output_file):
    # Specify the data types for the columns
    dtype_dict = {
        'Chromosome': 'str',
        'Source': 'str',
        'Type': 'str',
        'Start': 'str',
        'End': 'str',
        'Score': 'str',
        'Strand': 'str',
        'Phase': 'str',
        'Attributes': 'str',
        'ENSG': 'str',
        'PolyA_Start': 'str',  # Read as string initially to handle conversion
        'PolyA_End': 'str',    # Read as string initially to handle conversion
        'Transcript_Start': 'str',
        'Transcript_End': 'str'
    }

    # Read the data into a DataFrame with specified dtypes
    df = pd.read_csv(input_file, sep='\t', dtype=dtype_dict, low_memory=False)

    # Extract the Transcript_Name from the Attributes column
    df['Transcript_Name'] = df['Attributes'].apply(extract_transcript_name)

    # Clean the PolyA_Start and PolyA_End columns to remove non-numeric values
    df['PolyA_Start'] = pd.to_numeric(df['PolyA_Start'], errors='coerce')
    df['PolyA_End'] = pd.to_numeric(df['PolyA_End'], errors='coerce')
    df = df.dropna(subset=['PolyA_Start', 'PolyA_End', 'Transcript_Name'])

    # Group by 'Transcript_Name' and apply the filter function
    filtered_df = df.groupby('Transcript_Name').apply(filter_group).reset_index(drop=True)

    # Save the filtered DataFrame to a new CSV file with header
    filtered_df.to_csv(output_file, sep='\t', index=False, header=True)

    # Calculate differences and statistics
    differences = []
    ensg_largest_extension = None
    largest_difference = 0

    final_output_file = output_file.replace('.csv', '_final.csv')

    with open(output_file, mode='r') as infile, open(final_output_file, mode='w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')
        
        header = next(reader)  # Read header
        header.append('Difference')  # Add new column for differences
        writer.writerow(header)  # Write header to output file

        for row in reader:
            try:
                strand = row[6]
                if strand == '+':
                    selected_value = float(row[4])
                else:
                    selected_value = float(row[3])

                difference = abs(selected_value - float(row[10]))
                differences.append(difference)
                row.append(difference)  # Add difference to the row
                writer.writerow(row)  # Write row to output file

                if difference > largest_difference:
                    largest_difference = difference
                    ensg_largest_extension = row[9]

            except ValueError as e:
                print(f"Error converting row to float: {row} - {e}")
            except IndexError as e:
                print(f"Error with row indexing: {row} - {e}")

    # Calculate mean and median of the differences
    mean_difference = statistics.mean(differences)
    median_difference = statistics.median(differences)

    # Write statistics to a file
    stats_file = output_file.replace('.csv', '_threePrimeExtendStats.txt')
    plot_file = output_file.replace('.csv', '_threePrimeExtendPlot.png')
    with open(stats_file, 'w') as f:
        f.write(f"Mean of differences: {mean_difference}\n")
        f.write(f"Median of differences: {median_difference}\n")
        f.write(f"Largest extension: {largest_difference} (ENSG: {ensg_largest_extension})\n")

    # Create a distribution chart (histogram) of the differences
    plt.hist(differences, bins=30, edgecolor='black')
    plt.title('Distribution of Transcript Extensions 3 Prime')
    plt.xlabel('Transcript Extension Length')
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.xlim(left=0)  # Set the x-axis limit to start at 0
    plt.savefig(plot_file)
    plt.show()

    # Delete the intermediary output file
    os.remove(output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process PolyA entries from a CSV file.')
    parser.add_argument('input_file', help='Path to the input CSV file')
    parser.add_argument('output_file', help='Path to the output CSV file')
    args = parser.parse_args()
    
    main(args.input_file, args.output_file)