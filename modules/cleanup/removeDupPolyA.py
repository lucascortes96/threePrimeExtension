import pandas as pd
import argparse

# Function to filter the DataFrame
def filter_group(group):
    if group.iloc[0]['Strand'] == '+':
        return group.loc[group['PolyA_Start'].idxmax()]
    else:
        return group.loc[group['PolyA_Start'].idxmin()]

def main(input_file, output_file):
    # Specify the data types for the columns
    dtype_dict = {
        'Column1': 'str',  # Replace 'Column1' with the actual column name
        'Column2': 'str',  # Replace 'Column2' with the actual column name
        'PolyA_Start': 'str',  # Read as string initially to handle conversion
        'Strand': 'str',
        'Attributes': 'str'
        # Add other columns as needed
    }

    # Read the data into a DataFrame with specified dtypes
    df = pd.read_csv(input_file, sep='\t', dtype=dtype_dict, low_memory=False)

    # Clean the PolyA_Start column to remove non-numeric values
    df['PolyA_Start'] = pd.to_numeric(df['PolyA_Start'], errors='coerce')
    df = df.dropna(subset=['PolyA_Start'])

    # Group by 'ENSG' and apply the filter function
    filtered_df = df.groupby('ENSG').apply(filter_group).reset_index(drop=True)

    # Save the filtered DataFrame to a new CSV file
    filtered_df.to_csv(output_file, sep='\t', index=False, header=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Remove duplicate PolyA entries from a CSV file.')
    parser.add_argument('input_file', help='Path to the input CSV file')
    parser.add_argument('output_file', help='Path to the output CSV file')
    args = parser.parse_args()
    
    main(args.input_file, args.output_file)