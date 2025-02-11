"""
This script is used to process BigWig files and CSV files containing genomic data.
It extracts intervals from the BigWig file, checks for downstream values, and calculates drop-offs.

Functions:
- getIntervals: Extracts intervals from the BigWig file.
- getDownstreamValue: Checks for downstream values in the intervals.
- getDropOff: Calculates drop-offs in the intervals.
- main: The main function that orchestrates the execution of the script.
"""

import pyBigWig
import sys 
import pandas as pd 
import os

def getIntervals(bw, comp, chrom):
    """
    Extracts intervals from the BigWig file and updates the DataFrame with a boolean indicating if any interval value is greater than 10.
    
    Parameters:
    bw (pyBigWig object): The BigWig file.
    comp (DataFrame): The DataFrame containing genomic data.
    chrom (str): The chromosome number.
    
    Returns:
    DataFrame: The updated DataFrame with a new column 'longreadSupport'.
    """
    if comp is not None:
        comp['longreadSupport'] = False  # Initialize the new column with False
        for index, row in comp.iterrows():
            chromosome = row.iloc[0]
            start = row.iloc[9]
            
            end = row.iloc[10]
           
            strand = row.iloc[6]
            
            if strand == 'forward':
                print(start, end)
                strand = '+'
                end = end 
                start = start 
                
            else:
                strand = '-'
                start = start 
                end = end 
                print("BACKWARDS")
            try:
                interval_data = bw.intervals(chrom, start, end)
                print("INTTERVAL")
                if not interval_data:
                    interval_data = [(start, end, 0)]  # Default value if no intervals found
            except RuntimeError:
                interval_data = [(start, end, 0)]  # Default value in case of error

            # Check if any value in interval_data is greater than 10
            has_value_greater_than_10 = any(value > 10 for _, _, value in interval_data)
            print(has_value_greater_than_10)
            comp.at[index, 'longreadSupport'] = has_value_greater_than_10
            print(comp)
        
    return comp

def getDownstreamValue(intervals, comp):
    """
    Checks for downstream values in the intervals.
    
    Parameters:
    intervals (list): A list of tuples containing the chromosome, strand, and intervals.
    
    Returns:
    list/bool: The intervals if a downstream value > 5 is found, False otherwise.
    """
    matchedIntervals = []
    for i in range(len(intervals)):
        #print(intervals[i][2][0][2])
        if intervals[i][3][0][2] > 10:
            print(intervals[i][3][0][2])
            #print(intervals[i][2][0])
            matchedIntervals.append(intervals[i])
            comp.loc[intervals[i][0], 'longReadMatched'] = True
            print("TRUE")
    return matchedIntervals, comp
        
    
    


def getDropOff(peakHeight, binSize, chrom, intervals):
    """
    Calculates drop-offs in the intervals.
    
    Parameters:
    peakHeight (float): The peak height.
    binSize (int): The bin size.
    chrom (str): The chromosome number.
    intervals (list): A list of tuples containing the chromosome, strand, and intervals.
    
    Returns:
    bool: True if a drop-off is found, False otherwise.
    """
    for i in range(len(intervals)):
        polyAStart = intervals[i][0]
        vals = intervals.values(chrom, polyAStart, polyAStart + binSize)
        for i in range(len(vals)):
            if vals < 0.5*peakHeight:
                return True
            else:
                return False
    return False

def main():
    """
    The main function that orchestrates the execution of the script.
    """
    binSize = 10
    if len(sys.argv) < 2:
        print("Usage: python script.py <bigwig_file> <csv_file>")
        sys.exit(1)
    bw = sys.argv[1]
    chrom = bw.split('_')[-1].replace('.txt', '')
    chromosomes = list(range(1, 23)) + ['X', 'Y']
    bw = pyBigWig.open(bw)
    comp = None
    if len(sys.argv) > 2:
        comp = pd.read_csv(sys.argv[2], sep='\t')
        #comp['longReadMatched'] = False
        # Get the directory and base name of the input file
        input_dir, input_file = os.path.split(sys.argv[2])
        # Remove the extension from the base name
        base_name, _ = os.path.splitext(input_file)
    dsValues = []
    for i in chromosomes:
        chrom = 'chr' + str(i)
        comp = getIntervals(bw, comp, chrom)
        '''
        if intervals:
            print("intervals!")
            dsValuesNow, comp = getDownstreamValue(intervals, comp)
            for interval in dsValuesNow:
                dsValues.append([interval[0], interval[1], interval[2][0]])  # Appending as a list
    print("Length of matches",len(dsValues))
    
    # Assuming dsValues is a list of lists where each inner list represents a row
    df = pd.DataFrame(dsValues)

   # Create the output file names
    dsValues_file = os.path.join(input_dir, base_name + '_dsValues.csv')
    '''
    updated_comp_file = os.path.join(input_dir, base_name + '_updated.csv')

    # Save DataFrame to a CSV file
    #df.to_csv(dsValues_file, index=False)  
    comp.to_csv(updated_comp_file, index=False)

if __name__ == "__main__":
    main()