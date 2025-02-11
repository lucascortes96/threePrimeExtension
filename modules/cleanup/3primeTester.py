from ThreeprimeFilter import findIntersectingEnsgs

import pandas as pd
from io import StringIO
from bisect import bisect_left, bisect_right

# Create a test DataFrame
data = {
    "Chromosome": ["1", "1", "1", "1"],
    "Strand": ["forward", "forward", "reverse", "forward"],
    "ensembl_gene_id": ["ENSG000001980", "ENSG000001980", "ENSG00000188295", "ENSG00000188295"],
    "Start": [1709126, 1708476, 247101319, 247101319],
    "End": [1709477, 1708476, 247101100, 247101319],
    "PolyA_Start": [11122, 1709138, 14358, 247099027]
}

pipeline_res_df = pd.DataFrame(data)

'''
# Mock response from BioMart query
response = """ENSG00000290825\t1\t11121\t24894\t1
ENSG00000223972\t1\t12010\t13670\t1
ENSG00000237683\t1\t31109\t29554\t-1   
ENSG00068020\t1\t34554\t36081\t-1"""

ensembl_df = pd.read_csv(StringIO(response), sep="\t", 
                         names=["ensembl_gene_id", "chromosome_name", 
                                "start_position", "end_position", "strand"])

# Ensure numeric columns are of type int
ensembl_df["start_position"] = pd.to_numeric(ensembl_df["start_position"], errors='coerce')
ensembl_df["end_position"] = pd.to_numeric(ensembl_df["end_position"], errors='coerce')



pipeline_res_df = pd.DataFrame(data)

# Initialize a list to store results
matching_genes = []

for index, row in pipeline_res_df.iterrows():
    chromosome = row["Chromosome"]
    strand = "1" if row["Strand"] == "forward" else "-1"
    current_ensg = row["ensembl_gene_id"]
    # Define the range for this row
    polyA_start = row["PolyA_Start"]
    
    # Print the start_position and end_position for each row in ensembl_df
    for _, result_row in ensembl_df.iterrows():
        print(f"Result row start_position: {result_row['start_position']}, end_position: {result_row['end_position']}, polyA: {polyA_start}")
        
        # Filter the results for the current row's chromosome and range
        if (str(strand) == "1" and
            str(result_row["chromosome_name"]) == str(chromosome) and
            str(result_row["strand"]) == str(strand) and
            polyA_start >= result_row["start_position"] and
            polyA_start <= result_row["end_position"]):
            
            print(f"Row {index}:")
            print(f"Chromosome: {chromosome}, Strand: {strand}, PolyA_Start: {polyA_start}")
            print("Matching gene:")
            print(f"  Gene: {result_row['ensembl_gene_id']}, Start: {result_row['start_position']}, End: {result_row['end_position']}")
            # Append matching genes to the list
            #matching_genes.append(result_row["ensembl_gene_id"])
            if result_row["ensembl_gene_id"] != current_ensg:
                matching_genes.append((result_row["ensembl_gene_id"], current_ensg, chromosome, strand, result_row["start_position"], result_row["end_position"]))
            
        elif (str(strand) == "-1" and
            str(result_row["chromosome_name"]) == str(chromosome) and
            str(result_row["strand"]) == str(strand) and
            polyA_start <= result_row["start_position"] and
            polyA_start >= result_row["end_position"]):
            
            print(f"Row {index}:")
            print(f"Chromosome: {chromosome}, Strand: {strand}, PolyA_Start: {polyA_start}")
            print("Matching gene:")
            print(f"  Gene: {result_row['ensembl_gene_id']}, Start: {result_row['start_position']}, End: {result_row['end_position']}")
            
            # Append matching genes to the list
            if result_row["ensembl_gene_id"] != current_ensg:
                matching_genes.append((result_row["ensembl_gene_id"], current_ensg, chromosome, strand, result_row["start_position"], result_row["end_position"]))
        else:
            print(f"No match for row {index}:")
            print(f"Chromosome: {chromosome}, Strand: {strand}, PolyA_Start: {polyA_start}")
            print(f"Result row start_position: {result_row['start_position']}, end_position: {result_row['end_position']}")
# Print the final intersecting genes
print("Intersecting genes:")
print(matching_genes)
# Convert matching_genes to a DataFrame
matching_genes_df = pd.DataFrame(matching_genes, columns=["ref_gene_id", "result_ensg", "chromosome", "strand", "start_position", "end_position"])
print(matching_genes_df)
'''




result = findIntersectingEnsgs(pipeline_res_df)
print(result)
