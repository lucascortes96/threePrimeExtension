import re
import sys

def load_stable_ids(stable_id_file):
    stable_ids = set()
    with open(stable_id_file, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith("chromosome"):
                stable_id = line.strip().split('\t')[-1]
                stable_ids.add(stable_id)
                #print(stable_ids)
    return stable_ids

def filter_protein_coding_genes(input_file, stable_id_file, output_file):
    stable_ids = load_stable_ids(stable_id_file)
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        gene_block = []
        keep_block = False

        for line in infile:
            if line.startswith("###"):
                if keep_block:
                    outfile.write("".join(gene_block))
                    outfile.write(line)
                gene_block = []
                keep_block = False
            else:
                gene_block.append(line)
                if "biotype=protein_coding" in line:
                    attributes = line.split('\t')[8]
                    ensg_list = [attr.split('=')[1].replace('gene:', '') for attr in attributes.split(';') if attr.startswith("ID=gene:")]
                    #print(ensg_list)
                    if ensg_list and ensg_list[0] not in stable_ids:
                        #print(ensg_list)
                        keep_block = True
                    else:
                        print(ensg_list)

        # Write the last block if it should be kept
        if keep_block:
            outfile.write("".join(gene_block))

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python proteinCodingGeneFilterGff.py <input_file> <stable_id_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    stable_id_file = sys.argv[2]
    output_file = sys.argv[3]
    filter_protein_coding_genes(input_file, stable_id_file, output_file)