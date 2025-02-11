import pandas as pd
from Overlap import findIntersectingEnsgs

import pandas as pd

pipeline_data = {
    "Chromosome": ["1", "1", "2", "2", "3", "3"],
    "Source": ["source1", "source2", "source3", "source4", "source5", "source6"],
    "Type": ["type1", "type2", "type3", "type4", "type5", "type6"],
    "start_position": [100, 200, 300, 400, 500, 600],
    "end_position": [150, 250, 350, 450, 550, 650],
    "Score": [0, 0, 0, 0, 0, 0],
    "Strand": ["forward", "reverse", "forward", "reverse", "forward", "reverse"],
    "Phase": [".", ".", ".", ".", ".", "."],
    "Attributes": ["attr1", "attr2", "attr3", "attr4", "attr5", "attr6"],
    "ensembl_gene_id": ["ENSG000001", "ENSG000002", "ENSG000003", "ENSG000004", "ENSG000005", "ENSG000006"],
    "PolyA_Start": [120, 220, 320, 420, 520, 620],
    "PolyA_End": [130, 230, 330, 430, 530, 630],
    "Fantom_Start": [110, 210, 310, 410, 510, 610],
    "Fantom_End": [140, 240, 340, 440, 540, 640],
    "Fantom_Name": ["Fantom1", "Fantom2", "Fantom3", "Fantom4", "Fantom5", "Fantom6"]
}

pipeline_df = pd.DataFrame(pipeline_data)

# Create a test GFF DataFrame
gff_data = {
    "chromosome_name": ["1", "1", "2", "2", "3", "3"],
    "source": ["ensembl", "ensembl", "ensembl", "ensembl", "ensembl", "ensembl"],
    "feature": ["gene", "gene", "gene", "gene", "gene", "gene"],
    "start": [90, 190, 290, 390, 490, 590],
    "end": [121, 260, 360, 460, 560, 660],
    "start_position": [90, 190, 290, 390, 490, 590],
    "end_position": [160, 260, 360, 460, 560, 660],
    "score": [".", ".", ".", ".", ".", "."],
    "strand": ["1", "-1", "1", "-1", "1", "-1"],
    "frame": [".", ".", ".", ".", ".", "."],
    "attribute": [
        "gene_id \"ENSG000001\";", "gene_id \"ENSG000002\";", 
        "gene_id \"ENSG000003\";", "gene_id \"ENSG000004\";", 
        "gene_id \"ENSG000005\";", "gene_id \"ENSG000006\";"
    ],
    "ensembl_gene_id": [
        "ENSG000002","ENSG000003", 
        "ENSG000004", "ENSG000005", 
        "ENSG000006", "ENSG000007"
    ]
}


gff_df = pd.DataFrame(gff_data)


test = findIntersectingEnsgs(pipeline_df, gff_df)
print(test)