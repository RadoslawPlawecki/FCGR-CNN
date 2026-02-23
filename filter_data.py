"""
@author: Radoslaw Plawecki
"""

import pandas as pd

df = pd.read_csv("clinvar_result.txt", sep="\t", low_memory=False, index_col=0)

columns_to_keep = [
    "Gene(s)",
    "GRCh37Chromosome",
    "GRCh37Location",
    "GRCh38Chromosome",
    "GRCh38Location",
    "Variant type",
    "Molecular consequence",
    "Germline classification",
    "Germline review status",
    "VariationID"
]

df = df[columns_to_keep]

# only SNV
df = df[df["Variant type"] == "single nucleotide variant"]

# only missense
df = df[df["Molecular consequence"].str.contains("missense", na=False)]

# only clean classes
allowed_classes = [
    "Pathogenic",
    "Likely pathogenic",
    "Benign",
    "Likely benign",
]
df = df[df["Germline classification"].isin(allowed_classes)]

# only reliable reviews
allowed_review = [
    "practice guideline",
    "reviewed by expert panel",
    "criteria provided, multiple submitters, no conflicts"
]
df = df[df["Germline review status"].isin(allowed_review)]

# binary label
df["label"] = df["Germline classification"].apply(
    lambda x: 1 if "pathogenic" in x.lower() else 0
)

# if a variant affects several genes, use first gene
df["PrimaryGene"] = df["Gene(s)"].str.split("|").str[0]

# df.to_csv("cleaned_data_clinvar.csv", sep=';', index=True)
print(df["label"].value_counts())
