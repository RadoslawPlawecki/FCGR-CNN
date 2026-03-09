"""
@author: Radosław Pławecki
"""

from Bio import SeqIO

import pandas as pd

records = list(SeqIO.parse("data/genome/01_GCF_000001405.26_GRCh38_genomic.fna", "fasta"))

accessions, sequences = [], []
for i, record in enumerate(records):
    if record.id.startswith("NC_"):
        accessions.append(record.id)
        sequences.append(record.seq)

data = {
    "GenAcc": accessions,
    "RefSeq": sequences
}

df = pd.DataFrame(data)
df.to_csv("data/genome/02_genome_nc_only.csv", sep=';', index=True)

