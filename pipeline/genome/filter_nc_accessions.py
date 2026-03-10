"""
@author: Radosław Pławecki, Filip Mirski
"""

from Bio import SeqIO
import pandas as pd
import os

if not os.path.isdir("data/genome/chr/"):
    os.makedirs("data/genome/chr/", exist_ok=True)

for record in SeqIO.parse("data/genome/01_GCF_000001405.26_GRCh38_genomic.fna", "fasta"):
    if record.id.startswith("NC_"):
        print(f"[INFO] Processing {record.id}...")
        with open(f"data/genome/chr/{record.id}.txt", 'w') as f:
            f.write(str(record.seq).upper())
