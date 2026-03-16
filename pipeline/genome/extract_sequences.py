"""
@author: Radosław Pławecki
"""

import pandas as pd
from tqdm import tqdm
from pipeline.utils import get_seq
from pipeline.clinvar.validators import base_validator


#TODO: Add: function to mutate a sequence, extraction to CSV

def sequence_extractor(accession: str, loc: int, seq_len: int) -> str | None:
    seq = get_seq(accession)

    if loc == -1:
        return None

    half = seq_len // 2
    start = max(0, loc - half)
    end = min(len(seq), loc + half)

    return seq[start:end]

data = pd.read_csv("data/clinvar/04_updated_loc.csv", delimiter=';', usecols=['accession', 'loc', 'ref'])
seq_len = 500

matches, mismatches = 0, 0

print(f"[INFO] Extracting sequences of length {seq_len}")

for row in tqdm(data.itertuples(index=False), total=len(data)):
    ref_sequence = sequence_extractor(row.accession, int(row.loc), seq_len=seq_len)
    if ref_sequence is None:
        continue
    if base_validator(sequence=ref_sequence, loc=seq_len // 2, ref=row.ref):
        matches += 1
    else:
        mismatches += 1

match_rate = matches / (matches + mismatches) * 100

print("[INFO] Validation summary")
print(f"Matches: {matches}")
print(f"Mismatches: {mismatches}")
print(f"Success rate: {match_rate:.2f}%")
