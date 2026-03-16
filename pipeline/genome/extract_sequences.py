"""
@author: Radosław Pławecki
"""

import pandas as pd
from tqdm import tqdm
from pipeline.genome.sequence_extractor import SequenceExtractor
from pipeline.clinvar.validators import base_validator

input_path = "data/clinvar/05_updated_loc_sample.csv"
output_path = "data/fcgr/01_fcgr_sequences.csv"

df = pd.read_csv(input_path, delimiter=';')

seq_lens = [500, 1000, 1500, 2000]

matches, mismatches = 0, 0

print("--- PROGRAM STARTED ---")
for seq_len in seq_lens:
    print(f"[INFO] Extracting sequences of length {seq_len}")
    ref_col = []
    mut_col = []
    for row in tqdm(df.itertuples(index=False), total=len(df)):
        ref, mut = SequenceExtractor(row.accession, row.loc, row.mut, seq_len).get_sequences()
        ref_col.append(ref)
        mut_col.append(mut)
        if base_validator(sequence=ref, loc=seq_len // 2, ref=row.ref):
            matches += 1
        else:
            mismatches += 1
        if base_validator(sequence=mut, loc=seq_len // 2, ref=row.mut):
            matches += 1
        else:
            mismatches += 1
    df[f"ref_{seq_len}"] = ref_col
    df[f"mut_{seq_len}"] = mut_col

match_rate = matches / (matches + mismatches) * 100

print("[INFO] Validation summary")
print(f"Matches: {matches}")
print(f"Mismatches: {mismatches}")
print(f"Success rate: {match_rate:.2f}%")

df.to_csv(output_path, sep=';', index=False)
