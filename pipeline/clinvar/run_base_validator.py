"""
@author: Radosław Pławecki
"""

import pandas as pd
from tqdm import tqdm
from pipeline.clinvar.validators import base_validator_acc


data = pd.read_csv("data/clinvar/04_updated_loc.csv", delimiter=';', usecols=['accession', 'loc', 'ref'])

print("[INFO] Running validator")

matches, mismatches = 0, 0

for row in tqdm(data.itertuples(index=False), total=len(data)):
    if base_validator_acc(row.accession, int(row.loc), row.ref):
        matches += 1
    else:
        mismatches += 1

match_rate = matches / (matches + mismatches) * 100

print("[INFO] Validation summary")
print(f"Matches: {matches}")
print(f"Mismatches: {mismatches}")
print(f"Success rate: {match_rate:.2f}%")