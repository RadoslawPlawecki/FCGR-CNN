"""
@author: Radosław Pławecki
"""

import pandas as pd
from utils import get_seq
from tqdm import tqdm
from pipeline.clinvar.validators import base_validator


def sequence_extractor(accession: str, loc: int, seq_len: int) -> str | None:
    seq = get_seq(accession)

    if loc == -1:
        return None

    half = seq_len // 2
    start = max(0, loc - half)
    end = min(len(seq), loc + half)

    return seq[start:end]


data = pd.read_csv("data/clinvar/05_updated_loc_sample.csv", delimiter=';', usecols=['accession', 'loc', 'ref'])
seq_len = 500

print(len(data))

"""print(f"[INFO] Extracting sequences of length {seq_len}")
for row in tqdm(data.itertuples(index=False), total=len(data)):
    ref_sequence = sequence_extractor(row.accession, int(row.loc), seq_len=seq_len)
    if ref_sequence is None:
        print(ref_sequence)
        continue
    if base_validator(row.accession, int(row.loc), row.ref, seq=ref_sequence):
        print("[INFO] Match!")"""
