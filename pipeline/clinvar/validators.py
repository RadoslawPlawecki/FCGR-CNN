"""
@author: Radosław Pławecki
"""

import pandas as pd
from tqdm import tqdm


def grch38loc_validator(filename):
    df = pd.read_csv(filename, delimiter=';')

    # extract positions for g. or m. references
    df['RefPos'] = df['GenomicReference'].str.extract(r'(?:g|m)\.(\d+)')

    df['GRCh38Location_int'] = pd.to_numeric(df['GRCh38Location'], errors='coerce').astype('Int64')
    df['RefPos_int'] = pd.to_numeric(df['RefPos'], errors='coerce').astype('Int64')

    # compare
    df['Match'] = df['GRCh38Location_int'] == df['RefPos_int']

    # validation summary
    total_rows = len(df)
    matches = df['Match'].sum()
    mismatches = total_rows - matches
    match_rate = matches / total_rows * 100

    print("=== Validation Summary [GRCh38 Location] ===")
    print(f"Total rows scanned: {total_rows}")
    print(f"Successful matches: {matches}")
    print(f"Mismatches / errors: {mismatches}")
    print(f"Success rate: {match_rate:.2f}%")

    if mismatches > 0:
        print("\n--- Example mismatches ---")
        print(df.loc[~df['Match'], ['GRCh38Location', 'GenomicReference', 'RefPos']].head(10))
        

genome_cache = {}


def get_seq(accession):
    if accession not in genome_cache:
        with open(f"data/genome/chr/preprocessed/{accession}.txt") as f:
            genome_cache[accession] = f.read()
    return genome_cache[accession]


def base_validator(accession: str, loc: int, ref: str) -> bool:
    seq = get_seq(accession)
    return seq[loc] == ref



data = pd.read_csv("data/clinvar/04_updated_loc.csv", delimiter=';', usecols=['accession', 'loc', 'ref'])

print("[INFO] Running validator")

matches, mismatches = 0, 0

for row in tqdm(data.itertuples(index=False), total=len(data)):
    if base_validator(row.accession, int(row.loc), row.ref):
        matches += 1
    else:
        mismatches += 1

match_rate = matches / (matches + mismatches) * 100

print("[INFO] Validation summary")
print(f"Matches: {matches}")
print(f"Mismatches: {mismatches}")
print(f"Success rate: {match_rate:.2f}%")
