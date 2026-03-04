"""
@author: Radosław Pławecki
"""

import pandas as pd

filename = "final_data_clinvar.csv"

# load data
df = pd.read_csv(f"data/{filename}", delimiter=';')

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

print("=== Validation Summary ===")
print(f"Total rows scanned: {total_rows}")
print(f"Successful matches: {matches}")
print(f"Mismatches / errors: {mismatches}")
print(f"Success rate: {match_rate:.2f}%")

if mismatches > 0:
    print("\n--- Example mismatches ---")
    print(df.loc[~df['Match'], ['GRCh38Location', 'GenomicReference', 'RefPos']].head(10))
    