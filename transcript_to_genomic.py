"""
@author: Radoslaw Plawecki
"""

import pandas as pd
import re

from hgvs.dataproviders.uta import connect
import hgvs.parser
import hgvs.variantmapper

# import data and convert it to DataFrame
filename = "cleaned_data_clinvar.csv"
print(f"=== [1] Loading file {filename} ===")

data = pd.read_csv(f"data/{filename}", delimiter=';')
df = pd.DataFrame(data)
print(f"[INFO] {len(df)} data rows loaded")

names = df['Name']
print(f"[INFO] Number of names to process: {len(names)}")

cleaned_names = []
print("=== [2] Cleaning variants names ===")

for i, name in enumerate(names):
    print(f"[2.{i}] Processing: {name}")
    cleaned_name = re.sub(r"[\[\(\{][^\]\)\}]*[\]\)\}]", "", name)
    cleaned_names.append(cleaned_name.strip())
    print(f"   -> {cleaned_name}")

print(f"[INFO] Number of names after cleaning: {len(cleaned_names)}")

# initialize parser, data provider and variant mapper
print("=== [3] Parser initialization and UTA connection ===")
hp = hgvs.parser.Parser()
hdp = connect()
vm = hgvs.variantmapper.VariantMapper(hdp)
print("[INFO] Initialization finished")

print("=== [4] Start mapping c -> g ===")
for i, cleaned_name in enumerate(cleaned_names):
    print(f"\n[4.{i}] Variant processing: {cleaned_name}")
    try:
        if cleaned_name.startswith("NC_"):
            print(f"   -> No change: {cleaned_name}")
            continue

        c_var = hp.parse_hgvs_variant(cleaned_name)
        print(f"   -> Parsed: {c_var}")

        # pobranie akcesji transkryptu
        tx_ac = c_var.ac
        print(f"   -> Transcript: {tx_ac}")

        # pobranie mapowań
        mappings = hdp.get_tx_mapping_options(tx_ac)
        print(f"   -> Number of available mappings: {len(mappings)}")

        alt_ac = mappings[1][1]
        print(f"   -> Chosen genomic accession: {alt_ac}")

        # Konwersja
        g_var = vm.c_to_g(c_var, alt_ac)
        print(f"   -> Genomic result: {g_var}")

    except Exception as e:
        print(f"   !! Found an error while processing {cleaned_name}: {e}")

print("\n=== PROGRAM FINISHED ===")
