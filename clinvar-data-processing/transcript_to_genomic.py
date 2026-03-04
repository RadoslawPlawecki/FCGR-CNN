"""
@author: Radoslaw Plawecki
"""

import pandas as pd
import re

from tqdm import tqdm

from hgvs.dataproviders.uta import connect
import hgvs.parser
import hgvs.variantmapper

# import data and convert it to DataFrame
input_filename = "cleaned_data_clinvar.csv"
print(f"=== [1] Loading file {input_filename} ===")

data = pd.read_csv(f"data/{input_filename}", delimiter=';')
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

g_vars = []

with tqdm(total=len(cleaned_names), desc="Mapping c -> g") as pbar:
    for i, cleaned_name in enumerate(cleaned_names):
        print(f"\n[4.{i}] Variant processing: {cleaned_name}")
        try:
            if cleaned_name.startswith("NC_"):
                g_vars.append(cleaned_name)
                print(f"   -> No change: {cleaned_name}")
            else:
                c_var = hp.parse_hgvs_variant(cleaned_name)
                print(f"   -> Parsed: {c_var}")

                # download genomic accession
                tx_ac = c_var.ac
                print(f"   -> Transcript: {tx_ac}")

                # download mappings
                mappings = hdp.get_tx_mapping_options(tx_ac)
                print(f"   -> Number of available mappings: {len(mappings)}")

                alt_ac = mappings[1][1]
                print(f"   -> Chosen genomic accession: {alt_ac}")

                # conversion
                g_var = vm.c_to_g(c_var, alt_ac)
                print(f"   -> Genomic result: {g_var}")
                g_vars.append(g_var)

        except Exception as e:
            print(f"   !! Found an error while processing {cleaned_name}: {e}")

        pbar.update(1)

output_filename = "final_data_clinvar.csv"
print(f"=== [5] Saving to {output_filename} ===")

df["GenomicReference"] = g_vars
df.to_csv(f"data/{output_filename}", sep=';', index=False)

print("\n=== PROGRAM FINISHED ===")