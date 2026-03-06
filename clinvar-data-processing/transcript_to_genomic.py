"""
@author: Radosław Pławecki
"""

import pandas as pd
import numpy as np
import logging
import re

from tqdm import tqdm
from multiprocessing import Pool, cpu_count

from hgvs.dataproviders.uta import connect
import hgvs.parser
import hgvs.variantmapper

logging.getLogger("bioutils.seqfetcher").setLevel(logging.ERROR)

input_filename = "cleaned_data_clinvar.csv"
print(f"=== [1] Loading file {input_filename} ===")

data = pd.read_csv(f"clinvar-data-processing/data/{input_filename}", delimiter=';') 
df = pd.DataFrame(data) 
print(f"[INFO] {len(df)} data rows loaded")

names = df['Name']
print(f"[INFO] {len(names)} variants loaded")

cleaned_names = [] 
print("=== [2] Cleaning variants names ===")

for i, name in enumerate(names): 
    print(f"[2.{i}] Processing: {name}") 
    cleaned_name = re.sub(r"[\[\(\{][^\]\)\}]*[\]\)\}]", "", name) 
    cleaned_names.append(cleaned_name.strip()) 
    print(f" -> {cleaned_name}")

print(f"[INFO] Number of names after cleaning: {len(cleaned_names)}")

hp = None
hdp = None
vm = None

mapping_cache = {}

print("=== [3] Parser initialization and UTA connection ===")

def init_worker():
    global hp, hdp, vm
    hp = hgvs.parser.Parser()
    hdp = connect()
    vm = hgvs.variantmapper.VariantMapper(hdp)

def process_variant(cleaned_name):
    global hp, hdp, vm, mapping_cache
    try:
        if cleaned_name.startswith("NC_"):
            return cleaned_name
        c_var = hp.parse_hgvs_variant(cleaned_name)
        tx_ac = c_var.ac
        if tx_ac not in mapping_cache:
            mappings = hdp.get_tx_mapping_options(tx_ac)
            if len(mappings) == 0:
                return np.nan
            mapping_cache[tx_ac] = mappings[1][1]
        alt_ac = mapping_cache[tx_ac]
        g_var = vm.c_to_g(c_var, alt_ac)
        return str(g_var)
    except Exception:
        return np.nan

print("=== [4] Mapping c -> g (parallel) ===")

workers = max(cpu_count() - 1, 1)

with Pool(workers, initializer=init_worker) as pool:
    g_vars = list(
        tqdm(
            pool.imap(process_variant, cleaned_names),
            total=len(cleaned_names)
        )
    )

output_filename = "final_data_clinvar.csv"
print(f"=== [5] Saving to {output_filename} ===")

df["GenomicReference"] = g_vars
df.to_csv(f"clinvar-data-processing/data/{output_filename}", sep=";", index=False)

print("=== PROGRAM FINISHED ===")
