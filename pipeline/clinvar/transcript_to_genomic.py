"""
@author: Radosław Pławecki
"""

import pandas as pd
import numpy as np
import logging
import re

from tqdm import tqdm
from multiprocessing import Pool, cpu_count, TimeoutError

from hgvs.dataproviders.uta import connect
import hgvs.parser
import hgvs.variantmapper

logging.getLogger("bioutils.seqfetcher").setLevel(logging.ERROR)

input_filename = "data/clinvar/cleaned_data_clinvar.csv"

print(f"=== Loading file {input_filename} ===")
data = pd.read_csv(f"{input_filename}", delimiter=';')
df = pd.DataFrame(data)

print(f"[INFO] {len(df)} data rows loaded")
hp = None
hdp = None
vm = None


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
        # cache mappings
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


chunk_size = int(len(df) * 0.005)  # 0.5%
if chunk_size == 0:
    chunk_size = 1

n_chunks = int(np.ceil(len(df) / chunk_size))

print(f"[INFO] Data divided into {n_chunks} chunks (~0.5% each)")
for chunk_id in range(n_chunks):
    start = chunk_id * chunk_size
    end = start + chunk_size
    chunk = df.iloc[start:end].copy()

    print(f"\n=== Processing chunk {chunk_id+1}/{n_chunks} ===")
    names = chunk["Name"]
    cleaned_names = []

    print("=== Cleaning variants names ===")
    for name in names:
        cleaned_name = re.sub(r"[\[\(\{][^\]\)\}]*[\]\)\}]", "", name)
        cleaned_names.append(cleaned_name.strip())
    mapping_cache = {}

    print("=== Mapping c -> g (parallel) ===")
    workers = max(cpu_count() - 1, 1)
    g_vars = []
    skipped = 0
    with Pool(workers, initializer=init_worker) as pool:
        async_results = [
            pool.apply_async(process_variant, (name,))
            for name in cleaned_names
        ]
        for r in tqdm(async_results, total=len(async_results)):
            try:
                g_vars.append(r.get(timeout=20))  # 20 second timeout
            except TimeoutError:
                skipped += 1
                (f"Timeout Error: Record was skipped. Number of skipped records: {skipped}.")
                g_vars.append(np.nan)
    chunk["GenomicReference"] = g_vars
    output_filename = f"data/clinvar/chunks/CH_{chunk_id+1}.csv"

    print(f"=== Saving chunk to {output_filename} ===")
    chunk.to_csv(output_filename,sep=";",index=False)

print("=== PROGRAM FINISHED ===")
