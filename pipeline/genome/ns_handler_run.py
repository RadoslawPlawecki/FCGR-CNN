"""
@author: Radosław Pławecki
"""

import os
import pandas as pd
from ns_handler import NsHandler
from tqdm import tqdm

data_path = "data/genome/chr/raw/"
out_path = "data/genome/chr/preprocessed/"

first = True

for chr_path in tqdm(os.listdir(data_path)):
    with open(os.path.join(data_path, chr_path)) as f:
        seq = f.read()
    ns_handler = NsHandler(sequence=seq)
    df = pd.DataFrame({
        "accession": [os.path.splitext(chr_path)[0]],
        "leading_ns": [ns_handler.count_leading_ns()],
        "core_seq": [len(ns_handler.get_core_sequence())],
        "internal_ns": [ns_handler.count_internal_ns()],
        "trailing_ns": [ns_handler.count_trailing_ns()]
    })
    df.to_csv("data/genome/02_info_ns.csv", mode="w" if first else "a", header=first, index=False, sep=";")
    first = False
    with open(os.path.join(out_path, chr_path), 'w') as f:
        f.write(ns_handler.get_core_sequence())
