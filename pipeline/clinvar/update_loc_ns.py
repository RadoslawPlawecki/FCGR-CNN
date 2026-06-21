"""
@author: Radosław Pławecki
"""

import os
import pandas as pd


ns_path = "data/genome/02_info_ns.csv"
var_path = "data/clinvar/03_preprocessed_clinvar.csv"

out_path = "data/clinvar/04_updated_loc.csv"

df_ns = pd.read_csv(ns_path, delimiter=';', usecols=['accession', 'leading_ns'])
df_var = pd.read_csv(var_path, delimiter=';')

df = df_var.merge(df_ns[['accession', 'leading_ns']], on='accession', how='left')
# subtract leading N's and 1 (nucleotides are counted from 1)
df['loc'] = df['loc'] - df['leading_ns'] - 1

df = df.dropna(subset=['loc'])
df = df.drop(columns='leading_ns')

df.to_csv(out_path, sep=';', index=False)
