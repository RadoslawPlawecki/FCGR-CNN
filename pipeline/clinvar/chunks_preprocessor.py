"""
@author: Radosław Pławecki
"""

import os
import pandas as pd
from validators import grch38loc_validator

class ChunksPreprocessor:
    def __init__(self, dir_path):
        if not os.path.isdir(dir_path):
            raise ValueError("The path should be a directory with data chunks!")
        self.dir_path = dir_path

    def _load_chunk(self, file_path):
        return pd.read_csv(file_path, delimiter=';', usecols=["GRCh38Chromosome","GRCh38Location", "label", "GenomicReference"])

    def _clean_chunk(self, df):
        df['GenomicReference'] = df['GenomicReference'].astype(str)
        df = df[df["GenomicReference"].notna() & (df["GenomicReference"] != "")]
        ref_pos = df['GenomicReference'].str.extract(r'(?:g|m)\.(\d+)')[0]
        ref_pos = pd.to_numeric(ref_pos, errors='coerce').astype('Int64')
        loc = pd.to_numeric(df["GRCh38Location"], errors='coerce').astype('Int64')
        return df.loc[loc == ref_pos]

    def preprocess_chunks(self, output_path=None):
        first = True
        for chunk_path in os.listdir(self.dir_path):
            file_path = os.path.join(self.dir_path, chunk_path)
            if not os.path.isfile(file_path):
                continue
            df_chunk = self._load_chunk(file_path)
            df_chunk = self._clean_chunk(df_chunk)

            if output_path:
                df_chunk.to_csv(output_path, mode='w' if first else 'a', header=first, index=False, sep=';')

            first = False
        print(f"[INFO] Data was exported! Running validators...")
        grch38loc_validator(filename=output_path)
