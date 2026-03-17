"""
@author: Radosław Pławecki
"""

import re
from pathlib import Path
import polars as pl


def is_valid_dna(sequence):
        """
        Function to check if a DNA sequence is valid.

        :param sequence: DNA sequence.
        :return bool: True, if a DNA sequence is valid, otherwise - False.
        """
        return bool(re.fullmatch(r"[ACGTacgt]+", sequence))
        

genome_cache = {}


def get_seq(accession):
    if accession not in genome_cache:
        with open(f"data/genome/chr/preprocessed/{accession}.txt") as f:
            genome_cache[accession] = f.read()
    return genome_cache[accession]

def load_fcgr_seq(path: Path) -> pl.DataFrame:
    try:
        df = pl.read_csv(path, separator=";")
    except FileNotFoundError:
        print(f"File not found: {path}")
        raise
    except Exception as e:
        print(f"Error loading FCGR sequence data: {e}")
        raise
    print(f"[INFO] Successfully loaded FCGR sequence data from {path}")
    return df
