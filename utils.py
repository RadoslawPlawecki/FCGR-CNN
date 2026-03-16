"""
@author: Radosław Pławecki
"""

import re

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
        