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
        