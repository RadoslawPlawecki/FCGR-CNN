"""
@author: Radosław Pławecki
"""

import re
import os
import pandas as pd
from tqdm import tqdm

class NsHandler:
    def __init__(self, sequence: str):
        self.sequence = sequence
        self.n = len(sequence)
        self.start = 0
        self.end = self.n

    def count_leading_ns(self) -> int:
        start = 0
        while start < self.n and self.sequence[start] == 'N':
            start += 1
        self.start = start
        return start

    def count_trailing_ns(self) -> int:
        end = self.n - 1
        while end >= self.start and self.sequence[end] == 'N':
            end -= 1
        self.end = end + 1
        return self.n - 1 - end

    def count_internal_ns(self) -> int:
        return self.sequence[self.start:self.end].count('N')

    def get_core_sequence(self) -> str | None:
        if self.start <= self.end:
            return self.sequence[self.start:self.end]
        return None


data_path = "data/genome/chr/raw/"
out_path = "data/genome/02_info_ns_check.csv"

first = True

for chr_path in tqdm(os.listdir(data_path)):
    with open(os.path.join(data_path, chr_path)) as f:
        seq = f.read()
    ns_handler = NsHandler(sequence=seq)
    df = pd.DataFrame({
        "accession": [os.path.splitext(chr_path)[0]],
        "leading_ns": [ns_handler.count_leading_ns()],
        "core_seq": [ns_handler.get_core_sequence()],
        "internal_ns": [ns_handler.count_internal_ns()],
        "trailing_ns": [ns_handler.count_trailing_ns()]
    })
    df.to_csv(out_path, mode="w" if first else "a", header=first, index=False, sep=";")
    first = False
