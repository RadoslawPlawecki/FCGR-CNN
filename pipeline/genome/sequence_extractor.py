"""
@author: Radosław Pławecki
"""

from pipeline.utils import get_seq


class SequenceExtractor:
    """Class to handle sequence extraction from genome."""
    def __init__(self, accession: str, loc: int, mut: str, seq_len: int):
        self.seq = get_seq(accession)
        self.loc = int(loc)
        self.mut = mut
        self.seq_len = seq_len
        self.half_len = seq_len // 2

    def extract_sequence(self) -> str | None:
        """Method to extract sequences."""
        if self.loc == -1:
            return None
        start = max(0, self.loc - self.half_len)
        end = min(len(self.seq), self.loc + self.half_len)
        return self.seq[start:end]

    def mutate_sequence(self) -> str | None:
        """Method to introduce a missense mutation into a sequence."""
        seq = self.extract_sequence()
        if seq is None:
            return None
        mid_index = len(seq) // 2
        mutated_seq = seq[:mid_index] + self.mut + seq[mid_index + 1:]
        return mutated_seq

    def get_sequences(self) -> tuple[str | None, str | None]:
        """Method to get reference and mutated sequence."""
        return self.extract_sequence(), self.mutate_sequence()
