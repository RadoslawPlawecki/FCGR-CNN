"""
@author: Radosław Pławecki
"""

class NsHandler:
    """
    Class to analyse and trim leading/trailing N's in a sequence.
    """
    def __init__(self, sequence: str):
        """
        Initialize the class instance.
        """
        self.sequence = sequence
        self.n = len(sequence)
        self.start = 0
        self.end = self.n

    def count_leading_ns(self) -> int:
        """
        Count leading N's and get a new starting point.
        """
        start = 0
        while start < self.n and self.sequence[start] == 'N':
            start += 1
        self.start = start
        return start

    def count_trailing_ns(self) -> int:
        """
        Count trailing N's and get a new ending point.
        """
        end = self.n - 1
        while end >= self.start and self.sequence[end] == 'N':
            end -= 1
        self.end = end + 1
        return self.n - 1 - end

    def count_internal_ns(self) -> int:
        """
        Count internal N's.
        """
        return self.sequence[self.start:self.end].count('N')

    def get_core_sequence(self) -> str | None:
        """
        Trim leading and trailing N's.
        """
        if self.start <= self.end:
            return self.sequence[self.start:self.end]
        return None
