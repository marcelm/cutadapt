from collections import defaultdict, Counter
from typing import DefaultDict, Tuple


class ReadLengthStatistics:
    """
    Keep track of the lengths of written reads or read pairs
    """

    def __init__(self) -> None:
        # It would be more natural to use a Counter, but a
        # defaultdict is much faster
        self._written_lengths1: DefaultDict[int, int] = defaultdict(int)
        self._written_lengths2: DefaultDict[int, int] = defaultdict(int)

    def update(self, read) -> None:
        """Add a single-end read to the statistics"""
        self._written_lengths1[len(read)] += 1

    def update2(self, read1, read2) -> None:
        """Add a paired-end read to the statistics"""
        self._written_lengths1[len(read1)] += 1
        self._written_lengths2[len(read2)] += 1

    def written_reads(self) -> int:
        """Return number of written reads or read pairs"""
        return sum(self._written_lengths1.values())

    def written_bp(self) -> Tuple[int, int]:
        return (
            self._compute_total_bp(self._written_lengths1),
            self._compute_total_bp(self._written_lengths2),
        )

    def written_lengths(self) -> Tuple[Counter, Counter]:
        return (Counter(self._written_lengths1), Counter(self._written_lengths2))

    @staticmethod
    def _compute_total_bp(counts: DefaultDict[int, int]) -> int:
        return sum(length * count for length, count in counts.items())

    def __iadd__(self, other):
        written_lengths1, written_lengths2 = other.written_lengths()
        for length, count in written_lengths1.items():
            self._written_lengths1[length] += count
        for length, count in written_lengths2.items():
            self._written_lengths2[length] += count
        return self
