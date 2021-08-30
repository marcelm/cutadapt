"""
Filtering criteria
"""
from abc import ABC, abstractmethod

from .qualtrim import expected_errors
from .modifiers import ModificationInfo


class Predicate(ABC):
    @abstractmethod
    def __call__(self, read, info: ModificationInfo) -> bool:
        """
        Return True if the filtering criterion matches.
        """

    @classmethod
    def descriptive_identifier(cls) -> str:
        """
        Return a short name for this predicate based on the class name such as "too_long",
        "too_many_expected_errors".
        This is used as identifier in the JSON report.
        """
        return "".join(("_" + ch.lower() if ch.isupper() else ch) for ch in cls.__name__)[1:]


class TooShort(Predicate):
    def __init__(self, minimum_length: int):
        self.minimum_length = minimum_length

    def __repr__(self):
        return f"TooShort(minimum_length={self.minimum_length})"

    def __call__(self, read, info: ModificationInfo):
        return len(read) < self.minimum_length


class TooLong(Predicate):
    def __init__(self, maximum_length: int):
        self.maximum_length = maximum_length

    def __repr__(self):
        return f"TooLong(maximum_length={self.maximum_length})"

    def __call__(self, read, info: ModificationInfo):
        return len(read) > self.maximum_length


class TooManyExpectedErrors(Predicate):
    """
    Discard reads whose expected number of errors, according to the quality
    values, exceeds a threshold.

    The idea comes from usearch's -fastq_maxee parameter
    (http://drive5.com/usearch/).
    """
    def __init__(self, max_errors: float):
        self.max_errors = max_errors

    def __repr__(self):
        return f"TooManyExpectedErrors(max_errors={self.max_errors})"

    def __call__(self, read, info: ModificationInfo):
        """Return True when the read should be discarded"""
        return expected_errors(read.qualities) > self.max_errors


class TooManyN(Predicate):
    """
    Discard a read if it has too many 'N' bases. It handles both raw counts
    of Ns as well as proportions. Note, for raw counts, it is a 'greater than' comparison,
    so a cutoff of '1' will keep reads with a single N in it.
    """
    def __init__(self, count: float):
        """
        Count -- if it is below 1.0, it will be considered a proportion, and above and equal to
        1 will be considered as discarding reads with a number of N's greater than this cutoff.
        """
        assert count >= 0
        self.is_proportion = count < 1.0
        self.cutoff = count

    def __repr__(self):
        return f"TooManyN(cutoff={self.cutoff}, is_proportion={self.is_proportion})"

    def __call__(self, read, info: ModificationInfo):
        """Return True when the read should be discarded"""
        n_count = read.sequence.lower().count('n')
        if self.is_proportion:
            if len(read) == 0:
                return False
            return n_count / len(read) > self.cutoff
        else:
            return n_count > self.cutoff


class CasavaFiltered(Predicate):
    """
    Remove reads that fail the CASAVA filter. These have header lines that
    look like ``xxxx x:Y:x:x`` (with a ``Y``). Reads that pass the filter
    have an ``N`` instead of ``Y``.

    Reads with unrecognized headers are kept.
    """
    def __repr__(self):
        return "CasavaFiltered()"

    def __call__(self, read, info: ModificationInfo):
        _, _, right = read.name.partition(' ')
        return right[1:4] == ':Y:'  # discard if :Y: found


class DiscardUntrimmed(Predicate):
    """
    Return True if read is untrimmed.
    """
    def __repr__(self):
        return "DiscardUntrimmed()"

    def __call__(self, read, info: ModificationInfo):
        return not info.matches


class DiscardTrimmed(Predicate):
    """
    Return True if read is trimmed.
    """
    def __repr__(self):
        return "DiscardTrimmed()"

    def __call__(self, read, info: ModificationInfo):
        return bool(info.matches)
