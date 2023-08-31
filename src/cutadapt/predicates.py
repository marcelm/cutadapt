"""
Filtering criteria (predicates)
"""
from abc import ABC, abstractmethod

from .qualtrim import expected_errors
from .modifiers import ModificationInfo


class Predicate(ABC):
    @abstractmethod
    def test(self, read, info: ModificationInfo) -> bool:
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
        return "".join(
            ("_" + ch.lower() if ch.isupper() else ch) for ch in cls.__name__
        )[1:]


class TooShort(Predicate):
    """Select reads that are shorter than the specified minimum length"""

    def __init__(self, minimum_length: int):
        self.minimum_length = minimum_length

    def __repr__(self):
        return f"TooShort(minimum_length={self.minimum_length})"

    def test(self, read, info: ModificationInfo):
        return len(read) < self.minimum_length


class TooLong(Predicate):
    """Select reads that are longer than the specified maximum length"""

    def __init__(self, maximum_length: int):
        self.maximum_length = maximum_length

    def __repr__(self):
        return f"TooLong(maximum_length={self.maximum_length})"

    def test(self, read, info: ModificationInfo):
        return len(read) > self.maximum_length


class TooManyExpectedErrors(Predicate):
    """
    Select reads whose expected number of errors, according to the quality
    values, exceeds a threshold.

    The idea comes from usearch's -fastq_maxee parameter
    (http://drive5.com/usearch/).
    """

    def __init__(self, max_errors: float):
        self.max_errors = max_errors

    def __repr__(self):
        return f"TooManyExpectedErrors(max_errors={self.max_errors})"

    def test(self, read, info: ModificationInfo):
        return expected_errors(read.qualities) > self.max_errors


class TooHighAverageErrorRate(Predicate):
    """
    Select reads that have an average error rate above the threshold.
    This works better than TooManyExpectedErrors for reads that are expected to
    have varying lengths, such as for long read sequencing technologies.
    """

    def __init__(self, max_error_rate: float):
        if not 0.0 < max_error_rate < 1.0:
            raise ValueError(
                f"max_error_rate must be between 0.0 and 1.0, got {max_error_rate}."
            )
        self.max_error_rate = max_error_rate

    def __repr__(self):
        return f"TooHighAverageErrorRate(max_error_rate={self.max_error_rate}"

    def test(self, read, info: ModificationInfo):
        return (expected_errors(read.qualities) / len(read)) > self.max_error_rate


class TooManyN(Predicate):
    """
    Select reads that have too many 'N' bases.

    Both a raw count or a proportion (relative to the sequence length) can be used.
    """

    def __init__(self, count: float):
        """
        count -- if it is below 1.0, it will be considered a proportion, and above and equal to
        1 will be considered as discarding reads with a number of N's greater than this cutoff.
        """
        assert count >= 0
        self.is_proportion = count < 1.0
        self.cutoff = count

    def __repr__(self):
        return f"TooManyN(cutoff={self.cutoff}, is_proportion={self.is_proportion})"

    def test(self, read, info: ModificationInfo):
        n_count = read.sequence.lower().count("n")
        if self.is_proportion:
            if len(read) == 0:
                return False
            return n_count / len(read) > self.cutoff
        else:
            return n_count > self.cutoff


class CasavaFiltered(Predicate):
    """
    Select reads that have failed the CASAVA filter according to the read header.
    The headers look like ``xxxx x:Y:x:x`` (with a ``Y``). Reads that pass the filter
    have an ``N`` instead of ``Y``.

    Reads with unrecognized headers are not selected.
    """

    def __repr__(self):
        return "CasavaFiltered()"

    def test(self, read, info: ModificationInfo):
        _, _, right = read.name.partition(" ")
        return right[1:4] == ":Y:"  # discard if :Y: found


class DiscardUntrimmed(Predicate):
    """
    Select reads for which no adapter match was found
    """

    def __repr__(self):
        return "DiscardUntrimmed()"

    def test(self, read, info: ModificationInfo):
        return not info.matches


class DiscardTrimmed(Predicate):
    """
    Select reads for which at least one adapter match was found
    """

    def __repr__(self):
        return "DiscardTrimmed()"

    def test(self, read, info: ModificationInfo):
        return bool(info.matches)
