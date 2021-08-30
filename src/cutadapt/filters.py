"""
Filtering criteria
"""
from abc import ABC, abstractmethod

from .qualtrim import expected_errors
from .modifiers import ModificationInfo


# Constants used when returning from a Filter’s __call__ method to improve
# readability (it is unintuitive that "return True" means "discard the read").
DISCARD = True
KEEP = False


class Predicate(ABC):
    @abstractmethod
    def __call__(self, read, info: ModificationInfo) -> bool:
        """
        Return True if the filtering criterion matches.
        """


class TooShortReadFilter(Predicate):
    name: str = "too_short"

    def __init__(self, minimum_length):
        self.minimum_length = minimum_length

    def __repr__(self):
        return "TooShortReadFilter(minimum_length={})".format(self.minimum_length)

    def __call__(self, read, info: ModificationInfo):
        return len(read) < self.minimum_length


class TooLongReadFilter(Predicate):
    name: str = "too_long"

    def __init__(self, maximum_length):
        self.maximum_length = maximum_length

    def __repr__(self):
        return "TooLongReadFilter(maximum_length={})".format(self.maximum_length)

    def __call__(self, read, info: ModificationInfo):
        return len(read) > self.maximum_length


class MaximumExpectedErrorsFilter(Predicate):
    """
    Discard reads whose expected number of errors, according to the quality
    values, exceeds a threshold.

    The idea comes from usearch's -fastq_maxee parameter
    (http://drive5.com/usearch/).
    """
    name: str = "too_many_expected_errors"

    def __init__(self, max_errors):
        self.max_errors = max_errors

    def __repr__(self):
        return "MaximumExpectedErrorsFilter(max_errors={})".format(self.max_errors)

    def __call__(self, read, info: ModificationInfo):
        """Return True when the read should be discarded"""
        return expected_errors(read.qualities) > self.max_errors


class NContentFilter(Predicate):
    """
    Discard a read if it has too many 'N' bases. It handles both raw counts
    of Ns as well as proportions. Note, for raw counts, it is a 'greater than' comparison,
    so a cutoff of '1' will keep reads with a single N in it.
    """
    name: str = "too_many_n"

    def __init__(self, count):
        """
        Count -- if it is below 1.0, it will be considered a proportion, and above and equal to
        1 will be considered as discarding reads with a number of N's greater than this cutoff.
        """
        assert count >= 0
        self.is_proportion = count < 1.0
        self.cutoff = count

    def __repr__(self):
        return "NContentFilter(cutoff={}, is_proportion={})".format(
            self.cutoff, self.is_proportion)

    def __call__(self, read, info: ModificationInfo):
        """Return True when the read should be discarded"""
        n_count = read.sequence.lower().count('n')
        if self.is_proportion:
            if len(read) == 0:
                return False
            return n_count / len(read) > self.cutoff
        else:
            return n_count > self.cutoff


class CasavaFilter(Predicate):
    """
    Remove reads that fail the CASAVA filter. These have header lines that
    look like ``xxxx x:Y:x:x`` (with a ``Y``). Reads that pass the filter
    have an ``N`` instead of ``Y``.

    Reads with unrecognized headers are kept.
    """
    name: str = "casava_filtered"

    def __repr__(self):
        return "CasavaFilter()"

    def __call__(self, read, info: ModificationInfo):
        _, _, right = read.name.partition(' ')
        return right[1:4] == ':Y:'  # discard if :Y: found


class DiscardUntrimmedFilter(Predicate):
    """
    Return True if read is untrimmed.
    """
    name: str = "discard_untrimmed"

    def __repr__(self):
        return "DiscardUntrimmedFilter()"

    def __call__(self, read, info: ModificationInfo):
        return not info.matches


class DiscardTrimmedFilter(Predicate):
    """
    Return True if read is trimmed.
    """
    name: str = "discard_trimmed"

    def __repr__(self):
        return "DiscardTrimmedFilter()"

    def __call__(self, read, info: ModificationInfo):
        return bool(info.matches)
