"""
Read processing steps (filters)

Classes for writing and filtering of processed reads.

A Filter is a callable that has the read as its only argument. If it is called,
it returns True if the read should be filtered (discarded), and False if not.

To be used, a filter needs to be wrapped in one of the redirector classes.
They are called so because they can redirect filtered reads to a file if so
desired. They also keep statistics.

To determine what happens to a read, a list of redirectors with different
filters is created and each redirector is called in turn until one returns True.
The read is then assumed to have been "consumed", that is, either written
somewhere or filtered (should be discarded).
"""
from abc import ABC, abstractmethod
from collections import defaultdict, Counter
from typing import DefaultDict, Tuple, Dict, Optional, Any

from .modifiers import ModificationInfo
from .utils import reverse_complemented_sequence

# Constants used when returning from a Filter’s __call__ method to improve
# readability (it is unintuitive that "return True" means "discard the read").
DISCARD = True
KEEP = False


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


class SingleEndStep(ABC):
    @abstractmethod
    def __call__(self, read, info: ModificationInfo) -> bool:
        """
        Process a single read
        """


class PairedEndStep(ABC):
    @abstractmethod
    def __call__(self, read1, read2, info1: ModificationInfo, info2: ModificationInfo) -> bool:
        """
        Process read pair (read1, read2)
        """


class SingleEndFinalStep(SingleEndStep, ABC):
    def __init__(self):
        super().__init__()
        self.statistics = ReadLengthStatistics()

    def update_statistics(self, read) -> None:
        self.statistics.update(read)


class PairedEndFinalStep(PairedEndStep, ABC):
    def __init__(self):
        super().__init__()
        self.statistics = ReadLengthStatistics()

    def update_statistics(self, read1, read2):
        self.statistics.update2(read1, read2)


class NoFilter(SingleEndFinalStep):
    """
    No filtering, just send each read to the given writer.
    """
    def __init__(self, writer):
        super().__init__()
        self.writer = writer

    def __repr__(self):
        return "NoFilter({})".format(self.writer)

    def __call__(self, read, info: ModificationInfo):
        self.writer.write(read)
        self.update_statistics(read)
        return DISCARD


class PairedNoFilter(PairedEndFinalStep):
    """
    No filtering, just send each paired-end read to the given writer.
    """
    def __init__(self, writer):
        super().__init__()
        self.writer = writer

    def __repr__(self):
        return "PairedNoFilter({})".format(self.writer)

    def __call__(self, read1, read2, info1: ModificationInfo, info2: ModificationInfo):
        self.writer.write(read1, read2)
        self.update_statistics(read1, read2)
        return DISCARD


class Redirector(SingleEndStep):
    """
    Redirect discarded reads to the given writer. This is for single-end reads.
    """
    def __init__(self, writer, filter: SingleEndStep, filter2=None):
        super().__init__()
        # TODO filter2 should really not be here
        self.filtered = 0
        self.writer = writer
        self.filter = filter

    def __repr__(self):
        return "Redirector(writer={}, filter={})".format(self.writer, self.filter)

    def __call__(self, read, info: ModificationInfo):
        if self.filter(read, info):
            self.filtered += 1
            if self.writer is not None:
                self.writer.write(read)
            return DISCARD
        return KEEP


class PairedRedirector(PairedEndStep):
    """
    Redirect paired-end reads matching a filtering criterion to a writer.
    Different filtering styles are supported, differing by which of the
    two reads in a pair have to fulfill the filtering criterion.
    """
    def __init__(self, writer, filter, filter2, pair_filter_mode='any'):
        """
        pair_filter_mode -- these values are allowed:
            'any': The pair is discarded if any read matches.
            'both': The pair is discarded if both reads match.
            'first': The pair is discarded if the first read matches.
        """
        super().__init__()
        if pair_filter_mode not in ('any', 'both', 'first'):
            raise ValueError("pair_filter_mode must be 'any', 'both' or 'first'")
        self._pair_filter_mode = pair_filter_mode
        self.filtered = 0
        self.writer = writer
        self.filter = filter
        self.filter2 = filter2
        if filter2 is None:
            self._is_filtered = self._is_filtered_first
        elif filter is None:
            self._is_filtered = self._is_filtered_second
        elif pair_filter_mode == 'any':
            self._is_filtered = self._is_filtered_any
        elif pair_filter_mode == 'both':
            self._is_filtered = self._is_filtered_both
        else:
            self._is_filtered = self._is_filtered_first

    def __repr__(self):
        return "PairedRedirector(writer={}, filter={}, filter2={}, pair_filter_mode='{}')".format(
            self.writer, self.filter, self.filter2, self._pair_filter_mode)

    def _is_filtered_any(self, read1, read2, info1: ModificationInfo, info2: ModificationInfo):
        return self.filter(read1, info1) or self.filter2(read2, info2)

    def _is_filtered_both(self, read1, read2, info1: ModificationInfo, info2: ModificationInfo):
        return self.filter(read1, info1) and self.filter2(read2, info2)

    def _is_filtered_first(self, read1, read2, info1: ModificationInfo, info2: ModificationInfo):
        return self.filter(read1, info1)

    def _is_filtered_second(self, read1, read2, info1: ModificationInfo, info2: ModificationInfo):
        return self.filter2(read2, info2)

    def __call__(self, read1, read2, info1: ModificationInfo, info2: ModificationInfo):
        if self._is_filtered(read1, read2, info1, info2):
            self.filtered += 1
            if self.writer is not None:
                self.writer.write(read1, read2)
            return DISCARD
        return KEEP


class RestFileWriter(SingleEndStep):
    def __init__(self, file):
        self.file = file

    def __call__(self, read, info):
        # TODO this fails with linked adapters
        if info.matches:
            rest = info.matches[-1].rest()
            if len(rest) > 0:
                print(rest, read.name, file=self.file)
        return KEEP


class WildcardFileWriter(SingleEndStep):
    def __init__(self, file):
        self.file = file

    def __call__(self, read, info):
        # TODO this fails with linked adapters
        if info.matches:
            print(info.matches[-1].wildcards(), read.name, file=self.file)
        return KEEP


class InfoFileWriter(SingleEndStep):
    RC_MAP = {None: "", True: "1", False: "0"}

    def __init__(self, file):
        self.file = file

    def __call__(self, read, info: ModificationInfo):
        current_read = info.original_read
        if info.is_rc:
            current_read = reverse_complemented_sequence(current_read)
        if info.matches:
            for match in info.matches:
                for info_record in match.get_info_records(current_read):
                    # info_record[0] is the read name suffix
                    print(read.name + info_record[0], *info_record[1:], self.RC_MAP[info.is_rc],
                        sep='\t', file=self.file)
                current_read = match.trimmed(current_read)
        else:
            seq = read.sequence
            qualities = read.qualities if read.qualities is not None else ''
            print(read.name, -1, seq, qualities, sep='\t', file=self.file)

        return KEEP


class Demultiplexer(SingleEndFinalStep):
    """
    Demultiplex trimmed reads. Reads are written to different output files
    depending on which adapter matches.

    Untrimmed reads are sent to writers[None] if that key exists.
    """
    def __init__(self, writers: Dict[Optional[str], Any]):
        """
        writers maps an adapter name to a writer
        """
        super().__init__()
        self._writers = writers
        self._untrimmed_writer = self._writers.get(None, None)

    def __call__(self, read, info):
        """
        Write the read to the proper output file according to the most recent match
        """
        if info.matches:
            name = info.matches[-1].adapter.name
            self.update_statistics(read)
            self._writers[name].write(read)
        elif self._untrimmed_writer is not None:
            self.update_statistics(read)
            self._untrimmed_writer.write(read)
        return DISCARD


class PairedDemultiplexer(PairedEndFinalStep):
    """
    Demultiplex trimmed paired-end reads. Reads are written to different output files
    depending on which adapter (in read 1) matches.
    """
    def __init__(self, writers: Dict[Optional[str], Any]):
        super().__init__()
        self._writers = writers
        self._untrimmed_writer = self._writers.get(None, None)

    def __call__(self, read1, read2, info1: ModificationInfo, info2: ModificationInfo):
        assert read2 is not None
        if info1.matches:
            name = info1.matches[-1].adapter.name  # type: ignore
            self.update_statistics(read1, read2)
            self._writers[name].write(read1, read2)
        elif self._untrimmed_writer is not None:
            self.update_statistics(read1, read2)
            self._untrimmed_writer.write(read1, read2)
        return DISCARD


class CombinatorialDemultiplexer(PairedEndFinalStep):
    """
    Demultiplex paired-end reads depending on which adapter matches, taking into account
    matches on R1 and R2.
    """
    def __init__(self, writers: Dict[Tuple[Optional[str], Optional[str]], Any]):
        """
        Adapter names of the matches on R1 and R2 will be used to look up the writer in the
        writers dict. If there is no match on a read, None is used in the lookup instead
        of the name. Missing dictionary keys are ignored and can be used to discard
        read pairs.
        """
        super().__init__()
        self._writers = writers

    def __call__(self, read1, read2, info1, info2):
        """
        Write the read to the proper output file according to the most recent matches both on
        R1 and R2
        """
        assert read2 is not None
        name1 = info1.matches[-1].adapter.name if info1.matches else None
        name2 = info2.matches[-1].adapter.name if info2.matches else None
        key = (name1, name2)
        if key in self._writers:
            self.update_statistics(read1, read2)
            self._writers[key].write(read1, read2)
        return DISCARD
