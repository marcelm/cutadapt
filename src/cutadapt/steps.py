"""
Steps of the read output pipeline

After all read modifications have been done, a read is written to at
most one output file. For this, a pipeline represented as a list of "steps"
(SingleEndSteps or PairedEndSteps) is used. Each pipeline step can consume
(discard) a read or pass it on to the next step.

Steps are added to the pipeline in a certain order:

1. First RestFileWriter, InfoFileWriter, WildcardFileWriter because
   they should see all reads before filtering.
2. Filters come next. These are implemented as SingleEndFilter or PairedEndFilter
   instances with an appropriate Predicate. Filters can optionally send each
   consumed/filtered read to an output file.
3. The last pipeline step should be one of the "Sinks", which consume all reads.
   Demultiplexers are sinks, for example.
"""

from abc import ABC, abstractmethod
from typing import Tuple, Dict, Optional, Any

from dnaio import SequenceRecord

from .predicates import Predicate
from .modifiers import ModificationInfo
from .statistics import ReadLengthStatistics

RecordPair = Tuple[SequenceRecord, SequenceRecord]


class SingleEndStep(ABC):
    @abstractmethod
    def __call__(self, read, info: ModificationInfo) -> Optional[SequenceRecord]:
        """
        Process a single read. Return the processed read or None to indicate that
        the read has been consumed and should thus not be passed on to subsequent
        steps.
        """


class PairedEndStep(ABC):
    @abstractmethod
    def __call__(
        self, read1, read2, info1: ModificationInfo, info2: ModificationInfo
    ) -> Optional[RecordPair]:
        """
        Process (read1, read2). Return the processed read pair or None if
        the read pair has been "consumed" (filtered or written to an output file)
        and should thus not be passed on to subsequent steps.
        """


class HasStatistics(ABC):
    """
    Used for the final steps (sinks), which need to keep track of read length statistics
    """

    @abstractmethod
    def get_statistics(self) -> ReadLengthStatistics:
        pass


class HasFilterStatistics(ABC):
    @abstractmethod
    def filtered(self) -> int:
        """Return number of filtered reads or read pairs"""

    @abstractmethod
    def descriptive_identifier(self) -> str:
        """Name used in statistics"""


class SingleEndFilter(SingleEndStep, HasFilterStatistics):
    """
    A pipeline step that can filter reads, can redirect filtered ones to a writer, and
    counts how many were filtered.
    """

    def __init__(self, predicate: Predicate, writer):
        super().__init__()
        self._filtered = 0
        self._predicate = predicate
        self._writer = writer

    def __repr__(self):
        return f"SingleEndFilter(predicate={self._predicate}, writer={self._writer})"

    def descriptive_identifier(self) -> str:
        return self._predicate.descriptive_identifier()

    def filtered(self) -> int:
        return self._filtered

    def __call__(self, read, info: ModificationInfo) -> Optional[SequenceRecord]:
        if self._predicate.test(read, info):
            self._filtered += 1
            if self._writer is not None:
                self._writer.write(read)
            return None
        return read


class PairedEndFilter(PairedEndStep, HasFilterStatistics):
    """
    A pipeline step that can filter paired-end reads, redirect them to a file, and counts
    how many read pairs were filtered.

    Different filtering styles are supported, differing by which of the
    two reads in a pair have to fulfill the filtering criterion.
    """

    def __init__(
        self,
        predicate1: Optional[Predicate],
        predicate2: Optional[Predicate],
        writer,
        pair_filter_mode="any",
    ):
        """
        pair_filter_mode -- these values are allowed:
            'any': The pair is discarded if any read matches.
            'both': The pair is discarded if both reads match.
            'first': The pair is discarded if the first read matches.
        """
        super().__init__()
        if pair_filter_mode not in ("any", "both", "first"):
            raise ValueError("pair_filter_mode must be 'any', 'both' or 'first'")
        self._pair_filter_mode = pair_filter_mode
        self._filtered = 0
        self.predicate1 = predicate1
        self.predicate2 = predicate2
        self.writer = writer
        self._is_filtered: Any
        if predicate2 is None:
            self._is_filtered = self._is_filtered_first
        elif predicate1 is None:
            self._is_filtered = self._is_filtered_second
        elif pair_filter_mode == "any":
            self._is_filtered = self._is_filtered_any
        elif pair_filter_mode == "both":
            self._is_filtered = self._is_filtered_both
        else:
            self._is_filtered = self._is_filtered_first

    def __repr__(self):
        return (
            f"PairedEndFilter(predicate1={self.predicate1}, predicate2={self.predicate2}, "
            f"writer={self.writer}, pair_filter_mode='{self._pair_filter_mode}')"
        )

    def descriptive_identifier(self) -> str:
        if self.predicate1 is not None:
            return self.predicate1.descriptive_identifier()
        else:
            assert self.predicate2 is not None
            return self.predicate2.descriptive_identifier()

    def filtered(self) -> int:
        return self._filtered

    def _is_filtered_any(
        self, read1, read2, info1: ModificationInfo, info2: ModificationInfo
    ) -> bool:
        return self.predicate1.test(read1, info1) or self.predicate2.test(read2, info2)  # type: ignore

    def _is_filtered_both(
        self, read1, read2, info1: ModificationInfo, info2: ModificationInfo
    ) -> bool:
        return self.predicate1.test(read1, info1) and self.predicate2.test(read2, info2)  # type: ignore

    def _is_filtered_first(
        self, read1, read2, info1: ModificationInfo, info2: ModificationInfo
    ) -> bool:
        return self.predicate1.test(read1, info1)  # type: ignore

    def _is_filtered_second(
        self, read1, read2, info1: ModificationInfo, info2: ModificationInfo
    ) -> bool:
        return self.predicate2.test(read2, info2)  # type: ignore

    def __call__(
        self, read1, read2, info1: ModificationInfo, info2: ModificationInfo
    ) -> Optional[RecordPair]:
        if self._is_filtered(read1, read2, info1, info2):
            self._filtered += 1
            if self.writer is not None:
                self.writer.write(read1, read2)
            return None
        return (read1, read2)


class RestFileWriter(SingleEndStep):
    def __init__(self, file):
        self._file = file

    def __repr__(self):
        return f"RestFileWriter(file={self._file})"

    def __call__(self, read, info) -> Optional[SequenceRecord]:
        # TODO this fails with linked adapters
        if info.matches:
            rest = info.matches[-1].rest()
            if len(rest) > 0:
                print(rest, read.name, file=self._file)
        return read


class WildcardFileWriter(SingleEndStep):
    def __init__(self, file):
        self._file = file

    def __repr__(self):
        return f"WildcardFileWriter(file={self._file})"

    def __call__(self, read, info) -> Optional[SequenceRecord]:
        # TODO this fails with linked adapters
        if info.matches:
            print(info.matches[-1].wildcards(), read.name, file=self._file)
        return read


class InfoFileWriter(SingleEndStep):
    RC_MAP = {None: "", True: "1", False: "0"}

    def __init__(self, file):
        self._file = file

    def __repr__(self):
        return f"InfoFileWriter(file={self._file})"

    def __call__(self, read, info: ModificationInfo) -> Optional[SequenceRecord]:
        current_read = info.original_read
        if info.is_rc:
            current_read = current_read.reverse_complement()
        if info.matches:
            for match in info.matches:
                for info_record in match.get_info_records(current_read):
                    # info_record[0] is the read name suffix
                    print(
                        read.name + info_record[0],
                        *info_record[1:],
                        self.RC_MAP[info.is_rc],
                        sep="\t",
                        file=self._file,
                    )
                current_read = match.trimmed(current_read)
        else:
            seq = read.sequence
            qualities = read.qualities if read.qualities is not None else ""
            print(read.name, -1, seq, qualities, sep="\t", file=self._file)

        return read


class PairedSingleEndStep(PairedEndStep):
    """
    Wrap a SingleEndStep as a PairedEndStep

    The wrapped step is called with the first read
    """

    def __init__(self, step: SingleEndStep):
        self._step = step

    def __repr__(self):
        return f"PairedSingleEndStep(step={self._step})"

    def __call__(self, read1, read2, info1, info2) -> Optional[RecordPair]:
        _ = read2  # intentionally ignored
        _ = info2
        result = self._step(read1, info1)
        if result is None:
            return None
        return (result, read2)


# The following steps are used as final step in a pipeline.
# They send each read or read pair to the final intended output file,
# and they all track the lengths of written reads.


class SingleEndSink(SingleEndStep, HasStatistics):
    """
    Send each read to a writer and keep read length statistics.
    This is used as the last step in a pipeline.
    """

    def __init__(self, writer):
        super().__init__()
        self.writer = writer
        self._statistics = ReadLengthStatistics()

    def __repr__(self):
        return f"NoFilter({self.writer})"

    def __call__(self, read, info: ModificationInfo) -> Optional[SequenceRecord]:
        self.writer.write(read)
        self._statistics.update(read)
        return None

    def get_statistics(self) -> ReadLengthStatistics:
        return self._statistics


class PairedEndSink(PairedEndStep, HasStatistics):
    """
    Send each read pair to a writer and keep read length statistics.
    This is used as the last step in a pipeline.
    """

    def __init__(self, writer):
        super().__init__()
        self.writer = writer
        self._statistics = ReadLengthStatistics()

    def __repr__(self):
        return f"PairedNoFilter({self.writer})"

    def __call__(
        self, read1, read2, info1: ModificationInfo, info2: ModificationInfo
    ) -> Optional[RecordPair]:
        self.writer.write(read1, read2)
        self._statistics.update2(read1, read2)
        return None

    def get_statistics(self) -> ReadLengthStatistics:
        return self._statistics


class Demultiplexer(SingleEndStep, HasStatistics, HasFilterStatistics):
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
        self._statistics = ReadLengthStatistics()
        self._filtered = 0

    def __repr__(self):
        return f"<Demultiplexer len(writers)={len(self._writers)}>"

    def __call__(self, read, info) -> Optional[SequenceRecord]:
        """
        Write the read to the proper output file according to the most recent match
        """
        if info.matches:
            name = info.matches[-1].adapter.name
            self._statistics.update(read)
            self._writers[name].write(read)
        elif self._untrimmed_writer is not None:
            self._statistics.update(read)
            self._untrimmed_writer.write(read)
        else:
            self._filtered += 1
        return None

    def descriptive_identifier(self) -> str:
        return "discard_untrimmed"

    def get_statistics(self) -> ReadLengthStatistics:
        return self._statistics

    def filtered(self) -> int:
        return self._filtered


class PairedDemultiplexer(PairedEndStep, HasStatistics, HasFilterStatistics):
    """
    Demultiplex trimmed paired-end reads. Reads are written to different output files
    depending on which adapter (in read 1) matches.
    """

    def __init__(self, writers: Dict[Optional[str], Any]):
        super().__init__()
        self._writers = writers
        self._untrimmed_writer = self._writers.get(None, None)
        self._statistics = ReadLengthStatistics()
        self._filtered = 0

    def __call__(
        self, read1, read2, info1: ModificationInfo, info2: ModificationInfo
    ) -> Optional[RecordPair]:
        assert read2 is not None
        if info1.matches:
            name = info1.matches[-1].adapter.name  # type: ignore
            self._statistics.update2(read1, read2)
            self._writers[name].write(read1, read2)
        elif self._untrimmed_writer is not None:
            self._statistics.update2(read1, read2)
            self._untrimmed_writer.write(read1, read2)
        else:
            self._filtered += 1
        return None

    def descriptive_identifier(self) -> str:
        return "discard_untrimmed"

    def get_statistics(self) -> ReadLengthStatistics:
        return self._statistics

    def filtered(self) -> int:
        return self._filtered


class CombinatorialDemultiplexer(PairedEndStep, HasStatistics):
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
        self._statistics = ReadLengthStatistics()

    def __call__(self, read1, read2, info1, info2) -> Optional[RecordPair]:
        """
        Write the read to the proper output file according to the most recent matches both on
        R1 and R2
        """
        assert read2 is not None
        name1 = info1.matches[-1].adapter.name if info1.matches else None
        name2 = info2.matches[-1].adapter.name if info2.matches else None
        key = (name1, name2)
        if key in self._writers:
            self._statistics.update2(read1, read2)
            self._writers[key].write(read1, read2)
        return None

    def get_statistics(self) -> ReadLengthStatistics:
        return self._statistics
