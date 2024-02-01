import logging
from abc import ABC, abstractmethod
from pathlib import Path
from typing import List, Optional, Any, Tuple, Union, TextIO, BinaryIO

import dnaio

from .files import InputFiles, OutputFiles, open_raise_limit, FileFormat
from .utils import Progress
from .modifiers import (
    SingleEndModifier,
    PairedEndModifier,
    PairedEndModifierWrapper,
    ModificationInfo,
)
from .predicates import DiscardUntrimmed, Predicate
from .steps import PairedEndSink, PairedEndFilter, SingleEndStep

logger = logging.getLogger()


class Pipeline(ABC):
    """
    Processing pipeline that loops over reads and applies modifiers and filters
    """

    n_adapters = 0
    paired = False

    def __init__(self) -> None:
        self._steps: List[Any] = []
        self._static_steps: List[Any] = []
        self._input_file_format: Optional[FileFormat] = None
        self._infiles: Optional[InputFiles] = None
        self._outfiles: Optional[OutputFiles] = None
        self._demultiplexer = None
        self._textiowrappers: List[TextIO] = []

        # Filter settings
        self.discard_trimmed = False
        self.discard_untrimmed = False

    def _open_writer(
        self,
        *files: Optional[BinaryIO],
        force_fasta: Optional[bool] = None,
    ):
        # The files must already be file-like objects because we don’t want to
        # take care of threads and compression levels here.
        for f in files:
            assert not isinstance(f, (str, bytes, Path))
        if len(files) == 2 and files[1] is None:
            files = files[:1]
            interleaved = True
        else:
            interleaved = False
        assert self._input_file_format is not None
        return open_raise_limit(
            dnaio.open,
            *files,
            mode="w",
            qualities=self._input_file_format.has_qualities(),
            fileformat="fasta" if force_fasta else None,
            interleaved=interleaved,
        )

    def flush(self) -> None:
        for f in self._textiowrappers:
            f.flush()

    def close(self) -> None:
        self._close_input()
        self._close_output()

    def _close_input(self) -> None:
        if self._infiles is not None:
            self._infiles.close()

    def _close_output(self) -> None:
        for f in self._textiowrappers:
            f.close()
        # Closing a TextIOWrapper also closes the underlying file, so
        # this closes some files a second time.
        if self._outfiles is not None:
            self._outfiles.close()

    @abstractmethod
    def process_reads(
        self,
        infiles: InputFiles,
        progress: Optional[Progress] = None,
    ) -> Tuple[int, int, Optional[int]]:
        pass


class SingleEndPipeline(Pipeline):
    """
    Processing pipeline for single-end reads
    """

    def __init__(
        self,
        modifiers: List[SingleEndModifier],
        steps: List[SingleEndStep],
    ):
        super().__init__()
        self._modifiers: List[SingleEndModifier] = modifiers
        self._steps = steps

    def process_reads(
        self,
        infiles: InputFiles,
        progress: Optional[Progress] = None,
    ) -> Tuple[int, int, Optional[int]]:
        """Run the pipeline. Return statistics"""
        self._infiles = infiles
        self._reader = infiles.open()
        for i, step in enumerate(self._steps, 1):
            logger.debug("Pipeline step %d: %s", i, step)

        n = 0  # no. of processed reads
        total_bp = 0
        for read in self._reader:
            n += 1
            if n % 10000 == 0 and progress is not None:
                progress.update(10000)
            total_bp += len(read)
            info = ModificationInfo(read)
            for modifier in self._modifiers:
                read = modifier(read, info)
            for filter_ in self._static_steps + self._steps:
                read = filter_(read, info)
                if read is None:
                    break
        if progress is not None:
            progress.update(n % 10000)
        return (n, total_bp, None)


class PairedEndPipeline(Pipeline):
    """
    Processing pipeline for paired-end reads.
    """

    paired = True

    def __init__(
        self,
        input_file_format: FileFormat,
        modifiers: List[
            Union[
                PairedEndModifier,
                Tuple[Optional[SingleEndModifier], Optional[SingleEndModifier]],
            ]
        ],
        pair_filter_mode: str,
        steps,
    ):
        super().__init__()
        self._input_file_format = input_file_format
        self._modifiers: List[PairedEndModifier] = []
        self._steps = steps
        self._pair_filter_mode = pair_filter_mode
        self._reader = None
        # Whether to ignore pair_filter mode for discard-untrimmed filter
        self.override_untrimmed_pair_filter = False
        self._add_modifiers(modifiers)

    def _add_modifiers(self, modifiers):
        for modifier in modifiers:
            if isinstance(modifier, tuple):
                self._add_two_single_modifiers(*modifier)
            else:
                self._add_modifier(modifier)

    def _add_two_single_modifiers(
        self,
        modifier1: Optional[SingleEndModifier],
        modifier2: Optional[SingleEndModifier],
    ) -> None:
        """
        Add two single-end modifiers that modify R1 and R2, respectively.
        One of them can be None, in which case the modifier
        is only applied to the respective other read.
        """
        if modifier1 is None and modifier2 is None:
            raise ValueError("Not both modifiers can be None")
        self._modifiers.append(PairedEndModifierWrapper(modifier1, modifier2))

    def _add_modifier(self, modifier: PairedEndModifier) -> None:
        """Add a Modifier (without wrapping it in a PairedEndModifierWrapper)"""
        self._modifiers.append(modifier)

    def process_reads(
        self,
        infiles: InputFiles,
        progress: Optional[Progress] = None,
    ) -> Tuple[int, int, Optional[int]]:
        self._infiles = infiles
        self._reader = infiles.open()
        n = 0  # no. of processed reads
        total1_bp = 0
        total2_bp = 0
        assert self._reader is not None
        for reads in self._reader:
            n += 1
            if n % 10000 == 0 and progress is not None:
                progress.update(10000)
            read1, read2 = reads
            total1_bp += len(read1)
            total2_bp += len(read2)
            info1 = ModificationInfo(read1)
            info2 = ModificationInfo(read2)
            for modifier in self._modifiers:
                reads = modifier(*reads, info1, info2)  # type: ignore
            for filter_ in self._static_steps + self._steps:
                reads = filter_(*reads, info1, info2)
                if reads is None:
                    break
        if progress is not None:
            progress.update(n % 10000)
        return (n, total1_bp, total2_bp)

    def _make_filter(
        self,
        predicate1: Optional[Predicate],
        predicate2: Optional[Predicate],
        writer,
        pair_filter_mode=None,
    ):
        if pair_filter_mode is None:
            pair_filter_mode = self._pair_filter_mode
        return PairedEndFilter(
            predicate1, predicate2, writer, pair_filter_mode=pair_filter_mode
        )

    def _make_untrimmed_filter(self, writer):
        """
        Return a different filter wrapper when adapters were given only for R1
        or only for R2 (then override_untrimmed_pair_filter will be set)
        """
        return self._make_filter(
            DiscardUntrimmed(),
            DiscardUntrimmed(),
            writer,
            pair_filter_mode="both" if self.override_untrimmed_pair_filter else None,
        )

    def _final_filter(self, outfiles):
        writer = self._open_writer(
            outfiles.out,
            outfiles.out2,
            force_fasta=outfiles.force_fasta,
        )
        return PairedEndSink(writer)
