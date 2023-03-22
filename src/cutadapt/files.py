import errno
import io
import sys
from typing import BinaryIO, Optional, Dict, Tuple

import dnaio
from xopen import xopen

from cutadapt.utils import logger

try:
    import resource
except ImportError:
    # Windows
    resource = None  # type: ignore


def xopen_rb_raise_limit(path: str):
    """
    Open a (possibly compressed) file for reading in binary mode, trying to avoid the
    "Too many open files" problem using `open_raise_limit`.
    """
    mode = "rb"
    f = open_raise_limit(xopen, path, mode, threads=0)
    logger.debug("Opening '%s', mode '%s' with xopen resulted in %s", path, mode, f)
    return f


def open_raise_limit(func, *args, **kwargs):
    """
    Run 'func' (which should be some kind of open() function) and return its result.
    If "Too many open files" occurs, increase limit and try again.
    """
    try:
        f = func(*args, **kwargs)
    except OSError as e:
        if e.errno == errno.EMFILE:  # Too many open files
            logger.debug("Too many open files, attempting to raise soft limit")
            raise_open_files_limit(8)
            f = func(*args, **kwargs)
        else:
            raise
    return f


def raise_open_files_limit(n):
    if resource is None:
        return
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    soft = min(soft + n, hard)
    resource.setrlimit(resource.RLIMIT_NOFILE, (soft, hard))


class FileOpener:
    def __init__(self, compression_level: int = 6, threads: Optional[int] = None):
        """
        threads -- no. of external compression threads.
            0: write in-process
            None: min(cpu_count(), 4)
        """
        self.compression_level = compression_level
        self.threads = threads

    def xopen(self, path, mode):
        threads = self.threads if "w" in mode else 0
        f = open_raise_limit(
            xopen, path, mode, compresslevel=self.compression_level, threads=threads
        )
        logger.debug("Opening '%s', mode '%s' with xopen resulted in %s", path, mode, f)
        return f

    def xopen_or_none(self, path, mode):
        """Return opened file or None if the path is None"""
        if path is None:
            return None
        return self.xopen(path, mode)

    def xopen_pair(self, path1: str, path2: Optional[str], mode):
        if path1 is None and path2 is not None:
            raise ValueError(
                "When giving paths for paired-end files, only providing the second"
                " file is not supported"
            )
        file1 = self.xopen_or_none(path1, mode)
        file2 = self.xopen_or_none(path2, mode)
        return file1, file2

    def dnaio_open(self, *args, **kwargs):
        kwargs["opener"] = self.xopen
        f = dnaio.open(*args, **kwargs)
        if not isinstance(args[0], io.BytesIO):
            logger.debug(
                "Opening %r, mode '%s' with dnaio resulted in %s",
                args[0],
                kwargs["mode"],
                f,
            )
        return f

    def dnaio_open_raise_limit(self, *args, **kwargs):
        """
        Open a FASTA/FASTQ file for writing. If it fails because the number of open files
        would be exceeded, try to raise the soft limit and re-try.
        """
        return open_raise_limit(self.dnaio_open, *args, **kwargs)


class InputFiles:
    def __init__(
        self,
        *files: BinaryIO,
        interleaved: bool = False,
    ):
        self._files = files
        self.interleaved = interleaved
        for f in self._files:
            assert f is not None

    def open(self):
        return dnaio.open(*self._files, interleaved=self.interleaved, mode="r")

    def close(self) -> None:
        for file in self._files:
            file.close()


class InputPaths:
    def __init__(self, *paths: str, interleaved: bool = False):
        self.paths = paths
        self.interleaved = interleaved

    def open(self) -> InputFiles:
        files = [xopen_rb_raise_limit(path) for path in self.paths]
        return InputFiles(*files, interleaved=self.interleaved)


class OutputFiles:
    """
    The attributes are either None or open file-like objects except for demultiplex_out
    and demultiplex_out2, which are dictionaries that map an adapter name
    to a file-like object.
    """

    def __init__(
        self,
        out: Optional[BinaryIO] = None,
        out2: Optional[BinaryIO] = None,
        untrimmed: Optional[BinaryIO] = None,
        untrimmed2: Optional[BinaryIO] = None,
        too_short: Optional[BinaryIO] = None,
        too_short2: Optional[BinaryIO] = None,
        too_long: Optional[BinaryIO] = None,
        too_long2: Optional[BinaryIO] = None,
        info: Optional[BinaryIO] = None,
        rest: Optional[BinaryIO] = None,
        wildcard: Optional[BinaryIO] = None,
        demultiplex_out: Optional[Dict[str, BinaryIO]] = None,
        demultiplex_out2: Optional[Dict[str, BinaryIO]] = None,
        combinatorial_out: Optional[Dict[Tuple[str, str], BinaryIO]] = None,
        combinatorial_out2: Optional[Dict[Tuple[str, str], BinaryIO]] = None,
        force_fasta: Optional[bool] = None,
    ):
        self.out = out
        self.out2 = out2
        self.untrimmed = untrimmed
        self.untrimmed2 = untrimmed2
        self.too_short = too_short
        self.too_short2 = too_short2
        self.too_long = too_long
        self.too_long2 = too_long2
        self.info = info
        self.rest = rest
        self.wildcard = wildcard
        self.demultiplex_out = demultiplex_out
        self.demultiplex_out2 = demultiplex_out2
        self.combinatorial_out = combinatorial_out
        self.combinatorial_out2 = combinatorial_out2
        self.force_fasta = force_fasta

    def __iter__(self):
        for f in [
            self.out,
            self.out2,
            self.untrimmed,
            self.untrimmed2,
            self.too_short,
            self.too_short2,
            self.too_long,
            self.too_long2,
            self.info,
            self.rest,
            self.wildcard,
        ]:
            if f is not None:
                yield f
        for outs in (
            self.demultiplex_out,
            self.demultiplex_out2,
            self.combinatorial_out,
            self.combinatorial_out2,
        ):
            if outs is not None:
                for f in outs.values():
                    assert f is not None
                    yield f

    def as_bytesio(self) -> "OutputFiles":
        """
        Create a new OutputFiles instance that has BytesIO instances for each non-None output file
        """
        result = OutputFiles(force_fasta=self.force_fasta)
        for attr in (
            "out",
            "out2",
            "untrimmed",
            "untrimmed2",
            "too_short",
            "too_short2",
            "too_long",
            "too_long2",
            "info",
            "rest",
            "wildcard",
        ):
            if getattr(self, attr) is not None:
                setattr(result, attr, io.BytesIO())
        for attr in (
            "demultiplex_out",
            "demultiplex_out2",
            "combinatorial_out",
            "combinatorial_out2",
        ):
            if getattr(self, attr) is not None:
                setattr(result, attr, dict())
                for k, v in getattr(self, attr).items():
                    getattr(result, attr)[k] = io.BytesIO()
        return result

    def close(self) -> None:
        """Close all output files that are not stdout"""
        for f in self:
            if f is sys.stdout or f is sys.stdout.buffer:
                continue
            f.close()
