import errno
import io
import sys
import os
from abc import ABC, abstractmethod
from enum import Enum
from typing import BinaryIO, Optional, Dict, List, TextIO, Any
import dnaio
from xopen import xopen
import smart_open
import logging
import json
from cutadapt.utils import logger

try:
    import resource
except ImportError:
    # Windows
    resource = None  # type: ignore


def xopen_rb_raise_limit(path: str, transport_params: str = ""):
    """
    Open a (possibly compressed) file for reading in binary mode, trying to avoid the
    "Too many open files" problem using `open_raise_limit`.
    """
    mode = "rb"
    # transfer options : string of key=value,key=value or a file path: convert to dictionary
    transport_params = get_transport_params(path, transport_params)
    # smart_open to automatically open remote files, disable auto-compression
    f = open_raise_limit(
        smart_open.open,
        path,
        mode,
        compression="disable",
        transport_params=transport_params,
    )
    logging.getLogger("smart_open").setLevel(logging.WARNING)
    # pass through to xopen to handle compression
    f = open_raise_limit(xopen, f, mode, threads=4)
    logger.debug("Opening '%s', mode '%s' with xopen resulted in %s", path, mode, f)
    return f


def get_transport_params(path, transport_params):
    if not transport_params:
        return {}
    # load from json file
    if os.path.isfile(transport_params):
        with open(transport_params) as f:
            transport_params = json.load(f)
    else:
        transport_params = json.loads(transport_params)

    return transport_params


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
    def __init__(
        self,
        compression_level: int = 1,
        threads: Optional[int] = None,
        transport_params: str = "",
    ):
        """
        threads -- no. of external compression threads.
            0: write in-process
            None: min(cpu_count(), 4)
        """
        self.compression_level = compression_level
        self.threads = threads
        self.transport_params = transport_params

    def smart_open(self, path, mode):
        # get transport params for smart_open
        transport_params = get_transport_params(path, self.transport_params)
        # smart_open to automatically open remote files, disable auto-compression
        f = open_raise_limit(
            smart_open.open,
            path,
            mode,
            compression="disable",
            transport_params=transport_params,
        )
        logging.getLogger("smart_open").setLevel(logging.ERROR)
        logger.debug(
            "Opening output '%s', mode '%s' with smart_open resulted in %s",
            path,
            mode,
            f,
        )
        return f

    def xopen(self, path, mode):
        threads = self.threads if "w" in mode else 0
        # smart open to handle remote files
        f = self.smart_open(path, mode)
        # xopen to handle compression
        f = open_raise_limit(
            xopen, f, mode, compresslevel=self.compression_level, threads=threads
        )
        if "w" in mode:
            extra = f" (compression level {self.compression_level}, {threads} threads)"
        else:
            extra = ""
        logger.debug(
            "Opening output '%s', mode '%s'%s with xopen resulted in %s",
            path,
            mode,
            extra,
            f,
        )
        return f

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
    def __init__(
        self, *paths: str, interleaved: bool = False, transport_params: str = ""
    ):
        self.paths = paths
        self.interleaved = interleaved
        self.transport_params = transport_params

    def open(self) -> InputFiles:
        files = [
            xopen_rb_raise_limit(path, self.transport_params) for path in self.paths
        ]
        return InputFiles(*files, interleaved=self.interleaved)


class ProxyWriter(ABC):
    @abstractmethod
    def drain(self) -> List[bytes]:
        pass


class ProxyTextFile(ProxyWriter):
    """
    A file object for writing in text mode that is backed by a BytesIO object
    """

    def __init__(self):
        self._buffer = io.BytesIO()
        self._file = io.TextIOWrapper(self._buffer)

    def write(self, text):
        self._file.write(text)

    def drain(self) -> List[bytes]:
        self._file.flush()
        chunk = self._buffer.getvalue()
        self._buffer.seek(0)
        self._buffer.truncate()
        return [chunk]

    def __getstate__(self):
        """TextIOWrapper cannot be pickled. Just don’t include our state."""
        return True  # ensure __setstate__ is called

    def __setstate__(self, state):
        self.__init__()


class ProxyRecordWriter(ProxyWriter):
    """
    A writer for FASTA, FASTQ records etc. that is backed by a BytesIO object
    """

    def __init__(self, n_files: int, **kwargs):
        self._n_files = n_files
        self._kwargs = kwargs
        self._buffers = [io.BytesIO() for _ in range(n_files)]
        self._writer = open_raise_limit(dnaio.open, *self._buffers, mode="w", **kwargs)

    def write(self, *args, **kwargs):
        self._writer.write(*args, **kwargs)

    def drain(self) -> List[bytes]:
        chunks = [buf.getvalue() for buf in self._buffers]
        for buf in self._buffers:
            buf.seek(0)
            buf.truncate()
        return chunks

    def __getstate__(self):
        """Exclude the dnaio Reader class from the state"""
        return (self._n_files, self._kwargs)

    def __setstate__(self, state):
        n_files, kwargs = state
        self.__init__(n_files, **kwargs)


class OutputFiles:
    def __init__(
        self,
        *,
        proxied: bool,
        qualities: bool,
        interleaved: bool,
        file_opener: Optional[FileOpener] = None,
    ):
        self._file_opener: FileOpener = (
            file_opener if file_opener is not None else FileOpener()
        )
        self._binary_files: List[BinaryIO] = []
        self._binary_files_to_close: List[BinaryIO] = []
        self._text_files: List[TextIO] = []
        self._writers: List[Any] = []
        self._proxy_files: List[ProxyWriter] = []
        self._proxied = proxied
        self._to_close: List[BinaryIO] = []
        self._qualities = qualities
        self._interleaved = interleaved

    def open_text(self, path):
        # TODO
        # - serial runner needs only text_file
        # - parallel runner needs binary_file and proxy_file
        # split into SerialOutputFiles and ParallelOutputFiles?
        if self._proxied:
            binary_file = self._file_opener.xopen(path, "wb")
            self._binary_files.append(binary_file)
            self._binary_files_to_close.append(binary_file)
            proxy_file = ProxyTextFile()
            self._proxy_files.append(proxy_file)
            return proxy_file
        else:
            text_file = self._file_opener.xopen(path, "wt")
            self._text_files.append(text_file)
            return text_file

    def open_record_writer(
        self, *paths, interleaved: bool = False, force_fasta: bool = False
    ):
        kwargs: Dict[str, Any] = dict(
            qualities=self._qualities, interleaved=interleaved
        )
        if len(paths) not in (1, 2):
            raise ValueError("Expected one or two paths")
        if interleaved and len(paths) != 1:
            raise ValueError("Cannot write to two files when interleaved is True")
        if len(paths) == 1 and paths[0] == "-" and force_fasta:
            kwargs["fileformat"] = "fasta"
        if paths == (None,):
            paths = ("-",)
        for path in paths:
            assert path is not None
        binary_files = []
        for path in paths:
            binary_file = self._file_opener.xopen(path, "wb")
            binary_files.append(binary_file)
            self._binary_files.append(binary_file)
            self._binary_files_to_close.append(binary_file)
        if self._proxied:
            proxy_writer = ProxyRecordWriter(len(paths), **kwargs)
            self._proxy_files.append(proxy_writer)
            return proxy_writer
        else:
            writer = self._file_opener.dnaio_open(*binary_files, mode="w", **kwargs)
            self._writers.append(writer)
            return writer

    def open_stdout_record_writer(
        self, interleaved: bool = False, force_fasta: bool = False
    ):
        self._binary_files.append(sys.stdout.buffer)
        kwargs: Dict[str, Any] = dict(
            qualities=self._qualities, interleaved=interleaved
        )
        if force_fasta:
            kwargs["fileformat"] = "fasta"
        if self._proxied:
            proxy_writer = ProxyRecordWriter(1, **kwargs)
            self._proxy_files.append(proxy_writer)
            return proxy_writer
        else:
            writer = self._file_opener.dnaio_open(sys.stdout.buffer, mode="w", **kwargs)
            self._writers.append(writer)
            return writer

    def binary_files(self) -> List[BinaryIO]:
        return self._binary_files[:]

    def proxy_files(self) -> List[ProxyWriter]:
        return self._proxy_files

    def close(self) -> None:
        """Close all output files that are not stdout"""
        if not self._proxied:
            for f in self._text_files:
                f.close()
            for f in self._writers:
                f.close()
        for bf in self._binary_files_to_close:
            bf.close()


class FileFormat(Enum):
    FASTA = 1
    FASTQ = 2
    BAM = 3

    def has_qualities(self) -> bool:
        return self is FileFormat.FASTQ or self is FileFormat.BAM  # TODO BAM?


# TODO copied and adjusted from dnaio; upstream this
def detect_file_format(file: BinaryIO) -> FileFormat:
    if file.seekable():
        original_position = file.tell()
        magic = file.read(4)
        file.seek(original_position)
    else:
        # We cannot always use peek() because BytesIO objects do not suppert it
        magic = file.peek(4)[0:4]  # type: ignore
    if magic.startswith(b"@") or magic == b"":
        # Pretend FASTQ for empty input
        return FileFormat.FASTQ
    elif magic.startswith(b">") or magic.startswith(b"#"):
        # Some FASTA variants allow comments
        return FileFormat.FASTA
    elif magic == b"BAM\1":
        return FileFormat.BAM
    raise dnaio.exceptions.UnknownFileFormat(
        f"Input file format not recognized. The file starts with {magic!r}, "
        "but files in supported formats start with '>' (FASTA), '@' (FASTQ) or 'BAM'"
    )
