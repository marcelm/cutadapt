import io
import re
import sys
import time
import errno
import multiprocessing
import logging
from typing import Optional

from xopen import xopen
import dnaio

try:
    import resource
except ImportError:
    # Windows
    resource = None  # type: ignore


logger = logging.getLogger(__name__)


def available_cpu_count():
    """
    Return the number of available virtual or physical CPUs on this system.
    The number of available CPUs can be smaller than the total number of CPUs
    when the cpuset(7) mechanism is in use, as is the case on some cluster
    systems.

    Adapted from http://stackoverflow.com/a/1006301/715090
    """
    try:
        with open("/proc/self/status") as f:
            status = f.read()
        m = re.search(r"(?m)^Cpus_allowed:\s*(.*)$", status)
        if m:
            res = bin(int(m.group(1).replace(",", ""), 16)).count("1")
            if res > 0:
                return min(res, multiprocessing.cpu_count())
    except OSError:
        pass

    return multiprocessing.cpu_count()


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


class Progress:
    """
    Print an animated progress report to sys.stderr
    """

    def __init__(self, every=1):
        """
        every: minimum time to wait in seconds between progress updates
        """
        self._every = every
        self._animation = self.scissors()
        self._n = 0
        self._start_time = time.time()
        # Time at which the progress was last updated
        self._last_time = self._start_time
        self._last_n = 0

    def __repr__(self):
        return (
            f"Progress(_n={self._n}, elapsed={self._last_time - self._start_time:.3f})"
        )

    @staticmethod
    def scissors(width=10):
        while True:
            for is_reverse, rang in [
                (False, range(width + 1)),
                (True, range(width + 1)),
            ]:
                for position in rang:
                    for is_open in (True, False):
                        left = " " * position
                        right = "-" * (width - position)
                        if is_reverse:
                            sc = ">8" if is_open else "=8"
                            left, right = right, left
                        else:
                            sc = "8<" if is_open else "8="
                        yield "[" + left + sc + right + "]"

    def update(self, increment, _final=False):
        self._n += increment
        current_time = time.time()
        if _final:
            time_delta = current_time - self._start_time
            delta = self._n
        else:
            time_delta = current_time - self._last_time
            delta = self._n - self._last_n
        if delta < 1:
            return
        if time_delta == 0:
            return
        if not _final:
            if time_delta < self._every:
                return

        t = current_time - self._start_time
        hours = int(t) // 3600
        minutes = (int(t) - hours * 3600) // 60
        seconds = int(t) % 60
        per_second = delta / time_delta
        per_item = time_delta / delta

        animation = next(self._animation)
        if _final:
            animation = "Done".ljust(len(animation))
        print(
            "\r"
            "{animation} {hours:02d}:{minutes:02d}:{seconds:02d} "
            "{total:13,d} reads @ {per_item:5.1F} µs/read; {per_minute:6.2F} M reads/minute"
            "".format(
                hours=hours,
                minutes=minutes,
                seconds=seconds,
                total=self._n,
                per_item=per_item * 1e6,
                per_minute=per_second * 60 / 1e6,
                animation=animation,
            ),
            end="",
            file=sys.stderr,
        )
        self._last_time = current_time
        self._last_n = self._n

    def close(self):
        """
        Print final progress reflecting the final total
        """
        self.update(0, _final=True)
        print(file=sys.stderr)


class DummyProgress(Progress):
    """
    Does not print anything
    """

    def update(self, increment, _final=False):
        pass

    def close(self):
        pass


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
