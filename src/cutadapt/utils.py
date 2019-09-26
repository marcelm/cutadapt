import re
import sys
import time
import multiprocessing

import dnaio


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
    except IOError:
        pass

    return multiprocessing.cpu_count()


def raise_open_files_limit(n):
    try:
        import resource
    except ImportError:
        # Windows
        return
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    soft += n
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
        self._n = 0
        # Time at which the progress was last updated
        self._time = time.time()
        self._start_time = self._time
        self._animation = self.scissors()

    @staticmethod
    def scissors(width=10):
        while True:
            for is_reverse, rang in [(False, range(width + 1)), (True, range(width + 1))]:
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

    def update(self, total, _final=False):
        current_time = time.time()
        if _final:
            time_delta = current_time - self._start_time
            delta = total
        else:
            time_delta = current_time - self._time
            delta = total - self._n
        if delta < 1:
            return
        if not _final:
            if time_delta < self._every:
                return
            if total <= self._n:
                return

        t = current_time - self._start_time
        hours = int(t) // 3600
        minutes = (int(t) - hours * 3600) // 60
        seconds = int(t) % 60
        per_second = delta / time_delta
        per_item = time_delta / delta

        print(
            "\r"
            "{animation} {hours:02d}:{minutes:02d}:{seconds:02d} "
            "{total:13,d} reads  @  {per_item:7.1F} µs/read; {per_minute:6.2F} M reads/minute"
            "".format(
                hours=hours, minutes=minutes, seconds=seconds,
                total=total, per_item=per_item * 1E6, per_minute=per_second * 60 / 1E6,
                animation=next(self._animation)),
            end="", file=sys.stderr
        )
        self._n = total
        self._time = current_time

    def stop(self, total):
        """
        Print final progress reflecting the final total
        """
        self.update(total, _final=True)
        print(file=sys.stderr)


class DummyProgress:
    """
    Has the same interface as Progress, but does not print anything
    """
    def update(self, total, _final=False):
        pass

    def stop(self, total):
        pass


_REVCOMP_TRANS = str.maketrans(
    "ACGTUMRWSYKVHDBNacgtumrwsykvhdbn",
    "TGCAAKYWSRMBDHVNtgcaakywsrmbdhvn",
)


def reverse_complement(s: str):
    return s.translate(_REVCOMP_TRANS)[::-1]


def reverse_complemented_sequence(sequence: dnaio.Sequence):
    if sequence.qualities is None:
        qualities = None
    else:
        qualities = sequence.qualities[::-1]
    return dnaio.Sequence(sequence.name, reverse_complement(sequence.sequence), qualities)
