import re
import sys
import time
import multiprocessing
import logging
import locale


logger = logging.getLogger(__name__)

try:
    "µ".encode(locale.getpreferredencoding())
    MICRO = "µ"
except UnicodeEncodeError:
    MICRO = "u"


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
            "{total:13,d} reads @ {per_item:5.1F} {micro}s/read; {per_minute:6.2F} M reads/minute"
            "".format(
                hours=hours,
                minutes=minutes,
                seconds=seconds,
                total=self._n,
                per_item=per_item * 1e6,
                micro=MICRO,
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
