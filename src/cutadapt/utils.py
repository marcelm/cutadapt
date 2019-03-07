import re
import sys
import time
import multiprocessing


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


class Progress:
    """
    Write progress
    """

    def __init__(self, every=1):
        """
        :param every: minimum time to wait in seconds between progress updates
        """
        if not sys.stderr.isatty():
            # Disable progress
            self.update = self.do_nothing
            self._show = False
        else:
            self._show = True

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

    def do_nothing(self, total, _final=False):
        pass

    def update(self, total, _final=False):
        """
        """
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
        if self._show:
            print(file=sys.stderr)
