import sys

import dnaio

from ._qc import QCMetrics

MAX_PHRED = 93

if __name__ == "__main__":  # pragma: no cover
    metrics = QCMetrics()
    with dnaio.open(sys.argv[1]) as reader:
        for read in reader:
            metrics.add_read(read)
    pass

