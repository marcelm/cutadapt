import sys

import dnaio

from ._qc import NUMBER_OF_NUCS, NUMBER_OF_PHREDS, QCMetrics

A, C, G, T, N = [ord(c) & 0b111 for c in "ACGTN"]
PHRED_TO_ERROR_RATE = [10 ** (-p /10) for p in range(NUMBER_OF_PHREDS)]


def analyse_metrics(v: memoryview):
    matrix = v.cast("Q")
    table_size = NUMBER_OF_NUCS * NUMBER_OF_PHREDS
    sequence_length = len(matrix) // table_size
    for i in range(0, len(matrix), table_size):
        table = matrix[i:i+table_size]
        for phred, offset in enumerate(range(0, table_size, NUMBER_OF_NUCS)):
            nucs = table[offset:offset+NUMBER_OF_NUCS]


if __name__ == "__main__":  # pragma: no cover
    metrics = QCMetrics()
    with dnaio.open(sys.argv[1]) as reader:
        for read in reader:
            metrics.add_read(read)
    analyse_metrics(metrics.count_table_view())

