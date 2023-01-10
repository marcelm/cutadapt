import sys

import dnaio

from ._qc import NUMBER_OF_NUCS, NUMBER_OF_PHREDS, QCMetrics

A, C, G, T, N = [ord(c) & 0b111 for c in "ACGTN"]
PHRED_TO_ERROR_RATE = [10 ** (-p /10) for p in range(NUMBER_OF_PHREDS)]


def analyse_metrics(metrics: QCMetrics):
    matrix = metrics.count_table_view().cast("Q")
    total_number_of_reads = metrics.number_of_reads
    table_size = NUMBER_OF_NUCS * NUMBER_OF_PHREDS
    sequence_length = len(matrix) // table_size
    length_table = [0 for _ in range(sequence_length)]
    for sequence_pos, table_offset in enumerate(range(0, len(matrix), table_size)):
        table = matrix[table_offset:table_offset+table_size]
        for phred, row_offset in enumerate(range(0, table_size, NUMBER_OF_NUCS)):
            nucs = table[row_offset:row_offset+NUMBER_OF_NUCS]
        length_table[sequence_pos] = sum(table)


if __name__ == "__main__":  # pragma: no cover
    metrics = QCMetrics()
    with dnaio.open(sys.argv[1]) as reader:
        for read in reader:
            metrics.add_read(read)
    analyse_metrics(metrics)

