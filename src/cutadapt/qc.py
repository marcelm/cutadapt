import math
import sys

import dnaio

from ._qc import NUMBER_OF_NUCS, NUMBER_OF_PHREDS, QCMetrics

N, A, C, G, T = 0, 1, 2, 3, 4
PHRED_TO_ERROR_RATE = [10 ** (-p /10) for p in range(NUMBER_OF_PHREDS)]


def analyse_metrics(metrics: QCMetrics):
    matrix = metrics.count_table_view().cast("Q")
    table_size = NUMBER_OF_NUCS * NUMBER_OF_PHREDS
    sequence_length = len(matrix) // table_size

    length_counts = [0 for _ in range(sequence_length + 1)]
    length_counts[0] = metrics.number_of_reads
    qualities = [[0.0 for _ in range(sequence_length)] for _ in range(5)]
    base_content = [[0.0 for _ in range(sequence_length)] for _ in range(5)]
    gc_content = [0.0 for _ in range(sequence_length)]
    grand_total_bases = 0
    grand_total_gc = 0
    grand_total_at = 0

    for sequence_pos, table_offset in enumerate(range(0, len(matrix), table_size)):
        error_rates = [0.0 for _ in range(5)]
        base_counts = [0 for _ in range(5)]
        table = matrix[table_offset:table_offset+table_size]
        for phred, row_offset in enumerate(range(0, table_size, NUMBER_OF_NUCS)):
            nucs = table[row_offset:row_offset+NUMBER_OF_NUCS]
            for n_index, nuc_count in enumerate(nucs):
                error_rates[n_index] += PHRED_TO_ERROR_RATE[phred]
                base_counts[n_index] += nuc_count
        total_bases = sum(base_counts)
        for i in range(5):
            base_count = base_counts[i]
            average_error = error_rates[i] / base_count
            qualities[i][sequence_pos] = -10 * math.log10(average_error)
            base_content[i][sequence_pos] = base_count / total_bases
        at = base_counts[A] + base_counts[T]
        gc = base_counts[G] + base_counts[C]
        # Don't include N in gc content calculation
        gc_content[sequence_pos] = gc / (at + gc)
        grand_total_at += at
        grand_total_gc += gc
        length_counts[sequence_pos] += total_bases
        grand_total_bases += total_bases


if __name__ == "__main__":  # pragma: no cover
    metrics = QCMetrics()
    with dnaio.open(sys.argv[1]) as reader:
        for read in reader:
            metrics.add_read(read)
    analyse_metrics(metrics)

