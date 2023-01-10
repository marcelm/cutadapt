import array
import sys

import dnaio

from ._qc import QCMetrics

MAX_PHRED = 93

def analyse_metrics(v: memoryview):
    matrix = v.cast("Q")
    table_width = MAX_PHRED * 8
    sequence_length = len(matrix) // table_width
    for i in range(0, len(matrix), table_width):
        table = matrix[i:i+table_width]
        for qual, column_start in enumerate(range(0, len(table), 8)):
            column = table[column_start:column_start+8]
            print(qual, column.tolist())
        break


if __name__ == "__main__":  # pragma: no cover
    metrics = QCMetrics()
    with dnaio.open(sys.argv[1]) as reader:
        for read in reader:
            metrics.add_read(read)
    analyse_metrics(metrics.count_table_view())

