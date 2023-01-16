import array
import math
import sys
from typing import Iterator, Iterable, List, Tuple

import dnaio

import pygal  # type: ignore

from ._qc import NUMBER_OF_NUCS, NUMBER_OF_PHREDS, TABLE_SIZE, QCMetrics

N, A, C, G, T = 0, 1, 2, 3, 4
PHRED_TO_ERROR_RATE = [
    sum(10 ** (-p / 10) for p in range(start * 4, start * 4 + 4)) / 4
    for start in range(NUMBER_OF_PHREDS)]


def equidistant_ranges(length: int, parts: int) -> Iterator[Tuple[int, int]]:
    size = length // parts
    remainder = length % parts
    small_parts = parts - remainder
    start = 0
    for i in range(parts):
        part_size = size if i < small_parts else size + 1
        if part_size == 0:
            continue
        stop = start + part_size
        yield start, stop
        start = stop


class QCMetricsReport:
    raw_count_matrix: array.ArrayType
    aggregated_count_matrix = array.ArrayType
    raw_sequence_lengths = array.ArrayType
    _data_ranges: List[range]
    data_categories: List[str]
    max_length: int
    total_reads: int
    total_bases: int

    def __init__(self, metrics: QCMetrics, graph_resolution: int = 100):
        """Aggregate all data from a QCMetrics counter"""

        self.total_reads = metrics.number_of_reads
        self.max_length = metrics.max_length
        self.raw_count_matrix = array.array("Q")
        # Python will treat the memoryview as an iterable in the array constructor
        # use from_bytes instead for direct memcpy.
        self.raw_count_matrix.frombytes(metrics.count_table_view())
        # use bytes constructor to initialize the aggregated count matrix to 0.
        self.aggregated_count_matrix = array.array(
            "Q", bytes(8 * TABLE_SIZE * graph_resolution))

        self._data_ranges = [
            range(*r) for r in equidistant_ranges(metrics.max_length, graph_resolution)
        ]
        # Use one-based indexing for the graph categories. I.e. 1 is the first base.
        self.data_categories = [
            f"{r.start + 1}-{r.stop}"
            if r.start + 1 != r.stop else f"{r.start + 1}"
            for r in self._data_ranges
        ]

        matrix = memoryview(self.raw_count_matrix)
        categories_view = memoryview(self.aggregated_count_matrix)

        table_size = TABLE_SIZE
        for cat_index, category_range in enumerate(self._data_ranges):
            cat_offset = cat_index * table_size
            cat_view = categories_view[cat_offset: cat_offset + table_size]
            for table_index in category_range:
                offset = table_index * table_size
                table = matrix[offset: offset+table_size]
                for i, count in enumerate(table):
                    cat_view[i] += count
        self.total_bases = sum(self.aggregated_count_matrix)

        raw_sequence_lengths = array.array(
            "Q", bytes(8 * (self.max_length + 1)))
        raw_base_counts = array.array(
            "Q", bytes(8 * (self.max_length + 1)))
        # All reads have at least 0 bases
        raw_base_counts[0] = self.total_reads
        for i in range(self.max_length):
            table = matrix[i * 60:(i + 1)*60]
            raw_base_counts[i + 1] = sum(table)

        previous_count = 0
        for i in range(self.max_length, 0, -1):
            number_at_least = raw_base_counts[i]
            raw_sequence_lengths[i] = number_at_least - previous_count
            previous_count = number_at_least
        self.raw_sequence_lengths = raw_sequence_lengths

    def _tables(self) -> Iterator[memoryview]:
        category_view = memoryview(self.aggregated_count_matrix)
        for i in range(0, len(category_view), TABLE_SIZE):
            yield category_view[i:i + TABLE_SIZE]

    def base_content(self) -> List[List[float]]:
        content = [[0.0 for _ in range(len(self.data_categories))]
                        for _ in range(NUMBER_OF_NUCS)]
        for cat_index, table in enumerate(self._tables()):
            total = sum(table)
            for i in range(NUMBER_OF_NUCS):
                content[i][cat_index] = sum(table[i::NUMBER_OF_NUCS]) / total
        return content

    def total_gc_fraction(self) -> float:
        total_nucs = [
            sum(self.aggregated_count_matrix[i:len(self.aggregated_count_matrix):NUMBER_OF_NUCS])
            for i in range(NUMBER_OF_NUCS)
        ]
        at = total_nucs[A] + total_nucs[T]
        gc = total_nucs[G] + total_nucs[C]
        return gc / (at + gc)

    def q20_bases(self):
        q20s = 0
        for table in self._tables():
            q20s += sum(table[NUMBER_OF_NUCS * 5:])
        return q20s

    def q28_bases(self):
        q28s = 0
        for table in self._tables():
            q28s += sum(table[NUMBER_OF_NUCS * 7:])
        return q28s

    def mean_length(self) -> float:
        total_lengths = 0
        for length, number_of_reads in enumerate(self.raw_sequence_lengths):
            total_lengths += length * number_of_reads
        return total_lengths / self.total_reads

    def sequence_lengths(self):
        seqlength_view = memoryview(self.raw_sequence_lengths)[1:]
        lengths = [sum(seqlength_view[r.start:r.stop]) for r in self._data_ranges]
        return [self.raw_sequence_lengths[0]] + lengths

    def per_base_quality_plot(self) -> str:
        plot = pygal.Line(
            title="Per base sequence quality",
            dots_size=1,
            x_labels=list(range(1, self.max_length + 1)),
            x_labels_major=list(range(0, self.max_length, 10)),
            show_minor_x_labels=False,
            truncate_label=-1,
            width=1000,
            explicit_size=True,
            disable_xml_declaration=True,
        )
        plot.add("mean", self.mean_qualities)
        plot.add("A", self.per_base_qualities[A])
        plot.add("G", self.per_base_qualities[G])
        plot.add("C", self.per_base_qualities[C])
        plot.add("T", self.per_base_qualities[T])
        return plot.render(is_unicode=True)

    def sequence_length_distribution_plot(self) -> str:
        plot = pygal.Bar(
            title="Sequence length distribution",
            x_labels=["0"] + self.data_categories,
            truncate_label=-1,
            width=1000,
            explicit_size=True,
            disable_xml_declaration=True,
        )
        plot.add("Length", self.sequence_lengths())
        return plot.render(is_unicode=True)

    def base_content_plot(self) -> str:
        plot = pygal.Line(
            title="Base content",
            dots_size=1,
            x_labels=self.data_categories,
            truncate_label=-1,
            width=1000,
            explicit_size=True,
            disable_xml_declaration=True,
        )
        base_content = self.base_content()
        plot.add("A", base_content[A])
        plot.add("G", base_content[G])
        plot.add("C", base_content[C])
        plot.add("T", base_content[T])
        plot.add("N", base_content[N])
        return plot.render(is_unicode=True)

    def html_report(self):
        return f"""
        <html>
        <head>
            <meta http-equiv="content-type" content="text/html:charset=utf-8">
            <title>Cutadapt report</title>
        </head>
        <h1>Cutadapt report</h1>
        <h2>Summary</h2>
        <table>
        <tr><td>Mean length</td><td align="right">{self.mean_length()}</td></tr>
        <tr><td>total reads</td><td align="right">{self.total_reads}</td></tr>
        <tr><td>total bases</td><td align="right">{self.total_bases}</td></tr>
        <tr>
            <td>Q20 bases</td>
            <td align="right">{self.q20_bases()} ({self.q20_bases() * 100 / self.total_bases:.2f}%)</td>
        </tr>
        <tr>
            <td>Q28 bases</td>
            <td align="right">{self.q28_bases()} ({self.q28_bases() * 100 / self.total_bases:.2f}%)</td>
        </tr>
        <tr><td>GC content</td><td align="right">{self.total_gc_fraction() * 100:.2f}%</td></tr>
        </table>
        <h2>Quality scores</h2>
        {self.per_base_quality_plot()}
        </html>
        <h2>Sequence length distribution</h2>
        {self.sequence_length_distribution_plot()}
        <h2>Base content</h2>
        {self.base_content_plot()}
        </html>
        """


if __name__ == "__main__":  # pragma: no cover
    metrics = QCMetrics()
    with dnaio.open(sys.argv[1]) as reader:
        for read in reader:
            metrics.add_read(read)
    report = QCMetricsReport(metrics)
    print(report.html_report())
