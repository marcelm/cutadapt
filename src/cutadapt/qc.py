import math
import sys
from typing import List

import dnaio

import pygal  # type: ignore

from ._qc import NUMBER_OF_NUCS, NUMBER_OF_PHREDS, QCMetrics

N, A, C, G, T = 0, 1, 2, 3, 4
PHRED_TO_ERROR_RATE = [10 ** (-p / 10) for p in range(NUMBER_OF_PHREDS)]


class QCMetricsReport:
    total_reads: int
    total_bases: int
    q20_bases: int
    q30_bases: int
    _total_gc: int
    _total_at: int
    max_length: int
    sequence_lengths: List[int]
    per_base_qualities: List[List[float]]
    mean_qualities: List[float]
    base_content: List[List[float]]
    gc_content: List[float]

    @property
    def total_gc_fraction(self) -> float:
        return self._total_gc / (self._total_at + self._total_gc)

    def mean_length(self) -> float:
        total_lengths = 0
        for length, number_of_reads in enumerate(self.sequence_lengths):
            total_lengths += length * number_of_reads
        return total_lengths / self.total_reads

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
            x_labels=list(range(0, self.max_length + 1)),
            x_labels_major=list(range(0, self.max_length, 10)),
            show_minor_x_labels=False,
            truncate_label=-1,
            width=1000,
            explicit_size=True,
            disable_xml_declaration=True,
        )
        plot.add("Length", self.sequence_lengths)
        return plot.render(is_unicode=True)

    def base_content_plot(self) -> str:
        plot = pygal.Line(
            title="Base content",
            dots_size=1,
            x_labels=list(range(1, self.max_length + 1)),
            x_labels_major=list(range(0, self.max_length, 10)),
            show_minor_x_labels=False,
            truncate_label=-1,
            width=1000,
            explicit_size=True,
            disable_xml_declaration=True,
        )
        plot.add("GC", self.gc_content)
        plot.add("A", self.base_content[A])
        plot.add("G", self.base_content[G])
        plot.add("C", self.base_content[C])
        plot.add("T", self.base_content[T])
        plot.add("N", self.base_content[N])
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
            <td align="right">{self.q20_bases} ({self.q20_bases * 100 / self.total_bases:.2f}%)</td>
        </tr>
        <tr>
            <td>Q30 bases</td>
            <td align="right">{self.q30_bases} ({self.q30_bases * 100 / self.total_bases:.2f}%)</td>
        </tr>
        <tr><td>GC content</td><td align="right">{self.total_gc_fraction * 100:.2f}%</td></tr>
        </table>
        <h2>Quality scores</h2>
        {self.per_base_quality_plot()}
        </html>
        <h2>Sequence length distribution</h2>
        {self.sequence_length_distribution_plot()}
        <h2>Base content</h2>
        {self.base_content_plot()}
        """

    def __init__(self, metrics: QCMetrics):
        """Aggregate all data from a QCMetrics counter"""
        matrix = metrics.count_table_view().cast("Q")
        table_size = NUMBER_OF_NUCS * NUMBER_OF_PHREDS
        sequence_length = metrics.max_length
        length_counts = [0 for _ in range(sequence_length + 1)]
        length_counts[0] = metrics.number_of_reads
        qualities = [[0.0 for _ in range(sequence_length)] for _ in range(5)]
        mean_qualities = [0.0 for _ in range(sequence_length)]
        base_content = [[0.0 for _ in range(sequence_length)] for _ in range(5)]
        gc_content = [0.0 for _ in range(sequence_length)]
        grand_total_bases = 0
        grand_total_gc = 0
        grand_total_at = 0
        grand_total_q20 = 0
        grand_total_q30 = 0

        for sequence_pos, table_offset in enumerate(range(0, len(matrix), table_size)):
            error_rates = [0.0 for _ in range(5)]
            base_counts = [0 for _ in range(5)]
            table = matrix[table_offset : table_offset + table_size]
            for phred, row_offset in enumerate(range(0, table_size, NUMBER_OF_NUCS)):
                nucs = table[row_offset : row_offset + NUMBER_OF_NUCS]
                for n_index, nuc_count in enumerate(nucs):
                    if phred >= 20:
                        grand_total_q20 += nuc_count
                    if phred >= 30:
                        grand_total_q30 += nuc_count
                    error_rates[n_index] += PHRED_TO_ERROR_RATE[phred] * nuc_count
                    base_counts[n_index] += nuc_count
            total_bases = sum(base_counts)
            for i in range(NUMBER_OF_NUCS):
                base_count = base_counts[i]
                base_content[i][sequence_pos] = base_count / total_bases
                if base_count:
                    average_error = error_rates[i] / base_count
                    qualities[i][sequence_pos] = -10 * math.log10(average_error)
            mean_qualities[sequence_pos] = -10 * math.log10(
                sum(error_rates) / total_bases
            )
            at = base_counts[A] + base_counts[T]
            gc = base_counts[G] + base_counts[C]
            # Don't include N in gc content calculation
            gc_content[sequence_pos] = gc / (at + gc)
            grand_total_at += at
            grand_total_gc += gc
            # Bases at pos 0 are for a sequence of length 1
            length_counts[sequence_pos + 1] += total_bases
            grand_total_bases += total_bases
        sequence_lengths = [0 for _ in range(sequence_length + 1)]
        previous_count = 0
        for i in range(sequence_length, 0, -1):
            number_of_sequences_at_least_length = length_counts[i]
            number_of_sequences_exactly_length = (
                number_of_sequences_at_least_length - previous_count
            )
            sequence_lengths[i] = number_of_sequences_exactly_length
            previous_count = number_of_sequences_at_least_length

        self.total_reads = metrics.number_of_reads
        self.total_bases = grand_total_bases
        self.q20_bases = grand_total_q20
        self.q30_bases = grand_total_q30
        self._total_gc = grand_total_gc
        self._total_at = grand_total_at
        self.max_length = metrics.max_length
        self.sequence_lengths = sequence_lengths
        self.per_base_qualities = qualities
        self.mean_qualities = mean_qualities
        self.base_content = base_content
        self.gc_content = gc_content


if __name__ == "__main__":  # pragma: no cover
    metrics = QCMetrics()
    with dnaio.open(sys.argv[1]) as reader:
        for read in reader:
            metrics.add_read(read)
    report = QCMetricsReport(metrics)
    print(report.html_report())
