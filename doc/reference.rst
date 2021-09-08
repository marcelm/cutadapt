===============
Reference guide
===============

.. _json-report-format:

JSON report format
==================


* ``tag`` is always ``Cutadapt report``.
* ``schema_version`` is a tuple of two ints, the major and minor version.
  If additions are made to the schema, the minor version is increased. If backwards incompatible
  changes to the schema are made, the major version is increased.
* Keys only relevant for paired-end data or when using certain command-line options are
  always included, but when unused, get a value of ``null``.
* For adapters that allow partial matches, ``error_lengths`` describes the lengths up to which
  0, 1, 2 etc. errors are allowed. ``[9, 16]``: 0 errors up to a match of length 9, 1 error up to
  a match of length 16. The last number in this list is the length of the adapter sequence. For
  anchored adapter types, this is ``null``.
* ``adjacent_bases``: Statistics about adjacent bases are currently only kept for 3' adapters.
* ``dominant_adjacent_base`` is set to the appropriate nucleotide if the
  :ref:`report warns <warnbase>`: "The adapter is preceded by "x" extremely often."
  This is ``null`` if no such warning was printed.
* For paired-end data, numbers in the ``read_counts`` section are the number of *read pairs*.
* The adapter type is one of these strings:
     - ``"regular_five_prime"``
     - ``"regular_three_prime"``
     - ``"noninternal_five_prime"``
     - ``"noninternal_three_prime"``
     - ``"anchored_five_prime"``
     - ``"anchored_three_prime"``
     - ``"anywhere"``
     - ``"linked"``
* ``basepair_counts``

Example (slightly reformatted) ::

    {
      "tag": "Cutadapt report",
      "schema_version": [0, 1],
      "cutadapt_version": "3.5.dev121+g18b10ed.d20210908",
      "python_version": "3.8.10",
      "command_line_arguments": [
        "--json=out.cutadapt.json",
        "-m",
        "20",
        "--too-short-output=tooshort.fasta",
        "-a",
        "ACGTAC",
        "-q",
        "20",
        "--discard-trimmed",
        "-o",
        "/dev/null",
        "reads.fastq"
      ],
      "cores": 1,
      "input": {
        "path1": "reads.fastq",
        "path2": null,
        "paired": false,
        "interleaved": null
      },
      "read_counts": {
        "input": 100000,
        "filtered": {
          "too_short": 251,
          "too_long": null,
          "too_many_n": null,
          "too_many_expected_errors": null,
          "casava_filtered": null,
          "discard_trimmed": 2061,
          "discard_untrimmed": null
        },
        "output": 97688,
        "reverse_complemented": null,
        "read1_with_adapter": 2254,
        "read2_with_adapter": null
      },
      "basepair_counts": {
        "input": 10100000,
        "input_read1": 10100000,
        "input_read2": null,
        "quality_trimmed": 842048,
        "quality_trimmed_read1": 842048,
        "quality_trimmed_read2": null,
        "output": 9038081,
        "output_read1": 9038081,
        "output_read2": null
      },
      "adapters_read1": [
        {
          "name": "1",
          "type": "regular_three_prime",
          "specification": "ACGTAC",
          "total_matches": 2254,
          "on_reverse_complement": null,
          "five_prime_statistics": null,
          "three_prime_statistics": {
            "error_rate": 0.1,
            "error_lengths": [6],
            "matches": 2254,
            "adjacent_bases": {
              "A": 473,
              "C": 1240,
              "G": 328,
              "T": 207,
              "": 6
            },
            "dominant_adjacent_base": null,
            "trimmed_lengths": [
              [3, 1220, [1220]],
              [4, 319, [319]],
              [5, 30, [30]],
              [6, 4, [4]],
              [7, 5, [5]],
              [8, 7, [7]],
              [9, 4, [4]],
              [10, 7, [7]],
              [11, 7, [7]],
              [12, 6, [6]],
              [13, 8, [8]],
              [14, 1, [1]],
              [15, 2, [2]],
              [16, 3, [3]],
            ]
          }
        }
      ],
      "adapters_read2": null
    }


.. _info-file-format:

Info file format
================

When the ``--info-file`` command-line parameter is given, detailed
information about where adapters were found in each read are written
to the given file. It is a tab-separated text file that contains at
least one row per input read. Normally, there is exactly one row per
input read, but in the following cases, multiple rows may be output:

  * The option ``--times`` is in use.
  * A linked adapter is used.

A row is written for *all* input reads, even those that are discarded
from the final FASTA/FASTQ output due to filtering options.

Which fields are output in each row depends on whether an adapter match was
found in the read or not.

If an adapter match was found, these fields are output in a row:

1. Read name
2. Number of errors
3. 0-based start coordinate of the adapter match
4. 0-based end coordinate of the adapter match
5. Sequence of the read to the left of the adapter match (can be empty)
6. Sequence of the read that was matched to the adapter
7. Sequence of the read to the right of the adapter match (can be empty)
8. Name of the found adapter.
9. Quality values corresponding to sequence left of the adapter match (can be empty)
10. Quality values corresponding to sequence matched to the adapter (can be empty)
11. Quality values corresponding to sequence to the right of the adapter match (can be empty)
12. Flag indicating whether the read was reverse complemented: 1 if yes, 0 if not,
    and empty if ``--revcomp`` was not used.

The concatenation of the fields 5-7 yields the full read sequence. Column 8 identifies
the found adapter. `The section about named adapters <named-adapters>` describes
how to give a name to an adapter. Adapters without a name are numbered starting
from 1. Fields 9-11 are empty if quality values are not available.
Concatenating them yields the full sequence of quality values.

If the adapter match was found on the reverse complement of the read, fields 5 to 7
show the reverse-complemented sequence, and fields 9-11 contain the qualities in
reversed order.

If no adapter was found, the format is as follows:

1. Read name
2. The value -1 (use this to distinguish between match and non-match)
3. The read sequence
4. Quality values

When parsing the file, be aware that additional columns may be added in
the future. Also, some fields can be empty, resulting in
consecutive tabs within a line.

If the ``--times`` option is used and greater than 1, each read can appear
more than once in the info file. There will be one line for each found adapter,
all with identical read names. Only for the first of those lines will the
concatenation of columns 5-7 be identical to the original read sequence (and
accordingly for columns 9-11). For subsequent lines, the shown sequence are the
ones that were used in subsequent rounds of adapter trimming, that is, they get
successively shorter.

Linked adapters appear with up to two rows for each read, one for each constituent
adapter for which a match has been found. To be able to see which of the two
adapters a row describes, the adapter name in column 8 is modified: If the row
describes a match of the 5' adapter, the string ``;1`` is added. If it describes
a match of the 3' adapter, the string ``;2`` is added. If there are two rows, the
5' match always comes first.


.. versionadded:: 1.9
    Columns 9-11 were added.

.. versionadded:: 2.8
    Linked adapters in info files work.

.. versionadded:: 3.4
    Column 12 (revcomp flag) added
