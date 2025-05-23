===============
Reference guide
===============

.. _supported-file-formats:

Supported file formats
======================

Supported input formats are FASTA, FASTQ and unaligned BAM
(uBAM, only for single-end data at the moment).
Supported output formats are FASTA and FASTQ. Compression
:ref:`is supported in multiple formats and detected automatically <compressed-files>`.

The input file format is recognized from the file name extension. If the
extension was not recognized or when Cutadapt reads from standard input,
the contents are inspected instead.

The output file format is also recognized from the file name extension. If the
extensions was not recognized or when Cutadapt writes to standard output, the
same format as the input is used for the output.

When writing a FASTQ file, a second header (the text after the ``+`` on the
third line of a record) that possibly exists in the input is removed.
When writing a FASTA file, line breaks within the sequence are removed.

See also :ref:`file format conversion <file-format-conversion>`.

.. _compressed-files:

Compressed files
----------------

Cutadapt supports compressed input and output files. Whether an input file
needs to be decompressed or an output file needs to be compressed is detected
automatically by inspecting the file name: For example, if it ends in ``.gz``,
then gzip compression is assumed ::

    cutadapt -a AACCGGTT -o output.fastq.gz input.fastq.gz

All of Cutadapt's options that expect a file name support this.

The supported compression formats are gzip (``.gz``), bzip2 (``.bz2``),
xz (``.xz``) and zstandard (``.zst``).

The default compression level for gzip output is 4. Use option ``-Z`` to
change this to level 1. The files need more space, but it is faster and
therefore a good choice for short-lived intermediate files.

Cutadapt uses the `xopen library <https://github.com/pycompression/xopen>`_
to speed up reading and writing compressed files.

Concatenated gzip input files (multiblock gzips) are supported, as well as
concatenated bzip2 files.

.. _standard-input-output:

Standard input and output
-------------------------

If no output file is specified via the ``-o`` option, then the output is sent to
the standard output stream. Example::

    cutadapt -a AACCGGTT input.fastq > output.fastq

There is one difference in behavior if you use Cutadapt without ``-o``: The
report is sent to the standard error stream instead of standard output. You
can redirect it to a file like this::

    cutadapt -a AACCGGTT input.fastq > output.fastq 2> report.txt

Wherever Cutadapt expects a file name, you can also write a dash (``-``) in
order to specify that standard input or output should be used. For example::

    tail -n 4 input.fastq | cutadapt -a AACCGGTT - > output.fastq

The ``tail -n 4`` prints out only the last four lines of ``input.fastq``, which
are then piped into Cutadapt. Thus, Cutadapt will work only on the last read in
the input file.

In most cases, you should probably use ``-`` at most once for an input file and
at most once for an output file, in order not to get mixed output.

For the same reason, you should not use ``-`` for non-interleaved paired-end
data.

You cannot combine ``-`` and gzip compression since Cutadapt needs to know the
file name of the output or input file. if you want to have a gzip-compressed
output file, use ``-o`` with an explicit name.

One last "trick" is to use ``/dev/null`` as an output file name. This special
file discards everything you send into it. If you only want to see the
statistics output, for example, and do not care about the trimmed reads at all,
you could use something like this::

    cutadapt -a AACCGGTT -o /dev/null input.fastq.gz


.. _multicore:

Multi-core support
------------------

Cutadapt supports parallel processing, that is, it can use multiple CPU cores.
Multi-core is not enabled by default. To enable it, use the option ``-j N``
(or the spelled-out version ``--cores=N``), where ``N`` is the
number of cores to use.

To automatically detect the number of available cores, use ``-j 0``
(or ``--cores=0``). The detection takes into account resource restrictions
that may be in place. For example, if running Cutadapt as a batch job on a
cluster system, the actual number of cores assigned to the job will be used.
(This works if the cluster systems uses the cpuset(1) mechanism to impose
the resource limitation.)


.. versionadded:: 1.15

.. versionadded:: 1.18
    ``--cores=0`` for autodetection

.. versionadded:: 2.5
    Multicore works with ``--untrimmed/too-short/too-long-(paired)-output``

.. versionadded:: 2.7
    Multicore works with ``--info-file``, ``--rest-file``, ``--wildcard-file``

.. versionadded:: 3.0
    Multicore support for demultiplexing added.


.. _command-line-options:

Command-line options
====================

General options
---------------

``-h``, ``--help``
    Show help

``--version``
    Show version number and exit

``--debug``
    Print debug log. Use twice to also print the dynamic programming matrices
    computed when aligning an adapter against a read. This is highly verbose,
    it is recommended to use this only for a single read.

``-j CORES``, ``--cores CORES`` (default: 1)
    Run on the :ref:`given number of CPU cores <multicore>`.
    Use 0 to auto-detect the number of available cores.


Adapter-finding options
-----------------------

``-a ADAPTER``, ``--adapter ADAPTER``
    Specification of a :ref:`3' adapter <three-prime-adapters>`
    or a :ref:`linked adapter <linked-adapters>`.

``-g ADAPTER``, ``--front ADAPTER``
    Specification of a :ref:`5' adapter <five-prime-adapters>`
    or a :ref:`linked adapter <linked-adapters>`.

``-b ADAPTER``, ``--anywhere ADAPTER``
    Specification of an adapter that can be :ref:`5' or 3' ("anywhere") <anywhere-adapters>`.

``-e E``, ``--error-rate E``, ``--errors E`` (default: 0.1)
    This sets the :ref:`error tolerance <error-tolerance>` used when searching for adapters.

    If E is an integer >= 1, then E errors in a full-length adapter match are allowed.
    For each specified adapter, this is converted to a maximum allowed error rate.
    This allows proportionally fewer errors for shorter (partial) adapter matches.

    If E is a floating-point value with 0 <= E < 1, this sets the maximum allowed error rate
    directly.

``--no-indels`` (default: allow indels)
    Do not allow insertions and deletions when matching adapters against reads.

``-n COUNT``, ``--times COUNT`` (default: 1)
    Repeat the adapter finding and removal step up to COUNT times.
    :ref:`The default is to search for only one adapter in each read <more-than-one>`.

``-O MINLENGTH``, ``--overlap MINLENGTH`` (default: 3)
    Set the :ref:`minimum overlap <minimum-overlap>` to MINLENGTH.

``--match-read-wildcards``
    Interpret :ref:`IUPAC wildcards in reads <wildcards>` (such as ``N``).

``-N``, ``--no-match-adapter-wildcards``
    Do not interpret :ref:`IUPAC wildcards in adapters <wildcards>`.

``--action {trim,retain,mask,lowercase,none}`` (default: ``trim``)
    Specify what to do if an adapter match was found.

    ``trim``: Trim the adapter itself and up- or downstream sequence (depending on adapter type).

    ``retain``: Trim the up- or downstream sequence (depending on adapter type),
    but retain the adapter sequence itself.

    ``mask``: Replace the adapter sequence and up- or downstream sequence with 'N' characters

    ``lowercase``: Convert the adapter and up- or downstream sequence to lowercase.

    ``none``: Do not change the read. Found matches are still tracked and can be used for
    renaming the read or demultiplexing.

``--rc``, ``--revcomp``
    :ref:`Check both the read and its reverse complement for adapter matches <reverse-complement>`.
    If the reverse-complemented version yields a better match, output that one.

    For paired-end reads, the reverse complement is obtained by swapping R1 and R2.

    If the reverse-complemented version was chosen,
    the string `` rc`` is added to the read name.


Additional read modifications
-----------------------------

.. seealso::

   :ref:`Read modification order <read-modification-order>`

``-u LENGTH``, ``--cut LENGTH``
    :ref:`Remove a fixed number of bases from each read <cut-bases>`.
    If LENGTH is positive, remove bases from the beginning.
    If LENGTH is negative, remove bases from the end.
    Can be used twice if LENGTHs have different signs. This is
    applied *before* adapter trimming.

``-q [5'CUTOFF,]3'CUTOFF``, ``--quality-cutoff [5'CUTOFF,]3'CUTOFF``
    :ref:`Trim low-quality bases <quality-trimming>` from 5' and/or 3' ends of each
    read before adapter removal. This is applied to both reads if
    data is paired (use ``-Q`` to provide a different cutoff for R2).
    If one value is given, only the 3' end
    is trimmed. If two comma-separated cutoffs are given,
    the 5' end is trimmed with the first cutoff, the 3' end
    with the second.

    .. seealso:: :ref:`Description of the quality-trimming algorithm <quality-trimming-algorithm>`

``--nextseq-trim 3'CUTOFF``
    :ref:`NextSeq-specific quality trimming <nextseq-trim>` that
    also trims dark cycles appearing as high-quality G bases.

``--quality-base N`` (default: 33)
    Assume that quality values in FASTQ files are encoded as ascii(quality + N).
    This needs to be set to 64 for some very old Illumina FASTQ files.

``--poly-a``
    :ref:`Trim poly-A tails <poly-A>` from R1 and poly-T heads from R2.

``--length LENGTH``, ``-l LENGTH``
    Shorten reads to LENGTH, where LENGTH is an integer. Positive values remove bases at
    the end while negative ones remove bases at the beginning.

``--trim-n``
    Trim N's from 5' and 3' ends of reads. See: :ref:`Dealing with N bases <n-bases>`.

``--length-tag TAG``
    Search for TAG followed by a decimal number in the header of the FASTQ or FASTA record.
    Replace the decimal number with the correct length of the trimmed read.
    For example, use ``--length-tag 'length='`` to correct fields like 'length=123'.

``--strip-suffix SUFFIX``
    Remove this suffix from read names if present. Can be given multiple times.

``-x PREFIX``, ``--prefix PREFIX``
    Add this prefix to read names. Use ``{name}`` to insert the
    name of the matching adapter. Deprecated, use ``--rename`` instead.

``-y SUFFIX``, ``--suffix SUFFIX``
    Add this suffix to read names. Use ``{name}``` to insert the
    name of the matching adapter. Deprecated, use ``--rename`` instead.

``--rename TEMPLATE``
    :ref:`Rename reads <rename>` using the TEMPLATE, which can contain placeholders such as
    ``{id}``, ``{adapter_name}`` etc.

``--zero-cap``, ``-z``
    Change negative quality values to zero.

Filtering of processed reads
----------------------------

Filters are applied after above read modifications. Paired-end reads are
always discarded pairwise (see also ``--pair-filter``). The default is to not apply any filters.

``-m LEN[:LEN2]``, ``--minimum-length LEN[:LEN2]``
    Discard reads shorter than LEN. If LEN2 is given for paired-end data, it is applied to R2.

``-M LEN[:LEN2]``, ``--maximum-length LEN[:LEN2]``
    Discard reads longer than LEN. If LEN2 is given for paired-end data, it is applied to R2.

``--max-n COUNT``
    Discard reads with more than COUNT 'N' bases.
    If COUNT is a number between 0 and 1,
    it is interpreted as a fraction of the read length. See :ref:`Dealing with N bases <n-bases>`.

``--max-expected-errors E``, ``--max-ee E``
    Discard reads whose :ref:`expected number of errors <expected-errors>` exceeds the value *E*.

``--discard-trimmed``, ``--discard``
    Discard reads in which an adapter match was found.
    Use also ``-O`` to avoid discarding too many randomly matching reads.

``--discard-untrimmed``, ``--trimmed-only``
    Discard reads in which no adapter match was found.

``--discard-casava``
    Discard reads that did not pass CASAVA filtering (that is, the record header has ``:Y:``).

Output
------

``-o FILE``, ``--output FILE``
    Write processed output to FILE (FASTA or FASTQ).
    :ref:`Compressed file formats are supported <compressed-files>`.
    Including the special placeholder string ``{name}`` in the file name activates
    :ref:`demultiplexing`.
    Including ``{name1}`` and ``{name2}`` activates
    :ref:`combinatorial demultiplexing <combinatorial-demultiplexing>`.

    For paired-end data, this option is typically combined with ``-p``.

    If this option is omitted, :ref:`processed reads are written to
    standard output <standard-input-output>`.

``--quiet``
    Print only error messages.

``--report {full,minimal}`` (default: full)
    Which type of report to print: 'full' or 'minimal'.

``--json FILE``
    Write :ref:`a report in JSON format <json-report-format>` to FILE.

``--fasta``
    :ref:`Force writing FASTA to standard output <force-fasta>`.
    This option is usually not needed as FASTA output can be selected by
    using an appropriate output file name (``.fasta``, ``.fasta.gz`` etc.) with the ``-o``
    option. However, when processing FASTQ files *and* when not using ``-o``,
    FASTQ format is written to standard output by default.
    Use this option to force FASTA even in such a case.

``-Z``
    **Deprecated**: This option has become the default in Cutadapt 4.10.

    Use compression level 1 for gzipped output files.
    This is a shorthand for ``--compression-level=1``.

    See: :ref:`speed-up tricks <speedup>`

``--info-file FILE``
    Write information about each read and its adapter matches to FILE.
    See: :ref:`Info file format <info-file-format>`.

``-r FILE``, ``--rest-file FILE``
    When the adapter matches in the middle of a read, write the "rest" to FILE.
    For 3' adapters, the "rest" is the part of the read after the adapter match.
    For 5' adapters, the "rest" is the part of the read before the adapter match.

``--wildcard-file FILE``
    When the adapter has N wildcard bases, write adapter bases matching wildcard positions to FILE.
    This is unreliable unless you also use ``--noindels``.
    Does not work with linked adapters.

``--too-short-output FILE``
    Write reads that are too short (according to the length specified by ``-m``) to FILE.
    Default: discard too short reads

``--too-long-output FILE``
    Write reads that are too long (according to length specified by -M) to FILE.
    Default: discard too long reads

``--untrimmed-output FILE``
    Write reads that do not contain any adapter to FILE.
    Default: output to the same file as trimmed reads.

Paired-end options
------------------

.. seealso:: :ref:`Trimming paired-end reads <paired-end>`

The ``-A``, ``-G``, ``-B``, ``-U``, ``-Q`` options work like their lowercase counterparts,
but are applied to the second read in each pair (R2).

``-A ADAPTER``
    3' adapter to be removed from R2

``-G ADAPTER``
    5' adapter to be removed from R2

``-B ADAPTER``
    5'/3 adapter to be removed from R2

``-U LENGTH``
    Remove LENGTH bases from R2

``-Q [5'CUTOFF,]3'CUTOFF``
    Quality-trimming cutoff for R2. Default: same as for R1

..
      -p FILE, --paired-output FILE
                            Write R2 to FILE.
      --pair-adapters       Treat adapters given with -a/-A etc. as pairs. Either
                            both or none are removed from each read pair.
      --pair-filter {any,both,first}
                            Which of the reads in a paired-end read have to match
                            the filtering criterion in order for the pair to be
                            filtered. Default: any
      --interleaved         Read and/or write interleaved paired-end reads.
      --untrimmed-paired-output FILE
                            Write second read in a pair to this FILE when no adapter
                            was found. Use with --untrimmed-output. Default: output
                            to same file as trimmed reads
      --too-short-paired-output FILE
                            Write second read in a pair to this file if pair is too
                            short.
      --too-long-paired-output FILE
                            Write second read in a pair to this file if pair is too
                            long.


(To Do: needs to be finished, see ``cutadapt --help`` for now)

.. _json-report-format:

JSON report format
==================

The JSON reported is generated if ``--json=filename.cutadapt.json`` is used. The file name
extension must be ``.cutadapt.json`` for the file to be recognized by log-parsing tools such
as `MultiQC <https://multiqc.info>`_. (However, at the time of writing, MultiQC does not support
Cutadapt’s JSON report format.)

See how to :ref:`extract information from the JSON report with jq <json-jq>`.

Example
-------

This example was reformatted to use less vertical space::

    {
      "tag": "Cutadapt report",
      "schema_version": [0, 3],
      "cutadapt_version": "4.5",
      "python_version": "3.8.10",
      "command_line_arguments": [
        "--json=out.cutadapt.json", "--poly-a", "-m", "20", "-a", "AACCGGTTACGTTGCA",
        "-q", "20", "--discard-trimmed", "-o", "out.fastq.gz", "reads.fastq"],
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
        "poly_a_trimmed": 1028,
        "poly_a_trimmed_read1": 1028,
        "poly_a_trimmed_read2": null,
        "output": 9037053,
        "output_read1": 9037053,
        "output_read2": null
      },
      "adapters_read1": [
        {
          "name": "1",
          "total_matches": 2254,
          "on_reverse_complement": null,
          "linked": false,
          "five_prime_end": null,
          "three_prime_end": {
            "type": "regular_three_prime",
            "sequence": "AACCGGTTACGTTGCA",
            "error_rate": 0.1,
            "indels": true,
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
              {"len": 3, "expect": 1562.5, "counts": [1220]},
              {"len": 4, "expect": 390.6, "counts": [319]},
              {"len": 5, "expect": 97.7, "counts": [30]},
              {"len": 6, "expect": 24.4, "counts": [4]},
              {"len": 7, "expect": 24.4, "counts": [5]},
              {"len": 8, "expect": 24.4, "counts": [7]},
              {"len": 9, "expect": 24.4, "counts": [4]},
              {"len": 10, "expect": 24.4, "counts": [7]},
              {"len": 11, "expect": 24.4, "counts": [7]},
              {"len": 12, "expect": 24.4, "counts": [6]},
              {"len": 13, "expect": 24.4, "counts": [8, 2]},
              {"len": 14, "expect": 24.4, "counts": [1, 1]},
              {"len": 15, "expect": 24.4, "counts": [2, 0]},
              {"len": 16, "expect": 24.4, "counts": [3, 1]},
            ]
          }
        }
      ],
      "adapters_read2": null,
      "poly_a_trimmed_read1": [
        {"len": 23, "count": 10},
        {"len": 42, "count": 19}
      ],
      "poly_a_trimmed_read2": null
    }


Schema
------

Some concepts used in the JSON file:

* Keys are always included. If a key is not applicable, its value is set to null.
* Single-end data appears as "paired-end data without read 2". That is, values for
  read 1 are filled in and values for read 2 are set to null.

The file defines the following keys. For nested objects (dictionaries), a dot notation is used,
as in "outer_key.inner_key".

tag : string
   Always ``"Cutadapt report"``. A marker so that this can be recognized as a file produced by
   Cutadapt.

schema_version : list of two integers
   Major and minor version of the schema.
   If additions are made to the schema, the minor version is increased.
   If backwards incompatible changes are made, the major version is increased.

   Example: ``[0, 1]``

cutadapt_version : str
   The version of Cutadapt that generated the report.

   Example: ``"4.4"``

python_version : str
   The Python version used to run Cutadapt.

   Example: ``"3.10"``

command_line_arguments : list of strings
   The command-line arguments for this invocation. Only for information, do not parse this.

   Example: ``["-a", "ACGT", "-o", "out.fastq", "input.fastq"]```

cores : int
   Number of cores used

input : dictionary
   Input files

input.path1 : str
   Path to the first input file.

   Example: ``"reads.1.fastq"``

input.path2 : str | null
   Path to the second input file if given, null otherwise.

input.paired : bool
   True if input was paired-end reads, false if input was single-end reads.
   If this is true and input.path2 is null, input was interleaved.

read_counts : dictionary
   Read count statistics

read_counts.input : int
   Number of reads (for single-end data) or read pairs (for paired-end data) in the input.

read_counts.filtered : dictionary
   Statistics about filtered reads. Keys of the dictionary correspond to a filter.
   If a filter was not used, its value is set to null.

read_counts.filtered.too_short : int | null
   Number of reads or read pairs that were filtered because they were too short

read_counts.filtered.too_long : int | null
   Number of reads or read pairs that were filtered because they were too long

read_counts.filtered.too_many_n : int | null
   Number of reads or read pairs that were filtered because they had too many N bases

read_counts.filtered.too_many_expected_errors : int | null
   Number of reads or read pairs that were filtered because they had too many expected errors

read_counts.filtered.casava_filtered : int | null
   Number of reads or read pairs that were filtered because the CASAVA filter was ``Y``

read_counts.filtered.discard_trimmed : int | null
   Number of reads or read pairs that were filtered because at least one adapter match was found for them

read_counts.filtered.discard_untrimmed : int | null
   Number of reads or read pairs that were filtered because no adapter match was found for them

read_counts.output : int
   Number of reads written to the final output.
   This plus the sum of all filtered reads/read will equal the number of input reads.

read_counts.reverse_complemented : int | null
   If ``--revcomp`` was used, the number of reads or read pairs that were output
   reverse-complemented, null otherwise.

read_counts.read1_with_adapter : int | null
   Number of R1 reads (or single-end reads) with at least one adapter match,
   null if no adapter trimming was done.

read_counts.read2_with_adapter : int | null
   Number of R2 reads with at least one adapter match, null if input is single end or no
   adapter trimming was done.

basepair_counts : dictionary
   Statistics about the number of basepairs.

basepair_counts.input : int
   Total number of basepairs in the input. (The sum of the lengths of all input reads.)

basepair_counts.input_read1 : int
   Number of basepairs in the input, read 1 only.

basepair_counts.input_read2 : int | null
   If paired-end, number of basepairs in the input counting read 2 only, null otherwise.

basepair_counts.quality_trimmed : int | null
   Total number of basepairs removed due to quality trimming or null if no quality trimming was done.

basepair_counts.quality_trimmed_read1 : int | null
   Number of basepairs removed from read 1 due to quality trimming or null if no quality trimming
   was done.

basepair_counts.quality_trimmed_read2 : int
   Number of basepairs removed from read 2 due to quality trimming or null if no quality trimming was
   done or if input was single end.

basepair_counts.poly_a_trimmed : int | null
   Total number of basepairs removed due to poly-A trimming or null if no poly-A trimming was done.

basepair_counts.poly_a_trimmed_read1 : int | null
   Number of basepairs removed from read 1 due to poly-A trimming or null if no poly-A trimming was
   done.

basepair_counts.poly_a_trimmed_read2 : int
   Number of basepairs removed from read 2 due to poly-T trimming or null if no poly-T trimming was
   done or if input was single end.

basepair_counts.output : int
   Total number of basepairs in the final output.

basepair_counts.output_read1 : int
   Number of basepairs written to the read 1 final output.

basepair_counts.output_read2 : int | null
   Number of basepairs written to the read 2 final output.

adapters_read1 : list of dictionaries
   A list with statistics about all adapters that were matched against read 1.
   The list is empty if no adapter trimming was done. The schema for the items in this list is
   described below.

adapters_read2 : list of dictionaries | null
   A list with statistics about all adapters that were matched against read 2.
   The list is empty if no adapter trimming was done against R2. The value is set to null if
   the input was single end reads. The schema for the items in this list is described below.

poly_a_trimmed_read1 : list of dictionaries | null
   A histogram of the lengths of poly-A tails removed from read 1. Each item in the list is a
   dictionary with keys ``len`` and ``count``. This value is null if no poly-A trimming was done.

poly_a_trimmed_read2 : list of dictionaries | null
   A histogram of the lengths of poly-T "heads" removed from read 2, see above. This value is null
   if no poly-A/poly-T trimming was done or the input was single-end reads.


Adapter statistics
------------------

The statistics about each adapter (items in the adapters_read1 and adapters_read2 list) are
dictionaries with the following keys.

name : str
   The adapter name. If no adapter name was given, a name is automatically generated as
   "1", "2", "3" etc.

total_matches : int
   Number of times this adapter was found on a read. If ``--times`` is used, multiple matches
   per read are possible.

on_reverse_complement : int | null
   If ``--revcomp`` was used, the number of times the adapter was found on the reverse-complemented
   read, null otherwise.

linked : bool
   Whether this is a linked adapter. If true, then both ``five_prime_end`` and ``three_prime_end``
   (below) are filled in and describe the 5' and 3' components, respectively, of the linked adapter.

five_prime_end : dictionary | null
   Statistics about matches of this adapter to the 5' end, that is, causing a prefix of the
   read to be removed.

   If the adapter is of type regular_five_prime, noninternal_five_prime or anchored_five_prime,
   all its matches are summarized here.

   If the adapter is a linked adapter (``linked`` is true), the matches of its 5' component are
   summarized here.

   If the adapter is of type "anywhere", the matches that were determined to be 5' matches are
   summarized here.

   This is null for the other adapter types.

three_prime_end : dictionary | null
   Statistics about matches of this adapter to the 3' end, that is, causing a suffix of the read
   to be removed.

   If the adapter is of type regular_three_prime, noninternal_three_prime or anchored_three_prime,
   all its matches are summarized here.

   If the adapter is a linked adapter (``linked`` is true), the matches of its 3' component are
   summarized here.

   If the adapter is of type "anywhere", the matches that were determined to be 3' matches are
   summarized here.

   This is null for the other adapter types.

three/five_prime_end.type : str
   Type of the adapter. One of these strings:
     - ``"regular_five_prime"``
     - ``"regular_three_prime"``
     - ``"noninternal_five_prime"``
     - ``"noninternal_three_prime"``
     - ``"anchored_five_prime"``
     - ``"anchored_three_prime"``
     - ``"anywhere"``

   For linked adapters, this is the type of its 5' or 3' component.

three/five_prime_end.sequence : str
   Sequence of this adapter. For linked adapters, this is the sequence of its 5' or 3' component.

   Example: ``"AACCGGTT"``

three/five_prime_end.error_rate : float
   Error rate for this adapter. For linked adapters, the error rate for the respective end.

three/five_prime_end.indels : bool
   Whether indels are allowed when matching this adapter against the read.

three/five_prime_end.error_lengths : list of ints
   If the adapter type allows partial matches, this lists the lengths up to which 0, 1, 2 etc.
   errors are allowed. Example: ``[9, 16]`` means: 0 errors allowed up to a match of length 9,
   1 error up to a match of length 16. The last number in this list is the length of the adapter
   sequence.

   For anchored adapter types, this is null.

three/five_prime_end.matches : int
   The number of matches of this adapter against the 5' or 3' end.

three/five_prime_end.adjacent_bases : dictionary | null
   For 3' adapter types, this shows which bases occurred adjacent to (upstream of) the 3' adapter
   match. It is a dictionary mapping the strings "A", "C", "G", "T" and "" (empty string) to
   the number of occurrences. The empty string covers those cases in which the adjacent base
   was not one of A, C, G or T or in which there was no adjacent base (3' adapter found at the
   beginning of the read).

   This is null for 5' adapters (adjacent base statistics are currently not tracked for those).

three/five_prime_end.dominant_adjacent_base : str | null
   This is set to the dominant adjacent base if adjacent_bases exist and were determined to be
   sufficiently skewed, corresponding to the :ref:`warning <warnbase>`:
   "The adapter is preceded by "x" extremely often."

   This is null otherwise.

three/five_prime_end.trimmed_lengths : list of dictionaries
   The histogram of the lengths of removed sequences. Each item in the list is a dictionary
   that describes how often a sequence of a certain length was removed,
   broken down by the number of errors in the adapter match.

   Example::

      "trimmed_lengths": [
        {"len": 4, "expect": 390.6, "counts": [319]},
        {"len": 5, "expect": 97.7, "counts": [30]},
        {"len": 6, "expect": 24.4, "counts": [4]},
        {"len": 7, "expect": 24.4, "counts": [5]},
        {"len": 15, "expect": 24.4, "counts": [2, 1]},
      ]

three/five_prime_end.trimmed_lengths.expect : float
   How often a sequence of length *len* would be expected to be removed due to random chance.

three/five_prime_end.trimmed_lengths.counts : list of int
   Element at index *i* in this list gives how often a sequence of length *len* was removed due to
   an adapter match with *i* errors. Sum these values to get the total count.

   Example (5 sequences had 0 errors in the adapter matches, 3 had 1 and 1 had 2)::

   [5, 3, 1]


.. _info-file-format:

Info file format
================

When the ``--info-file`` command-line parameter is given, detailed
information about where adapters were found in each read are written
to the given file. It is a tab-separated text file that contains at
least one row per input read. Normally, there is exactly one row per
input read, but in the following cases, multiple rows may be output:

 - The option ``--times`` is in use.
 - A linked adapter is used.

A row is written for *all* input reads, even those that are discarded
from the final FASTA/FASTQ output due to filtering options.

.. note:: Paired-end reads are not supported.
   The info file currently does not contain any info about read 2 when
   Cutadapt is run in paired-end mode.

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


.. _properly-paired-reads:

Properly paired reads
=====================

When reading paired-end reads, Cutadapt compares the read IDs of R1 and R2.
It prints an error message and aborts if they do not match.

Comments in the FASTQ or FASTA header are ignored when doing the comparison.
Also, if the read ID ends with ``1`` or ``2`` or ``3``, then that is also
ignored.

For example, two FASTQ headers that would be considered to denote properly paired reads are::

    @my_read/1 a comment

and::

    @my_read/2 another comment

This is an example for *improperly paired* reads::

    @my_read/1;1

and::

    @my_read/2;1

Since the ``1`` and ``2`` are ignored only if they occur at the end of the read
name, and since the ``;1`` is considered to be part of the read name, these
reads will not be considered to be properly paired.


.. _read-modification-order:

Read modification order
=======================

Read modifications are applied in the following order to each read.
Steps not requested on the command-line are skipped.

1. Unconditional base removal with ``--cut``
2. Quality trimming (``-q``)
3. Adapter trimming (``-a``, ``-b``, ``-g`` and uppercase versions)
4. Poly-A/poly-T trimming (``--poly-a``)
5. Read shortening (``--length``)
6. N-end trimming (``--trim-n``)
7. Length tag modification (``--length-tag``)
8. Read name suffix removal (``--strip-suffix``)
9. Addition of prefix and suffix to read name (``-x``/``--prefix`` and ``-y``/``--suffix``)
10. Read renaming according to ``--rename``
11. Replacing of negative quality values with zero (zero capping)
