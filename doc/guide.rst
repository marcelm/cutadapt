cutadapt
========

cutadapt removes adapter sequences from DNA high-throughput sequencing
data. This is necessary when the reads are longer than the molecule that
is sequenced, such as in microRNA data.

cutadapt is implemented in Python. It comes with an extension module
written in Cython that implements the alignment algorithm.

License
-------

cutadapt is licensed under the MIT license (see the file ``LICENSE``).

Project homepage
----------------

See http://code.google.com/p/cutadapt/ . Please use the Google code
issue tracker for bug reports and feature requests.

Galaxy
------

If you want to use cutadapt within the web-based Galaxy platform
(http://galaxy.psu.edu/), please see the README file in the galaxy/
subfolder. Galaxy support was contributed by Lance Parsons.

Basic usage
===========

The basic command-line for cutadapt is::

    cutadapt -a AACCGGTT input.fastq > output.fastq

The adapter sequence is given with the ``-a`` option. Replace
``AACCGGTT`` with your actual adapter sequence.

``input.fastq`` is a file with reads. The result will be written to
standard output. Use redirection with ``>`` (or the ``-o`` option) to
write the output to a file.

cutadapt prints out a report after it has finished processing the reads.
If you use ``-o``, the report is sent to standard output and to stderr
otherwise.

By default, the output file contains all reads, including those that did
not contain an adapter. (See also the ``--discard`` option.)

The following examples refer to basespace reads. See the "Colorspace"
section on how to use cutadapt with SOLiD reads.

Only a few command-line options are explained in this document. To see
all options, run::

    cutadapt --help

In particular, see the explanation for the different types of adapters
that are supported.

Trimming FASTQ files
--------------------

Assuming your sequencing data is available as a FASTQ file, use this
command line::

    cutadapt -a ADAPTER-SEQUENCE input.fastq > output.fastq

gz-compressed input is supported::

    cutadapt -a ADAPTER-SEQUENCE input.fastq.gz > output.fastq

gz-compressed output is also supported, but the -o parameter (output
file) needs to be used since gzip compression is auto-detected by inspecting
the file name::

    cutadapt -a ADAPTER-SEQUENCE -o output.fastq.gz input.fastq.gz

If your Python installation includes support for bzip2 compression, then
bzip2-compressed files are also supported and recognized by their
extension ``.bz2``.

Adapters in FASTA files
-----------------------

To read a list of adapter sequences from a FASTA file, specify the file
in the following way::

    cutadapt -a file:adapters.fasta input.fastq > output.fastq

All of the sequences in the file ``adapters.fasta`` will be used as 3'
adapters. As always, only the best matching adapter will be trimmed from
each read.

Named adapters
--------------

You can give names to the adapters. The names are shown in addition to
the sequences themselves in the statistics overview when the program has
finished trimming the reads. You can use it like this::

    cutadapt -a My_adapter=ACGTAA input.fastq > output.fastq

Here, the actual adapter sequence is ``ACGTAA`` and the name assigned to
it is ``My_adapter``. When adapters are read from a FASTA file, the
sequence header is used as the adapter name.

Wildcards
---------

The wildcard character ``N`` in the adapter sequence is supported. It matches
any nucleotide. This is useful for trimming adapters that have a variable
barcode embedded in them::

    cutadapt -a ACGTAANNNNTTAGC -o output.fastq input.fastq

Wildcard characters in the reads are also supported, but this must be
enabled with ``--match-read-wildcards``.

FASTA file
----------

Cut an adapter from reads given in a FASTA file. Try to remove an
adapter three times (this is usually not needed), use the default error
rate of 10%, write result to ``output.fa``::

    cutadapt -n 3 -a TGAGACACGCAACAGGGGAAAGGCAAGGCACACAGGGGATAGG input.fa > output.fa

Multiple adapters
-----------------

As many adapters as desired can be given to the program by using the
``-a``, ``-b`` or ``-g`` in any combination, for example, five ``-a``
adapters and two ``-g`` adapters. All adapters will be searched for, but
only the best matching one will be trimmed from each read (but see the
``--times`` option)::

    cutadapt -b TGAGACACGCA -g AGGCACACAGGG input.fastq > output.fastq

Quality trimming
----------------

The ``-q`` (or ``--trim-qualities``) parameter can be used to trim
low-quality ends from reads before adapter removal. For this to work
correctly, the quality values must be encoded as ascii(phred quality +
33). If they are encoded as ascii(phred quality + 64), you need to add
``--quality-base=64`` to the command line.

The trimming algorithm is the same as the one used by BWA. That is:
Subtract the given cutoff from all qualities; compute partial sums from
all indices to the end of the sequence; cut sequence at the index at
which the sum is minimal.

Removing bases from the beginning or end of each read
-----------------------------------------------------

By using the ``--cut`` or its abbreviation ``-u``, it is possible to
unconditionally remove bases from the beginning or end of each read. If
the given length is positive, the bases are removed from the beginning
of each read. If it is negative, the bases are removed from the end.

Remove the first seven bases of each read::

    cutadapt -u 7 -o trimmed.fastq reads.fastq

Remove the last seven bases of each read::

    cutadapt -u -7 -o trimmed.fastq reads.fastq

The ``-u``/``--cut`` option can be combined with the other options, but
the desired bases are removed *before* any adapter trimming.

Paired-end adapter trimming
---------------------------

Cutadapt supports paired-end trimming, but currently two passes over the
data are required.

Assume the input is in ``reads.1.fastq`` and ``reads.2.fastq`` and that
``ADAPTER_FWD`` should be trimmed from the forward reads (first file)
and ``ADAPTER_REV`` from the reverse reads (second file).

If you do not use any of the filtering options that discard reads, such
as ``--discard``, ``--minimum-length`` or ``--maximum-length``, then run
cutadapt on each file separately::

    cutadapt -a ADAPTER_FWD -o trimmed.1.fastq reads1.fastq
    cutadapt -a ADAPTER_REV -o trimmed.2.fastq reads2.fastq

You can use the options that are listed under 'Additional modifications'
in cutadapt's help output without problems. For example, if you want to
quality-trim the first read in each pair with a threshold of 10, and the
second read in each pair with a threshold of 15, then the commands could
be::

    cutadapt -q 10 -a ADAPTER_FWD -o trimmed.1.fastq reads1.fastq
    cutadapt -q 15 -a ADAPTER_REV -o trimmed.2.fastq reads2.fastq

However, if you use one of the filtering options that discard reads,
then you need to give both input read files to cutadapt and the
``--paired-output`` option is needed to keep the two files synchronized.
First trim the forward read, writing output to temporary files (we also
add some quality trimming)::

    cutadapt -q 10 -a ADAPTER_FWD --minimum-length 20 -o tmp.1.fastq -p tmp.2.fastq reads.1.fastq reads.2.fastq

The ``-p`` is an abbreviation for ``--paired-output``. Then trim the
reverse read, using the temporary files as input::

    cutadapt -q 15 -a ADAPTER_REV --minimum-length 20 -o trimmed.2.fastq -p trimmed.1.fastq tmp.2.fastq tmp.1.fastq

Finally, remove the temporary files::

    rm tmp.1.fastq tmp.2.fastq

In each call to cutadapt, the read-modifying options such as ``-q`` only
apply to the first file (first ``reads.1.fastq``, then ``tmp.2.fastq``
in this example). Reads in the second file are not affected by those
options, but by the filtering options: If a read in the first file is
discarded, then the matching read in the second file is also filtered
and not written to the output given by ``--paired-output`` in order to
keep both output files synchronized.

When you use ``-p``/``--paired-output``, then cutadapt also checks
whether the files are properly paired. An error is raised if one of the
files contains more reads than the other or if the read names in the two
files do not match. Only the part of the read name before the first
space is considered. If the read name ends with ``/1`` or ``/2``, then
that is also ignored. For example, two FASTQ headers that would be
considered to denote properly paired reads are::

    @my_read/1 a comment

and::

    @my_read/2 another comment

Illumina TruSeq
---------------

If you have reads containing Illumina TruSeq adapters, follow these
steps.

Trim read 1 with ``A`` + the “TruSeq Indexed Adapter”. Use only the
prefix of the adapter sequence that is common to all Indexed Adapter
sequences::

    cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o trimmed.1.fastq.gz reads.1.fastq.gz

Trim read 2 with the reverse complement of the “TruSeq Universal
Adapter”::

    cutadapt -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o trimmed.2.fastq.gz reads.2.fastq.gz

See also the section about paired-end adapter trimming above.

If you want to simplify this a bit, you can also use ``AGATCGGAAGAGC``
as the adapter sequence in both cases::

    cutadapt -a AGATCGGAAGAGC -o trimmed.1.fastq.gz reads.1.fastq.gz
    cutadapt -a AGATCGGAAGAGC -o trimmed.2.fastq.gz reads.2.fastq.gz

The adapter sequences can be found in the document `Illumina TruSeq
Adapters
De-Mystified <http://tucf-genomics.tufts.edu/documents/protocols/TUCF_Understanding_Illumina_TruSeq_Adapters.pdf>`__.

Adapters
========

These are some 454 adapters::

    A1:   5'- TCCATCTCATCCCTGCGTGTCCCATCTGTTCCCTCCCTGTCTCA
    A2:   5'- TGAGACAGGGAGGGAACAGATGGGACACGCAGGGATGAGATGGA
    B1:   5'- CCTATCCCCTGTGTGCCTTGCCTATCCCCTGTTGCGTGTCTCA
    B2:   5'- TGAGACACGCAACAGGGGAAAGGCAAGGCACACAGGGGATAGG

This is an AB SOLiD adapter (in color space) used in the SREK protocol::

    330201030313112312

Algorithm
=========

cutadapt uses a modified semi-global alignment algorithm. For speed, the
algorithm is implemented as a Cython extension module in ``_align.pyx``.

Cutadapt’s processing speed is currently not dominated by the alignment
algorithm, but by parsing the input and writing the output.

Partial adapter matches
-----------------------

Cutadapt correctly deals with partial adapter matches, and also with any
trailing sequences after the adapter. As an example, suppose your
adapter sequence is "ADAPTER" (specified via the ``-a`` or ``--adapter``
command-line parameter). If you have these input sequences::

    MYSEQUENCEADAPTER
    MYSEQUENCEADAP
    MYSEQUENCEADAPTERSOMETHINGELSE

All of them will be trimmed to "MYSEQUENCE". If the sequence starts with
an adapter, like this::

    ADAPTERSOMETHING

It will be empty after trimming.

When the allowed error rate is sufficiently high (set with parameter
``-e``), errors in the adapter sequence are allowed. For example,
``ADABTER`` (1 mismatch), ``ADAPTR`` (1 deletion), and ``ADAPPTER`` (1
insertion) will all be recognized if the error rate is set to 0.15.

Anchoring 5' adapters
---------------------

If you specify an adapter with the ``-g`` (``--front``) parameter, the
adapter may overlap the beginning of the read or occur anywhere within
it. If it appears within the read, the sequence that precedes it will
also be trimmed in addition to the adapter. For example, with
``-g ADAPTER``, these sequences::

    HELLOADAPTERTHERE
    APTERTHERE

will both be trimmed to ``THERE``. To avoid this, you can prefix the
adapter with the character ``^``. This will restrict the search, forcing
the adapter to be a prefix of the read. With ``-g ^ADAPTER``, only reads
like this will be trimmed::

    ADAPTERHELLO

Allowing adapters anywhere
--------------------------

Cutadapt assumes that any adapter specified via the ``-a`` (or
``--adapter``) parameter was ligated to the 3' end of the sequence. This
is the correct assumption for at least the SOLiD and Illumina small RNA
protocols and probably others. The assumption is enforced by the
alignment algorithm, which only finds the adapter when its starting
position is within the read. In other words, the 5' base of the adapter
must appear within the read. The adapter and all bases following it are
removed.

If, on the other hand, your adapter can also be ligated to the 5' end
(on purpose or by accident), you should tell cutadapt so by using the
``-b`` (or ``--anywhere``) parameter. It will then use a slightly
different alignment algorithm (so-called semiglobal alignment), which
allows any type of overlap between the adapter and the sequence. In
particular, the adapter may appear only partially in the beginning of
the read, like this::

    PTERMYSEQUENCE

The decision which part of the read to remove is made as follows: If
there is at least one base before the found adapter, then the adapter is
considered to be a 3' adapter and the adapter itself and everything
following it is removed. Otherwise, the adapter is considered to be a 5'
adapter and it is removed from the read.

Here are some examples, which may make this clearer (left: read, right:
trimmed read)::

    MYSEQUENCEADAPTER -> MYSEQUENCE (3' adapter)
    MADAPTER -> M (3' adapter)
    ADAPTERMYSEQUENCE -> MYSEQUENCE (5' adapter)
    PTERMYSEQUENCE -> MYSEQUENCE (5' adapter)

The regular algorithm (``-a``) would trim the first two examples in the
same way, but trim the third to an empty sequence and trim the fourth
not at all.

The ``-b`` parameter currently does not work with color space data.

Interpreting the statistics output
==================================

After every run, cutadapt prints out per-adapter statistics. The output
starts with something like this::

    Adapter 'ACGTACGTACGTTAGCTAGC', length 20, was trimmed 2402 times.

The meaning of this should be obvious.

The next piece of information is this::

    No. of allowed errors:
    0-9 bp: 0; 10-19 bp: 1; 20 bp: 2

The adapter has, as was conveniently shown above, a length of 20
characters. We are using the default error rate of 0.1. What this
implies is shown above: Matches up to a length of 9 bp are allowed to
have no errors. Matches of lengths 10-19 bp are allowd to have 1 error
and matches of length 20 can have 2 errors.

Finally, a table is output that gives more detailed information about
the lengths of the removed sequences. The following is only an excerpt;
some rows are left out::

    Overview of removed sequences
    length  count   expect  max.err error counts
    3       140     156.2   0       140
    4       57      39.1    0       57
    5       50      9.8     0       50
    6       35      2.4     0       35
    ...
    100     397     0.0     3       358 36 3

The first row tells us the following: Three bases were removed in 140
reads; randomly, one would expect this to occur 156.2 times; the maximum
number of errors at that match length is 0 (this is actually redundant
since we know already that no errors are allowed at lengths 0-9bp).

The last column shows the number of reads that had 0, 1, 2 ... errors.
In the last row, for example, 358 reads matched the adapter with zero
errors, 36 with 1 error, and 3 matched with 2 errors.

The "expect" column gives only a rough estimate of the number of
sequences that is expected to match randomly (it assumes a GC content of
50%, for example), but it can help to estimate whether the matches that
were found are true adapter matches or if they are due to chance. At
lengths 6, for example, only 2.4 reads are expected, but 35 do match,
which hints that most of these matches are due to actual adapters.

Note that the "length" column refers to the length of the removed
sequence. That is, the actual length of the match in the above row at
length 100 is 20 since that is the adapter length. Assuming the read
length is 100, the adapter was found in the beginning of 397 reads and
therefore those reads were trimmed to a length of zero.

The table may also be useful in case the given adapter sequence contains
an error. In that case, it may look like this::

    ...
    length  count   expect  max.err error counts
    10      53      0.0     1       51 2
    11      45      0.0     1       42 3
    12      51      0.0     1       48 3
    13      39      0.0     1       0 39
    14      40      0.0     1       0 40
    15      36      0.0     1       0 36
    ...

We can see that no matches longer than 12 have zero errors. In this
case, it indicates that the 13th base of the given adapter sequence is
incorrect.

Format of the info file
=======================

When the ``--info-file`` command-line parameter is given, detailed
information about the found adapters is written to the given file. The
output is a tab-separated text file. Each line corresponds to one read
of the input file. The fields are:

1. Read name
2. Number of errors
3. 0-based start coordinate of the adapter match
4. 0-based end coordinate of the adapter match
5. Sequence of the read to the left of the adapter match (can be empty)
6. Sequence of the read that was matched to the adapter
7. Sequence of the read to the right of the adapter match (can be empty)
8. Name of the found adapter.

The concatenation of the fields 5-7 yields the full read sequence. The
adapter name for column 8 can be given by writing ``-a name=sequence``
instead of just ``-a sequence``. Adapters without a name are numbered
starting from 1.

If no adapter was found, the format is as follows:

-  Read name
-  The value -1
-  The read sequence

When parsing that file, be aware that additional columns may be added in
the future. Note also that some fields can be empty, resulting in
consecutive tabs within a line. Also, in the current version, when the
``--times`` option is set to a value other than 1 (the default value),
multiple lines are written to the info file for each read.

Colorspace
==========

Cutadapt was designed to work with colorspace reads from the ABi SOLiD
sequencer. Colorspace trimming is activated by the ``--colorspace``
option (or use ``-c`` for short). The input reads can be given either:

-  in a FASTA file
-  in a FASTQ file
-  in a ``.csfasta`` and a ``.qual`` file (this is the native SOLiD
   format).

In all cases, the colors must be represented by the characters 0, 1, 2,
3. Example input files are in the cutadapt distribution at
``tests/data/solid.*``. The ``.csfasta``/``.qual`` file format is
automatically assumed if two input files are given to cutadapt.

In colorspace mode, the adapter sequences given to the ``-a``, ``-b``
and ``-g`` options can be given both as colors or as nucleotides. If
given as nucleotides, they will automatically be converted to
colorspace. For example, to trim an adapter from ``solid.csfasta`` and
``solid.qual``, use this command-line::

    cutadapt -c -a CGCCTTGGCCGTACAGCAG solid.csfasta solid.qual > output.fastq

In case you know the colorspace adapter sequence, you can also write
``330201030313112312`` instead of ``CGCCTTGGCCGTACAGCAG`` and the result
is the same.

Ambiguity in colorspace
-----------------------

The ambiguity of colorspace encoding leads to some effects to be aware
of when trimming 3' adapters from colorspace reads. For example, when
trimming the adapter ``AACTC``, cutadapt searches for its
colorspace-encoded version ``0122``. But also ``TTGAG``, ``CCAGA`` and
``GGTCT`` have an encoding of ``0122``. This means that effectively four
different adapter sequences are searched and trimmed at the same time.
There is no way around this, unless the decoded sequence were available,
but that is usually only the case after read mapping.

The effect should usually be quite small. The number of false positives
is multiplied by four, but with a sufficiently large overlap (3 or 4 is
already enough), this is still only around 0.2 bases lost per read on
average. If inspecting k-mer frequencies or using small overlaps, you
need to be aware of the effect, however.

Double-encoding, BWA and MAQ
----------------------------

The read mappers MAQ and BWA (and possibly others) need their colorspace
input reads to be in a so-called "double encoding". This simply means
that they cannot deal with the characters 0, 1, 2, 3 in the reads, but
require that the letters A, C, G, T be used for colors. For example, the
colorspace sequence ``0011321`` would be ``AACCTGC`` in double-encoded
form. This is not the same as conversion to basespace! The read is still
in colorspace, only letters are used instead of digits. If that sounds
confusing, that is because it is.

Note that MAQ is unmaintained and should not be used in new projects.

BWA’s colorspace support was dropped in versions more recent than 0.5.9,
but that version works well.

When you want to trim reads that will be mapped with BWA or MAQ, you can
use the ``--bwa`` option, which enables colorspace mode (``-c``),
double-encoding (``-d``), primer trimming (``-t``), all of which are
required for BWA, in addition to some other useful options.

The ``--maq`` option is an alias for ``--bwa``.

Colorspace examples
-------------------

To cut an adapter from SOLiD data given in ``solid.csfasta`` and
``solid.qual``, to produce MAQ- and BWA-compatible output, allow the
default of 10% errors and write the resulting FASTQ file to
output.fastq::

    cutadapt --bwa -a CGCCTTGGCCGTACAGCAG solid.csfasta solid.qual > output.fastq

Instead of redirecting standard output with ``>``, the ``-o`` option can
be used. This also shows that you can give the adapter in colorspace and
how to use a different error rate::

    cutadapt --bwa -e 0.15 -a 330201030313112312 -o output.fastq solid.csfasta solid.qual

This does the same as above, but produces BFAST-compatible output,
strips the \_F3 suffix from read names and adds the prefix "abc:" to
them::

    cutadapt -c -e 0.15 -a 330201030313112312 -x abc: --strip-f3 solid.csfasta solid.qual > output.fastq

Bowtie
------

Quality values of colorspace reads are sometimes negative. Bowtie gets
confused and prints this message:

    Encountered a space parsing the quality string for read xyz

BWA also has a problem with such data. Cutadapt therefore converts
negative quality values to zero in colorspace data. Use the option
``--no-zero-cap`` to turn this off.

To Do / Ideas
=============

-  show average error rate
-  In color space and probably also for Illumina data, gapped alignment
   is not necessary
-  use ``str.format`` instead of ``%``
-  allow to change scores at runtime (using command-line parameters)
-  multi-threading
-  ``--progress``
-  run pylint, pychecker
-  length histogram
-  refactor read\_sequences (use classes)
-  put write\_read into a Fast(a\|q)Writer class?
-  allow .txt input/output
-  check whether input is FASTQ although -f fasta is given
-  close on StopIteration
-  search for adapters in the order in which they are given on the
   command line
-  more tests for the alignment algorithm
