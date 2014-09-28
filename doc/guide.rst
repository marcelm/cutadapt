==========
User guide
==========

Basic usage
===========

If you just want to trim a 3' adapter, the basic command-line for cutadapt is::

    cutadapt -a AACCGGTT -o output.fastq input.fastq

The sequence of the adapter is given with the ``-a`` option. Of course, you
need to replace ``AACCGGTT`` with your actual adapter sequence. Reads are read
from the input file ``input.fastq`` and written to the output file
``output.fastq``.

Cutadapt searches for the adapter in all reads and removes it when it finds it.
All reads that were present in the input file will also be present in the output
file, some of them trimmed, some of them not. Even reads that were trimmed
entirely (because the adapter was found in the very beginning) are output. All
of this can be changed with command-line options, explained further down.

A report is printed after cutadapt has finished processing the reads.


Input and output file formats
-----------------------------

Input files for cutadapt need to be in one the these formats:

* FASTA (file name extensions: ``.fasta``, ``.fa``, ``.fna``, ``.csfasta``, ``.csfa``)
* FASTQ (extensions: ``.fastq``, ``.fq``)
* A pair of a FASTA file and a ``.(cs)qual`` file

The latter format is (or was) used for colorspace data from the SOLiD
instruments.

The input file format is recognized from the file name extension (given in
parentheses in the list above). You can also explicitly specify which format
the input has by using the ``--format`` option.

The output format is the same as the input format, except for the FASTA/QUAL
pairs -- those will always be converted to FASTQ. Also, cutadapt does not check
the output file name: If you input FASTQ data, but use ``-o output.fasta``, then
the output file will actually be in FASTQ format.


Compressed files
----------------

Cutadapt supports compressed input and output files. Whether an input file
needs to be decompressed or an output file needs to be compressed is detected
automatically by inspecting the file name: If it ends in ``.gz``, then gzip
compression is assumed. You can therefore run cutadapt like this and it works
as expected::

    cutadapt -a AACCGGTT -o output.fastq.gz input.fastq.gz

Note that the all of cutadapt's options that expect a file name support this.

If your Python installation includes support for bzip2 compression, then
bzip2-compressed files are also supported and recognized by their
extension ``.bz2``.

(If you are interested in LZMA compression (``.xz`` and ``.lzma``), contact me.)


Standard input and output
-------------------------

If no output file is specified via the ``-o`` option, then the output is sent to
the standard output stream. Instead of the example command line from above, you
can therefore also write::

    cutadapt -a AACCGGTT input.fastq > output.fastq

There is one difference in behavior if you use cutadapt without ``-o``: The
report is sent to the standard error stream instead of standard output. You
can redirect it to a file like this::

    cutadapt -a AACCGGTT input.fastq > output.fastq 2> report.txt

Anywhere cutadapt expects a file name, you can also write a dash (``-``) in
order to specify that standard input or output should be used. For example::

    tail -n 4 input.fastq | cutadapt -a AACCGGTT - > output.fastq

The ``tail -n 4`` prints out only the last four lines of ``input.fastq``, which
are then piped into cutadapt. Thus, cutadapt will work only on the last read in
the input file.

In most cases, you should probably use ``-`` at most once for an input file and
at most once for an output file, in order not to get mixed output.

You cannot combine ``-`` and gzip compression since cutadapt needs to know the 
file name of the output or input file. if you want to have a gzip-compressed 
output file, use ``-o`` with an explicit name.

One last "trick" is to use ``/dev/null`` as an output file name. This special
file discards everything you send into it. If you only want to see the
statistics output, for example, and do not care about the trimmed reads at all,
you could use something like this::

    cutadapt -a AACCGGTT -o /dev/null input.fastq


Trimming reads
==============

Cutadapt supports trimming of four different kinds of adapters:

=================================================== ===========================
Adapter type                                        Command-line option
=================================================== ===========================
:ref:`3' adapter <three-prime-adapters>`            ``-a ADAPTER``
:ref:`5' adapter <five-prime-adapters>`             ``-g ADAPTER``
:ref:`Anchored 5' adapter <anchored-adapters>`      ``-g ^ADAPTER``
:ref:`5' or 3' (both possible) <anywhere-adapters>` ``-b ADAPTER``
=================================================== ===========================

By default, all adapters :ref:`are searched error-tolerantly <error-tolerance>`.
Adapter sequences :ref:`may also contain the "N" wildcard
character <wildcards>`.

In addition, it is possible to :ref:`remove a fixed number of
bases <cut-bases>` from the beginning or end of each read, and to :ref:`remove
low-quality bases (quality trimming) <quality-trimming>` from the 3' end of
each read.


.. _three-prime-adapters:

3' adapters
-----------

A 3' adapter is a piece of DNA ligated to the 3' end of the DNA fragment you
are interested in. The sequencer starts the sequencing process at the 5' end of
the fragment and sequences into the adapter if the read is long enough.
The read that it outputs will then have a part of the adapter in the 
end. Or, if the adapter was short and the read length quite long, then the 
adapter will be somewhere within the read (followed by other bases).

For example, assume your fragment of interest is *MYSEQUENCE* and the adapter is
*ADAPTER*. Depending on the read length, you will get reads that look like this::

    MYSEQUEN
    MYSEQUENCEADAP
    MYSEQUENCEADAPTER
    MYSEQUENCEADAPTERSOMETHINGELSE

Use cutadapt's ``-a ADAPTER`` option to remove this type of adapter. This will
be the result::

    MYSEQUEN
    MYSEQUENCE
    MYSEQUENCE
    MYSEQUENCE

As can be seen, cutadapt correctly deals with partial adapter matches, and also
with any trailing sequences after the adapter. Cutadapt deals with 3' adapters
by removing the adapter itself and any sequence that may follow. If the sequence
starts with an adapter, like this::

    ADAPTERSOMETHING

Then the sequence will be empty after trimming. Note that, by default, empty
reads are not discarded and will appear in the output.


.. _five-prime-adapters:

5' adapters
-----------

.. note::
    Unless your adapter may also occur in a degraded form, you probably
    want to use an anchored 5' adapter, described in the next section.

A 5' adapter is a piece of DNA ligated to the 5' end of the DNA fragment of
interest. The adapter sequence is expected to appear at the start of the read,
but may be partially degraded. The sequence may also appear somewhere within
the read. In all cases, the adapter itself and the sequence preceding it is
removed.

Again, assume your fragment of interest is *MYSEQUENCE* and the adapter is
*ADAPTER*. The reads may look like this::

    ADAPTERMYSEQUENCE
    DAPTERMYSEQUENCE
    TERMYSEQUENCE
    SOMETHINGADAPTERMYSEQUENCE

All the above sequences are trimmed to ``MYSEQUENCE`` when you use `-g ADAPTER`.
As with 3' adapters, the resulting read may have a length of zero when the
sequence ends with the adapter. For example, the read ::

    SOMETHINGADAPTER

will be empty after trimming.


.. _anchored-adapters:

Anchored 5' adapters
--------------------

In many cases, the above behavior is not really what you want for trimming 5'
adapters. You may know, for example, that degradation does not occur and that
the adapter is also not expected to be within the read. Thus, you always expect
the read to look like the first example from above::

    ADAPTERSOMETHING

If you want to trim only this type of adapter, use ``-g ^ADAPTER``. The ``^`` is
supposed to indicate the the adapter is "anchored" at the beginning of the read.
In other words: The adapter is expected to be a prefix of the read. Note that
cases like these are also recognized::

    ADAPTER
    ADAPT
    ADA

The read will simply be empty after trimming.

Be aware that cutadapt still searches for adapters error-tolerantly and, in
particular, allows insertions. So if your maximum error rate is sufficiently
high, even this read will be trimmed::

    BADAPTERSOMETHING

The ``B`` in the beginnig is seen as an insertion. If you also want to prevent
this from happening, use the option ``--no-indels`` to disallow insertions and
deletions entirely.


.. _anywhere-adapters:

5' or 3' adapters
-----------------

The last type of adapter is a combination of the 5' and 3' adapter. You can use
it when your adapter is ligated to the 5' end for some reads and to the 3' end
in other reads. This probably does not happen very often, and this adapter type
was in fact originally implemented because the library preparation in an
experiment did not work as it was supposed to.

For this type of adapter, the sequence is specified with ``-b ADAPTER`` (or use
the longer spelling ``--anywhere ADAPTER``). The adapter may appear in the
beginning (even degraded), within the read, or at the end of the read (even
partially). The decision which part of the read to remove is made as follows: If
there is at least one base before the found adapter, then the adapter is
considered to be a 3' adapter and the adapter itself and everything
following it is removed. Otherwise, the adapter is considered to be a 5'
adapter and it is removed from the read, but the sequence after it it remains.

Here are some examples.

============================== =================== =====================
Read before trimming           Read after trimming Detected adapter type
============================== =================== =====================
``MYSEQUENCEADAPTERSOMETHING`` ``MYSEQUENCE``      3' adapter
``MYSEQUENCEADAPTER``          ``MYSEQUENCE``      3' adapter
``MYSEQUENCEADAP``             ``MYSEQUENCE``      3' adapter
``MADAPTER``                   ``M``               3' adapter
``ADAPTERMYSEQUENCE``          ``MYSEQUENCE``      5' adapter
``PTERMYSEQUENCE``             ``MYSEQUENCE``      5' adapter
``TERMYSEQUENCE``              ``MYSEQUENCE``      5' adapter
============================== =================== =====================

The ``-b`` option currently does not work with color space data.


.. _error-tolerance:

Error tolerance
---------------

All searches for adapter sequences are error tolerant. Allowed errors are
mismatches, insertions and deletions. For example, if you search for the
adapter sequence ``ADAPTER`` and the error tolerance is set appropriately
(as explained below), then also ``ADABTER`` will be found (with 1 mismatch),
as well as ``ADAPTR`` (with 1 deletion), and also ``ADAPPTER`` (with 1
insertion).

The level of error tolerance is adjusted by specifying a *maximum error rate*,
which is 0.1 (=10%) by default. Use the ``-e`` option to set a different value.
To determine the number of allowed errors, the maximum error rate is multiplied
by the length of the match (and then rounded off).

What does that mean?
Assume you have a long adapter ``LONGADAPTER`` and it appears in full somewhere
within the read. The length of the match is 11 characters since the full adapter
has a length of 11, therefore 11·0.1=1.1 errors are allowed with the default
maximum error rate of 0.1. This is rounded off to 1 allowed error. So the
adapter will be found within this read::

    SEQUENCELONGADUPTERSOMETHING

If the match is a bit shorter, however, the result is different::

    SEQUENCELONGADUPT

Only 9 characters of the adapter match: ``LONGADAPT`` matches ``LONGADUPT``
with one substitution. Therefore, only 9·0.1=0.9 errors are allowed. Since this
is rounded off to zero allowed errors, the adapter will not be found.

The number of errors allowed for a given adapter match length is also shown in
the report that cutadapt prints::

    Sequence: 'LONGADAPTER'; Length: 11; Trimmed: 2 times.

    No. of allowed errors:
    0-9 bp: 0; 10-11 bp: 1

This tells us what we now already know: For match lengths of 0-9 bases, zero
errors are allowed and for matches of length 10-11 bases, one error is allowed.

The reason for this behavior is to ensure that short matches are not favored
unfairly. For example, assume the adapter has 40 bases and the maximum error
rate is 0.1, which means that four errors are allowed for full-length matches.
If four errors were allowed even for a short match such as one with 10 bases, this would
mean that the error rate for such a case is 40%, which is clearly not what was
desired.

Insertions and deletions can be disallowed by using the option
``--no-indels``.

See also the :ref:`section on details of the alignment algorithm <algorithm>`.


.. _wildcards:

Wildcards
---------

The wildcard character ``N`` in the adapter sequence is supported. It matches
any nucleotide. This is useful for trimming adapters that have a variable
barcode embedded in them::

    cutadapt -a ACGTAANNNNTTAGC -o output.fastq input.fastq

Wildcard characters in the reads are also supported, but this must be
enabled with ``--match-read-wildcards``.


.. _cut-bases:

Removing a fixed numebr of bases
--------------------------------

By using the ``--cut`` option or its abbreviation ``-u``, it is possible to
unconditionally remove bases from the beginning or end of each read. If
the given length is positive, the bases are removed from the beginning
of each read. If it is negative, the bases are removed from the end.

For example, to remove the first five bases of each read::

    cutadapt -u 5 -o trimmed.fastq reads.fastq

To remove the last seven bases of each read::

    cutadapt -u -7 -o trimmed.fastq reads.fastq

The ``-u``/``--cut`` option can be combined with the other options, but
the desired bases are removed *before* any adapter trimming.


.. _quality-trimming:

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


Trimming paired-end reads
-------------------------

Cutadapt supports trimming of paired-end reads, but currently two passes over
the data are required.

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


Multiple adapters
=================

It is possible to specify more than one adapter sequence by using the options
``-a``, ``-b`` and ``-g`` more than once. Any combination is allowed, such as
five ``-a`` adapters and two ``-g`` adapters. Each read will be searched for
all given adapters, but **only the best matching adapter is removed**. (But it
is possible to :ref:`trim more than one adapter from each
read <more-than-one>`). This is how a command may look like to trim one of two
possible 3' adapters::

    cutadapt -a TGAGACACGCA -a AGGCACACAGGG -o output.fastq input.fastq

The adapter sequences can also be read from a FASTA file. Instead of giving an
explicit adapter sequence, you need to write ``file:`` followed by the name of
the FASTA file::

    cutadapt -a file:adapters.fasta -o output.fastq input.fastq

All of the sequences in the file ``adapters.fasta`` will be used as 3'
adapters. The other adapter options ``-b`` and ``-g`` also support this. Again,
only the best matching adapter is trimmed from each read.

When cutadapt has multiple adapter sequences to work with, either specified
explicitly on the command line or via a FASTA file, it decides in the
following way which adapter should be trimmed:

* All given adapter sequences are matched to the read.
* Adapter matches where the overlap length (see the ``-O`` parameter) is too
  small or where the error rate is too high (``-e``) are removed from further
  consideration.
* Among the remaining matches, the one with the **greatest number of matching
  bases** is chosen.
* If there is a tie, the first adapter wins. The order of adapters is the order
  in which they are given on the command line or in which they are found in the
  FASTA file.

If your adapter sequences are all similar and differ only by a variable barcode
sequence, you should use a single adapter sequence instead that
:ref:`contains wildcard characters <wildcards>`.


.. _named-adapters:

Named adapters
--------------

Cutadapt reports statistics for each adapter separately. To identify the
adapters, they are numbered and the adapter sequence is also printed::

    === Adapter 1 ===

    Sequence: AACCGGTT; Length 8; Trimmed: 5 times.

If you want this to look a bit nicer, you can give each adapter a name in this
way::

    cutadapt -a My_Adapter=AACCGGTT -o output.fastq input.fastq

The actual adapter sequence in this example is ``AACCGGTT`` and the name
assigned to it is ``My_Adapter``. The report will then contain this name in
addition to the other information::

    === Adapter 'My_Adapter' ===

    Sequence: TTAGACATATCTCCGTCG; Length 18; Trimmed: 5 times.

When adapters are read from a FASTA file, the sequence header is used as the
adapter name.

Adapter names are also used in column 8 of :ref:`info files <info-file>`.


.. _more-than-one:

Trimming more than one adapter from each read
---------------------------------------------

By default, at most one adapter sequence is removed from each read, even if
multiple adapter sequences were provided. This can be changed by using the
``--times`` option (or its abbreviated form ``-n``). Cutadapt will then search
for all the given adapter sequences repeatedly, either until no adapter match
was found or until the specified number of rounds was reached.

As an example, assume you have a protocol in which a 5' adapter gets ligated
to your DNA fragment, but it's possible that the adapter is ligated more than
once. So your sequence could look like this::

    ADAPTERADAPTERADAPTERMYSEQUENCE

To be on the safe side, you assume that there are at most 5 copies of the
adapter sequence. This command can be used to trim the reads correctly::

    cutadapt -g ^ADAPTER -n 5 -o output.fastq input.fastq

This feature can also be used to search for *5'/3' linked adapters*. For example,
when the 5' adapter is *FIRST* and the 3' adapter is *SECOND*, then the read
could look like this::

    FIRSTMYSEQUENCESECOND

That is, the sequence of interest is framed by the 5' and the 3' adapter. The
following command can be used to trim such a read::

    cutadapt -g ^FIRST -a SECOND -n 2 ...

Support for linked adapters is currently incomplete. For example, it is not
possible to specify that SECOND should only be trimmed when FIRST also occurs.
`See also this feature
request <https://code.google.com/p/cutadapt/issues/detail?id=34>`_, and
comment on it if you would like to see this implemented.


Illumina TruSeq
===============

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


Cutadapt's output
=================

How to read the report
----------------------

After every run, cutadapt prints out per-adapter statistics. The output
starts with something like this::

    Sequence: 'ACGTACGTACGTTAGCTAGC'; Length: 20; Trimmed: 2402 times.

The meaning of this should be obvious.

The next piece of information is this::

    No. of allowed errors:
    0-9 bp: 0; 10-19 bp: 1; 20 bp: 2

The adapter has, as was shown above, has a length of 20
characters. We are using the default error rate of 0.1. What this
implies is shown above: Matches up to a length of 9 bp are allowed to
have no errors. Matches of lengths 10-19 bp are allowd to have 1 error
and matches of length 20 can have 2 errors. See also :ref:`the section about
error-tolerant matching <error-tolerance>`.

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
since we know already that no errors are allowed at lengths 0-9 bp).

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


.. _info-file:

Format of the info file
-----------------------

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

The concatenation of the fields 5-7 yields the full read sequence. Column 8
identifies the found adapter. `The section about named
adapters <named-adapters>` describes how to give a name to an adapter. Adapters
without a name are numbered starting from 1.

If no adapter was found, the format is as follows:

1.  Read name
2.  The value -1
3.  The read sequence

When parsing that file, be aware that additional columns may be added in
the future. Note also that some fields can be empty, resulting in
consecutive tabs within a line. Also, in the current version, when the
``--times`` option is set to a value other than 1 (the default value),
multiple lines are written to the info file for each read.


.. _algorithm:

The alignment algorithm
=======================

Since the publication of the `EMBnet journal application note about
cutadapt <http://dx.doi.org/10.14806/ej.17.1.200>`_, the alignment algorithm
used for finding adapters has changed significantly. This new algorithm is
described in this section.
An even more detailed description of the algorithm is available in Chapter 2
of my PhD thesis `Algorithms and tools for the analysis of high-throughput DNA
sequencing data <http://hdl.handle.net/2003/31824>`_.

The algorithm is based on *semiglobal alignment*, also called *free-shift*,
*ends-free* or *overlap* alignment.


.. note::
    The documentation is currently being worked on. Text until here has
    been re-written. Text below may be not in the correct order or incomplete.


The alignment uses unit costs, which means that mismatches, insertions and deletions are
counted as one error.


An optimal alignment fulfills all of these criteria:

- its error_rate is at most max_error_rate
- Among those alignments with error_rate <= max_error_rate, the alignment contains
  a maximal number of matches (there is no alignment with more matches).
- If there are multiple alignments with the same no. of matches, then one that
  has minimal no. of errors is chosen.
- If there are still multiple candidates, choose the alignment that starts at the
  leftmost position within the read.


The new alignment algorithm checks the error rate while aligning and only
reports alignments that do not have too many errors.

maximizing matches -- example: read TCGTATGCCCTCC and adapter TCGTATGCCGTCTTC, max. error rate 20%.

Here is the alignment that one would expect:

TCGTATGCCCTCC   (read)
=========X==X
TCGTATGCCGTCTTC (adapter)

But it turns out that the actual alignment that is found is this one:

TCGTATGCCGTCTTC
=========X==XX=
TCGTATGCCCTC--C

Since it has length 15 and contains three errors, the error rate is 3/15=20%.

To understand why that alignment was chosen, one needs to know that cutadapt
doesn't really care about the actual number of errors as long as the alignment
has not more errors than the maximum error rate allows (20% in your case). The
alignment algorithm doesn't try to find an alignment with few errors, but
instead it tries to find one with many matches. And in this case, the first
alignment has 11 matches, and the second one has 12.

This is sometimes a bit surprising, but consistent with the way the algorithm
was designed to behave.

As an explanation, the problem lies in the maximum error rate:
Users should be able to specify that as a parameter because it is quite easy to
understand.
The problem then is that one cannot optimize for 'as few errors as possible'
since an alignment in which the read and the adapter do not overlap at all is
always optimal (since it has zero errors). Typically, this is solved by using
alignment scores, which are not as intuitive.
Instead, the number of matches is optimized, but still under the constraint that
the actual error rate may not go above the specified maximum error rate.

One other important point to note is that this particularity doesn't influence
the trimming result: Try to run cutadapt with the maximum error rate reduced to
0.16. This prevents the second alignment with 3 errors to be found and suddenly
the reported number of errors in the info file is down to 2.
