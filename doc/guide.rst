==========
User guide
==========

Basic usage
===========

To trim a 3' adapter, the basic command-line for Cutadapt is::

    cutadapt -a AACCGGTT -o output.fastq input.fastq

The sequence of the adapter is given with the ``-a`` option. You need to replace
``AACCGGTT`` with the correct adapter sequence. Reads are read from the input
file ``input.fastq`` and are written to the output file ``output.fastq``.

Compressed in- and output files are also supported::

    cutadapt -a AACCGGTT -o output.fastq.gz input.fastq.gz

Cutadapt searches for the adapter in all reads and removes it when it finds it.
Unless you use a filtering option, all reads that were present in the input file
will also be present in the output file, some of them trimmed, some of them not.
Even reads that were trimmed to a length of zero are output. All of this can be
changed with command-line options, explained further down.

:ref:`Trimming of paired-end data <paired-end>` is also supported.


Input and output file formats
-----------------------------

The supported input and output file formats are FASTA and FASTQ, with
optional compression.

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

The supported compression formats are gzip (``.gz``), bzip2 (``.bz2``)
and xz (``.xz``).

The default compression level for gzip output is 6. Use option ``-Z`` to
change this to level 1. The files need more space, but it is faster and
therefore a good choice for short-lived intermediate files.

If available, Cutadapt uses `pigz <https://zlib.net/pigz/>`_ to speed up
writing and reading of gzipped files.


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

    cutadapt -a AACCGGTT -o /dev/null input.fastq


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

Make also sure that you have ``pigz`` (parallel gzip) installed if you use
multiple cores and write to a ``.gz`` output file. Otherwise, compression of
the output will be done in a single thread and therefore be a bottleneck.


.. versionadded:: 1.15

.. versionadded:: 1.18
    ``--cores=0`` for autodetection

.. versionadded:: 2.5
    Multicore works with ``--untrimmed/too-short/too-long-(paired)-output``

.. versionadded:: 2.7
    Multicore works with ``--info-file``, ``--rest-file``, ``--wildcard-file``

.. versionadded:: 3.0
    Multicore support for demultiplexing added.


Speed-up tricks
---------------

There are several tricks for limiting wall-clock time while using Cutadapt.

``-Z`` (equivalent to ``--compression-level=1``) can be used to limit the
amount of CPU time which is spent on the compression of output files.
Alternatively, choosing filenames not ending with ``.gz``, ``.bz2`` or ``.xz``
will make sure no CPU time is spent on compression at all.  On systems
with slow I/O, it can actually be faster to set a higher compression-level
than 1.

Increasing the number of cores with ``-j`` will increase the number of reads per
minute at near-linear rate.

It is also possible to use pipes in order to bypass the filesystem and pipe
cutadapt's output into an aligner such as BWA. The ``mkfifo`` command allows
you to create named pipes in bash.

.. code-block::bash

    mkfifo R1.fastq R2.fastq
    cutadapt -a ${ADAPTER_R1} -A ${ADAPTER_R2} -o R1.fastq -p R2.fastq ${READ1} ${READ2} > cutadapt.report & \
    bwa mem ${INDEX} R1.fastq R2.fastq

This command will run cutadapt and BWA simultaneously, using Cutadapt’s output as
BWA’s input, and capturing Cutadapt’s report in ``cutadapt.report``.

Read processing stages
======================

Cutadapt can do a lot more in addition to removing adapters. There are various
command-line options that make it possible to modify and filter reads and to
redirect them to various output files. Each read is processed in the following
order:

1. :ref:`Read modification options <modifying-reads>` are applied. This includes
   :ref:`adapter removal <adapter-types>`,
   :ref:`quality trimming <quality-trimming>`, read name modifications etc. The
   order in which they are applied is the order in which they are listed in the
   help shown by ``cutadapt --help`` under the “Additional read modifications”
   heading. Adapter trimming itself does not appear in that list and is
   done after quality trimming and before length trimming (``--length``/``-l``).

2. :ref:`Filtering options <filtering>` are applied, such as removal of too
   short or untrimmed reads. Some of the filters also allow to redirect a read
   to a separate output file.  The filters are applied in the order in which
   they are listed in the help shown by ``cutadapt --help`` under the
   “Filtering of processed reads” heading.
3. If the read has passed all the filters, it is written to the output file.


.. _adapter-types:

Adapter types
=============

Cutadapt can detect multiple adapter types. 5' adapters preceed the sequence of
interest and 3' adapters follow it. Further distinctions are made according to
where in the read the adapter sequence is allowed to occur.

========================================================= =============================
Adapter type                                              Command-line option
========================================================= =============================
:ref:`Regular 3' adapter <three-prime-adapters>`          ``-a ADAPTER``
:ref:`Regular 5' adapter <five-prime-adapters>`           ``-g ADAPTER``
:ref:`Non-internal 3' adapter <non-internal>`             ``-a ADAPTERX``
:ref:`Non-internal 5' adapter <non-internal>`             ``-g XADAPTER``
:ref:`Anchored 3' adapter <anchored-3adapters>`           ``-a ADAPTER$``
:ref:`Anchored 5' adapter <anchored-5adapters>`           ``-g ^ADAPTER``
:ref:`5' or 3' (both possible) <anywhere-adapters>`       ``-b ADAPTER``
:ref:`Linked adapter <linked-adapters>`                   | ``-a ^ADAPTER1...ADAPTER2``
                                                          | ``-g ADAPTER1...ADAPTER2``
========================================================= =============================

By default, all adapters :ref:`are searched error-tolerantly <error-tolerance>`.
Adapter sequences :ref:`may also contain any IUPAC wildcard
character (degenerate bases) <wildcards>` (such as ``N``).

In addition, it is possible to :ref:`remove a fixed number of
bases <cut-bases>` from the beginning or end of each read, to :ref:`remove
low-quality bases (quality trimming) <quality-trimming>` from the 3' and 5' ends,
and to :ref:`search for adapters also in the reverse-complemented reads <reverse-complement>`.


Overview of adapter types
-------------------------

3' adapter types
~~~~~~~~~~~~~~~~

A 3' adapter is assumed to be ligated to the 3' end of your sequence of interest.
When such an adapter is found, the adapter sequence itself and the sequence
following it (if there is any) are trimmed. This table shows in which ways
the different 3' adapter types are allowed to occur in a read in order to be
recognized by the program.

================================== =================== ======================== ============================= =========================
Adapter location in read           Read layout         | Found by regular 3’    | Found by non-internal 3’    | Found by anchored 3’
                                                       | ``-a ADAPTER``         | ``-a ADAPTERX``             | ``-a ADAPTER$``
================================== =================== ======================== ============================= =========================
Full adapter sequence anywhere     acgtacgtADAPTERacgt                      yes                           no                         no
Partial adapter sequence at 3’ end acgtacgtacgtADAP                         yes                           yes                        no
Full adapter sequence at 3’ end    acgtacgtacgtADAPTER                      yes                           yes                       yes
================================== =================== ======================== ============================= =========================


5' adapter types
~~~~~~~~~~~~~~~~

A 5' adapter is assumed to be ligated to the 5' end of your sequence of interest.
When such an adapter is found, the adapter sequence itself and the sequence
preceding it (if there is any) are trimmed. This table shows in which ways
the different 5' adapter types are allowed to occur in a read in order to be
recognized by the program.

================================== =================== ======================== ============================= =========================
Adapter location in read           Read layout         | Found by regular 5’    | Found by non-internal 5’    | Found by anchored 5’
                                                       | ``-g ADAPTER``         | ``-g XADAPTER``             | ``-g ^ADAPTER``
================================== =================== ======================== ============================= =========================
Full adapter sequence anywhere     acgtADAPTERacgtacgt                      yes                           no                         no
Partial adapter sequence at 5’ end PTERacgtacgtacgt                         yes                           yes                        no
Full adapter sequence at 5’ end    ADAPTERacgtacgtacgt                      yes                           yes                       yes
================================== =================== ======================== ============================= =========================


.. _three-prime-adapters:

Regular 3' adapters
-------------------

A 3' adapter is a piece of DNA ligated to the 3' end of the DNA fragment of
interest. The sequencer starts the sequencing process at the 5' end of the
fragment. If the fragment is shorter than the read length, the sequencer
will sequence into the adapter and the reads will thus contain some part
of the adapter. Depending on how much longer the read is than the fragment
of interest, the adapter occurs 1) not at all, 2) partially or fully at the
end of the read (not followed by any other bases), or 3) in full somewhere
within the read, followed by some other bases.

Use Cutadapt’s ``-a`` option to find and trim such an adapter, allowing
both partial and full occurrences.

For example, assume your fragment of interest is *mysequence* and the adapter is
*ADAPTER*. Depending on the read length, you will get reads that look like this::

    mysequen
    mysequenceADAP
    mysequenceADAPTER
    mysequenceADAPTERsomethingelse

Using ``-a ADAPTER`` to remove this type of adapter, this will
be the result::

    mysequen
    mysequence
    mysequence
    mysequence

As this example shows, Cutadapt allows regular 3' adapters to occur in full
anywhere within the read (preceeded and/or succeeded by zero or more bases), and
also partially degraded at the 3' end. Cutadapt deals with 3' adapters
by removing the adapter itself and any sequence that may follow. As a consequence,
a sequence that starts with an adapter, like this, will be trimmed to an empty read::

    ADAPTERsomething

By default, empty reads are kept and will appear in the output. If you do not
want this, use the ``--minimum-length``/``-m`` :ref:`filtering option <filtering>`.


.. _five-prime-adapters:

Regular 5' adapters
-------------------

.. note::
    Unless your adapter may also occur in a degraded form, you probably
    want to use an :ref:`anchored 5' adapter <anchored-3adapters>`.

A 5' adapter is a piece of DNA ligated to the 5' end of the DNA fragment of
interest. For this type of adapter to be found, the adapter sequence needs to
either appear in full somewhere within the read (internal match) or at the
start (5' end) of it, where in the latter case also partial occurrences are
allowed. In all cases, the adapter itself and the sequence preceding it is
removed.

Assume your fragment of interest is *mysequence* and the adapter is
*ADAPTER*. The reads may look like this::

    ADAPTERmysequence
    DAPTERmysequence
    TERmysequence
    somethingADAPTERmysequence

All the above sequences are trimmed to ``mysequence`` when you use `-g ADAPTER`.
As with 3' adapters, the resulting read may have a length of zero when the
sequence ends with the adapter. For example, the read ::

    somethingADAPTER

will be empty after trimming.


.. _anchored-5adapters:

Anchored 5' adapters
--------------------

An anchored 5' adapter is an adapter that is expected to occur in full
length at the beginning of the read. Example::

    ADAPTERsomething

This is usually how forward PCR primers are found in the read in amplicon
sequencing, for instance. In Cutadapt’s terminology, this type of adapter
is called "anchored" to distinguish it from :ref:`"regular" 5'
adapters <anchored-3adapters>`, which are 5' adapters with a less strict
placement requirement.

If the adapter sequence is ``ADAPTER``, use ``-g ^ADAPTER`` to remove an
anchored 5' adapter. The ``^`` is meant to indicate the "anchoring" to the
beginning of the read. With this, the example read ``ADAPTERsomething`` is
trimmed to just ``something``.

An anchored 5' adapter must occur in full at the beginning of the read.
If the read happens to be shorter than the adapter, partial occurrences
such as ``ADAPT`` are not found.

The requirement for a full match at the beginning of the read is relaxed
when Cutadapt searches error-tolerantly, as it does by default. In
particular, insertions and deletions may allow reads such as these to be
trimmed, assuming the maximum error rate is sufficiently high::

    BADAPTERsomething
    ADAPTE

The ``B`` in the beginning is seen as an insertion, and the missing ``R``
as a deletion. If you also want to prevent this from happening, use the
option ``--no-indels``, which disallows insertions and deletions entirely.



.. _anchored-3adapters:

Anchored 3' adapters
--------------------

It is also possible to anchor 3' adapters to the end of the read. This is
useful, for example, if you work with merged overlapping paired-end
reads. Add the ``$`` character to the end of an
adapter sequence specified via ``-a`` in order to anchor the adapter to the
end of the read, such as ``-a ADAPTER$``. The adapter will only be found if it
occurs in full at the end of the read (that is, it must be a *suffix* of the
read.

The requirement for a full match exactly at the end of the read is relaxed when
Cutadapt searches error-tolerantly, as it does by default.
You can disable insertions and deletions with ``--no-indels``.

Anchored 3' adapters work as if you had reversed the sequence and used an
appropriate anchored 5' adapter.

As an example, assume you have these reads::

    mysequenceADAP
    mysequenceADAPTER
    mysequenceADAPTERsomethingelse

Using ``-a ADAPTER$`` will result in::

    mysequenceADAP
    mysequence
    mysequenceADAPTERsomethingelse

That is, only the middle read is trimmed at all.


.. _non-internal:

Non-internal 5' and 3' adapters
-------------------------------

The non-internal 5' and 3' adapter types disallow internal occurrences of the
adapter sequence. This is like a less strict version of anchoring: The
adapter must always be at one of the ends of the read, but - unlike anchored
adapters - partial occurrences are also ok.

Use ``-a ADAPTERX`` (replace ``ADAPTER`` with your actual adapter sequence, but
use a literal ``X``) to disallow internal matches for a 3' adapter. Use
``-g XADAPTER`` to disallow them for a 5' adapter.
Mnemonic: The ``X`` is not allowed to “shift into” the read.

Here are some examples for trimming reads with ``-a ADAPTERX``:

================================== ==================================
Input read                         Processed read
================================== ==================================
``mysequenceADAP``                 ``mysequence``
``mysequenceADAPTER``              ``mysequence``
``mysequenceADAPTERsomethingelse`` ``mysequenceADAPTERsomethingelse``
================================== ==================================

Here are some examples for trimming reads with ``-g XADAPTER``:

================================== ===================================
Input read                         Processed read
================================== ===================================
``APTERmysequence``                ``mysequence``
``ADAPTERmysequence``              ``mysequence``
``somethingelseADAPTERmysequence`` ``somethingelseADAPTERmysequence``
================================== ===================================

.. versionadded:: 1.17

.. _linked-adapters:

Linked adapters (combined 5' and 3' adapter)
--------------------------------------------

If your sequence of interest is surrounded by a 5' and a 3' adapter, and you want
to remove both adapters, then you can use a *linked adapter*. A linked
adapter combines a 5' and a 3' adapter. By default, the adapters are not anchored,
but in many cases, you should anchor the 5’ adapter by prefixing it with ``^``.

:ref:`See the previous sections <anchored-5adapters>` for what anchoring means.

.. note::
   Cutadapt versions before 2.0 anchored the 5’ adapter within linked adapters
   automatically even if the initial ``^`` was not specified. If you have scripts
   written for Cutadapt versions earlier than 2.0, please add the ``^`` so that
   the behavior does not change!

Linked adapters are specified as two sequences separated by ``...`` (three dots)::

    cutadapt -a ^ADAPTER1...ADAPTER2 -o out.fastq.gz in.fastq.gz

If you anchor an adapter, it will also become marked as being *required*. If a
required adapter cannot be found, the read will not be trimmed at all even if
the other adapter occurs. If an adapter is not required, it is *optional*.

Also, when you use the ``--discard-untrimmed`` option (or ``--trimmed-only``) with a
linked adapter, then a read is considered to be trimmed only if all required adapters
were found.

In the previous example, ``ADAPTER1`` was anchored and therefore required, but ``ADAPTER2``
was optional. Anchoring also ``ADAPTER2`` (and making it required as well) would look like this::

    cutadapt -a ^ADAPTER1...ADAPTER2$ -o out.fastq.gz in.fastq.gz

As an example, assume the 5' adapter is *FIRST*, the 3' adapter is *SECOND*
and you have these input reads::

    FIRSTmysequenceSECONDextrabases
    FIRSTmysequenceSEC
    FIRSTmyseque
    anotherreadSECOND

Trimming with ::

    cutadapt -a ^FIRST...SECOND -o output.fastq input.fastq

will result in ::

    mysequence
    mysequence
    myseque
    anotherreadSECOND

The 3' adapter in the last read is not trimmed because the anchored 5’ adapter is required, but
missing in the read.

Linked adapters do not work when used in combination with ``--info-file`` and ``--mask-adapter``.

To provide :ref:`adapter-search parameters <search-parameters>`
for linked adapters, they need to be set for each constituent adapter separately, as in
``-g "ADAPTER1;min_overlap=5...ADAPTER2;min_overlap=6"``.

.. versionadded:: 1.10

.. versionadded:: 1.13
   Ability to anchor the 3' adapter.

.. versionadded:: 2.0
   The 5’ adapter is no longer anchored by default.


.. _linked-override:

Changing which adapters are required
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As described, when you specify a linked adapter with ``-a``, the adapters that are anchored
become *required*, and the non-anchored adapters become *optional*. To change this, you can
instead use ``-g`` to specify a linked adapter. In that case, *both* adapters are required
(even if they are not anchored). This type of linked adapter type is especially suited for
trimming CRISPR screening reads. For example::

    cutadapt -g ADAPTER1...ADAPTER2 -o out.fastq.gz in.fastq.gz

Here, both ``ADAPTER1`` and ``ADAPTER2`` are not anchored, but they are required because ``-g``
was used.

The ``-g`` option does not cover all cases, so you can also mark each adapter explicitly as
required or optional using the :ref:`search parameters <search-parameters>`
``required`` and ``optional``. This is the only way to make an anchored adapter optional.
For example, to request that an anchored 5' adapter (here ``ADAPTER1``) should not be required,
you can specify it like this ::

    cutadapt -a "^ADAPTER1;optional...ADAPTER2" -o output.fastq.gz input.fastq.gz

.. versionadded:: 1.13
    Option ``-g`` added.

.. versionchanged:: 1.15
    Option ``-g`` requires both adapters.


Linked adapter statistics
~~~~~~~~~~~~~~~~~~~~~~~~~

For linked adapters, the statistics report contains a line like this::

    === Adapter 1 ===

    Sequence: AAAAAAAAA...TTTTTTTTTT; Type: linked; Length: 9+10; Trimmed: 3 times; Half matches: 2

The value for “Half matches” tells you how often only the 5'-side of the adapter was found, but not
the 3'-side of it. This applies only to linked adapters with regular (non-anchored) 3' adapters.


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
adapter and it is removed from the read, but the sequence after it remains.

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

.. _rightmost:

Multiple adapter occurrences within a single read
-------------------------------------------------

If a single read contains multiple copies of the same adapter, the basic rule is
that the leftmost match is used for both 5' and 3' adapters. For example, when
searching for a 3' adapter in ::

    cccccADAPTERgggggADAPTERttttt

the read will be trimmed to ::

    ccccc

When the adapter is a 5' adapter instead, the read will be trimmed to ::

    gggggADAPTERttttt

For 5' adapters, this can be changed so that the *rightmost* occurrence is found
by using the ``rightmost`` :ref:`search parameter <search-parameters>`, as in
``-g "ACGT;rightmost"``.


.. versionadded:: 4.1
   The ``rightmost`` search parameter


.. _trimming-parameters:

.. _search-parameters:

Adapter-search parameters
=========================

The adapter search algorithm has a few parameters specific to each adapter
that control how the adapter sequence is found. The command-line options ``-e``
and ``-O`` set the maximum error rate and minimum overlap parameters (see
details in the following sections) for all
adapters listed via the ``-a``/``-b``/``-g`` etc. options. When trimming more
than one adapter, it may be necessary to change search parameters for each
adapter individually. You can do so by adding a semicolon and ``parameter=value`` to the end
of the adapter sequence, as in ``-a "ADAPTER;max_error_rate=0.2"``. There are also "flags"
that enable certain behavior. These are written without the ``=value`` part.

Multiple parameters can be set, as in ``-a "ADAPTER;max_error_rate=0.2;min_overlap=5"``.
For linked adapters, search parameters need to be specified separately for each adapter
as in ``-g "ADAPTER1;min_overlap=5...ADAPTER2;min_overlap=6"``.

Remember to add the quotation marks; otherwise the shell will interpret the semicolon as a
separator between two commands.

The following parameters are supported:

======================================================= =============== ================================
Parameter                                               Global option   Adapter-specific parameter
======================================================= =============== ================================
Maximum error rate (default: 0.1)                       ``-e 0.2``      | ``ADAPTER;e=0.2`` or
                                                                        | ``ADAPTER;max_errors=0.2`` or
                                                                        | ``ADAPTER;max_error_rate=0.2``

Minimum overlap (default: 3)                            ``-O 5``        | ``ADAPTER;o=5`` or
                                                                        | ``ADAPTER;min_overlap=5``

Disallow indels                                         ``--no-indels`` ``ADAPTER;noindels``
Allow indels (this is the default)                                      ``ADAPTER;indels``
Allow matches anywhere                                                  ``ADAPTER;anywhere``

:ref:`Linked adapter required <linked-override>`                        ``ADAPTER;required``
:ref:`Linked adapter optional <linked-override>`                        ``ADAPTER;optional``
:ref:`Find rightmost 5' adapter occurrence <rightmost>`                 ``ADAPTER;rightmost``
======================================================= =============== ================================

The minimum overlap length cannot be set for anchored adapters as these always need to occur at full
length.

When using the ``file:`` notation to read in adapters from a FASTA file,
it is possible to specify file-specific search parameters::

    cutadapt -a "file:adapters.fa;min_overlap=5;noindels"

The individual adapter specifications in the FASTA file can also contain search parameters::

    >adapter1
    ^ACGT;min_overlap=3
    >adapter2
    AACCGGT;noindels

More specific parameters override less specific ones:

1. Adapter-specific parameters override the file-specific settings
2. File-specific search parameters override the global settings

.. versionadded:: 1.18
    Syntax for setting adapter-specific search parameters

.. versionadded:: 3.5
    The ``indels`` and ``noindels`` parameters.

.. versionadded:: 4.1
    Support file-specific search parameters (when using the ``file:`` notation)

.. versionadded: 4.1
    The ``rightmost`` search parameter

.. _error-tolerance:

Error tolerance
---------------

All searches for adapter sequences are error tolerant. Allowed errors are
mismatches, insertions and deletions. For example, if you search for the
adapter sequence ``ADAPTER`` and the error tolerance is set appropriately
(as explained below), then also ``ADABTER`` will be found (with 1 mismatch),
as well as ``ADAPTR`` (with 1 deletion), and also ``ADAPPTER`` (with 1
insertion). If insertions and deletions are disabled with ``--no-indels``,
then mismatches are the only type of errors.

The level of error tolerance is determined by a *maximum error rate*, which is
0.1 (=10%) by default. An adapter occurrence is only found if the actual
error rate of the match does not exceed the maximum error rate. The actual
error rate is computed as the *number of errors in the match*
divided by the *length of the matching part of the adapter*.

For example, an adapter match of length 8 containing 1 error has an error rate
of 1/8=0.125. At the default maximum error rate 0.1, it would not be found, but
a match of length 10 containing 1 error has an error rate of 1/10=0.1 and would
be found.

Relating the number of errros to the length of the matching part of the
adapter is important because Cutadapt allows for partial adapter
occurrences (for the non-anchored adapter types). If only the absolute
number of errors were used, shorter matches would be favored unfairly. For
example, assume an adapter has 30 bases and we allow three errors over that
length. If we allowed these three errors even for a partial occurrences of,
for example, four bases, we can immediately see that this results in
unexpected matches. Using the error rate as a criterion helps to keep
sensitivity and specificity roughly the same over the possible lengths of
the matches.

The ``-e`` option on the command line allows you to change the maximum error rate.
If the value is between 0 and 1 (but not 1 exactly), then this sets the maximum
error rate directly for all specified adapters. The default is ``-e 0.1``. You
can also use the adapter-specific parameter ``max_error_rate`` or ``max_errors``
or just ``e`` to override the default for a single adapter only.
Examples: ``-a "ADAPTER;max_error_rate=0.15"``, ``-a "ADAPTER;e=0.15"``
(the quotation marks are necessary).

Alternatively, you can also specify a value of 1 or greater as the number of
allowed errors, which is then converted to a maximum error rate for each adapter
individually. For example, with an adapter of length 10, using ``-e 2`` will
set the maximum error rate to 0.2 for an adapter of length 10.

The value does not have to be an integer, and if you use an adapter type
that allows partial matches, you may want to add 0.5 to the desired number of
errors, which achieves that even slightly shorter than full-lengths
matches will be allowed at the specified number of errors. In short, if you
want to allow two errors, use ``-e 2.5``.

This also works in the adapter-specific parameters.
Examples: ``-a "ADAPTER;e=1"``, ``-a "ADAPTER;max_errors=2.5"``. Note that
``e``, ``max_error_rate`` and ``max_errors`` are all equivalent and the
decision whether a rate or an absolute number is meant is based on
whether the given value is less than 1 or not.

The number of errors allowed for a given adapter match length is also shown under
the “No. of allowed errors” heading in the report that Cutadapt prints::

    Sequence: 'SOMEADAPTER'; Length: 11; Trimmed: 2 times.

    No. of allowed errors:
    0-9 bp: 0; 10-11 bp: 1

This tells us: For match lengths of 0-9 bases, zero errors are allowed and for
matches of length 10-11 bases, one error is allowed.

See also the :ref:`section on details of the alignment algorithm <adapter-alignment-algorithm>`.

.. versionadded:: 2.11
    Allow specifying the number of errors

N wildcard characters
~~~~~~~~~~~~~~~~~~~~~

Any ``N`` wildcard characters in the adapter sequence are skipped when
computing the error rate. That is, they do not contribute to the length of
a match. For example, the adapter sequence ``ACGTACNNNNNNNNGTACGT`` has a length
of 20, but only 12 non-``N``-characters. At a maximum error rate of 0.1, only
one error is allowed if this sequence is found in full in a read because
12·0.1=1.2, which is 1 when rounded down.

This is done because ``N`` bases cannot contribute to the number of errors.
In previous versions, ``N`` wildcard characters did contribute to the match
length, but this artificially inflates the number of allowed errors. For example,
an adapter like ``N{18}CC`` (18 ``N`` wildcards followed by ``CC``) would
effectively match anywhere because the default error rate of 0.1 would allow for
two errors, but there are only two non-``N`` bases in the particular adapter.

However, even in previous versions, the location with the greatest number of
matching bases is chosen as the best location for an adapter, so in many cases
the adapter would still be placed properly.

.. versionadded:: 2.0
    Ignore ``N`` wildcards when computing the error rate.


.. _minimum-overlap:
.. _random-matches:

Minimum overlap (reducing random matches)
-----------------------------------------

Since Cutadapt allows partial matches between the read and the adapter sequence
for most adapter types, short matches can occur by chance, leading to erroneously
trimmed bases. For
example, just by chance, we expect that roughly 25% of all reads end with a base
that is identical to the first base of the adapter. To reduce the number of
falsely trimmed bases,
the alignment algorithm requires that at least *three bases* of the adapter
are aligned to the read.

This minimum overlap length can be changed globally (for all adapters) with the parameter
``--overlap`` (or its short version ``-O``). The option is ignored for
anchored adapters since these do not allow partial matches.

Alternatively, use the adapter-specific
parameter ``min_overlap`` to change it for a single adapter only. Example:
``-a "ADAPTER;min_overlap=5"`` (the quotation marks are necessary).
For anchored adapters, attempting to set a minimum overlap this way will
result in an error.

In :ref:`linked adapters <linked-adapters>`, the minimum overlap length is applied
separately to the 5' and the 3' adapter.

If a read contains a partial adapter sequence shorter than the minimum overlap length,
no match will be found (and therefore no bases are trimmed).

Requiring at least three bases to match is quite conservative. Even if no
minimum overlap was required, we can compute that we lose only about 0.44 bases
per read on average, see `Section 2.3.3 in my
thesis <http://hdl.handle.net/2003/31824>`_. With the default minimum
overlap length of 3, only about 0.07 bases are lost per read.

When choosing an appropriate minimum overlap length, take into account that
true adapter matches are also lost when the overlap length is higher than
zero, reducing Cutadapt's sensitivity.

It is possible that fewer bases are removed from a read than the minimum
overlap length seems to imply. The overlap length is the number of bases
in the adapter that got aligned to the read, which means that if there are
deletions in the adapter, the corresponding part in the read will be shorter.
(This is only relevant when the maximum allowed error rate and/or the minimum
overlap length are changed such that at least one error is allowed over the
given length.)


Allowing partial matches at both ends
-------------------------------------

The regular 5' and 3' adapter types allow partial adapter occurrences only
at the 5' and 3' end of the read, respectively. To allow partial matches at both ends,
you can use the ``anywhere`` adapter-specific parameter.

A 3' adapter specified via ``-a ADAPTER`` will be found even
when it occurs partially at the 3' end, as in ``mysequenceADAPT``. However,
it will by default not be found if it occurs partially at the 5' end, as in
``APTERmysequence``. To find the adapter in both cases, specify
the adapter as ``-a "ADAPTER;anywhere"``.

Similarly, for a 5' adapter specified via ``-g ADAPTER``, partial matches at
the 3' end are not found, as in ``mysequenceADAPT``. To allow partial matches
at both ends, use ``-g "ADAPTER;anywhere"``.

.. note::
    With ``anywhere``, partial matches at the end that is usually not allowed
    to be matched will result in empty reads! This means that short random
    matches have a much greater detrimental effect and you should
    :ref:`increase the minimum overlap length <random-matches>`.


.. _reverse-complement:

Searching reverse complements
-----------------------------

By default, Cutadapt expects adapters to be given in the same orientation (5' to 3') as the reads.
That is, neither reads nor adapters are reverse-complemented.

To change this, use option ``--revcomp`` or its abbreviation ``--rc``. If given, Cutadapt searches
both the read and its reverse complement for adapters. If the reverse complemented read yields
a better match, then that version of the read is kept. That is, the output file will contain the
reverse-complemented sequence. This can be used to “normalize” read orientation/strandedness.

To determine which version of the read yields the better match, the full adapter search (possibly
multiple rounds if ``--times`` is used) is done independently on both versions, and the version that
results in the higher number of matching nucleotides is considered to be the better one.

The name of a reverse-complemented read is changed by adding a space and ``rc`` to it. (Please
file an issue if you would like this to be configurable.)

The report will show the number of reads that were reverse-complemented, like this::

    Total reads processed:  60
    Reads with adapters:    50 (83.3%)
    Reverse-complemented:   20 (33.3%)

Here, 20 reverse-complemented reads contain an adapter and 50 - 20 = 30 reads that did not need to
be reverse-complemented contain an adapter.

Option ``--revcomp`` is currently available only for single-end data.

.. versionadded:: 2.8


Specifying adapter sequences
============================

.. _wildcards:

Wildcards
---------

All `IUPAC nucleotide codes <http://www.bioinformatics.org/sms/iupac.html>`_
(wildcard characters, degenerate bases) are supported.
For example, use an ``N`` in the adapter
sequence to match any nucleotide in the read, or use ``-a YACGT`` for an adapter
that matches both ``CACGT`` and ``TACGT``. The wildcard character ``N`` is
useful for trimming adapters with an embedded variable barcode::

    cutadapt -a ACGTAANNNNTTAGC -o output.fastq input.fastq

Even the ``X`` wildcard that does not match any nucleotide is supported. If
used as in ``-a ADAPTERX`` or ``-g XADAPTER``, it acquires a special meaning for
the matching algorithm
:ref:`and disallows internal adapter matches <non-internal>`.

The character ``I``, used to encode the base inosine, is automatically
replaced with ``N`` within the adapter sequence.

Wildcard characters are by default only allowed in adapter sequences and
are not recognized when they occur in a read. This is to avoid matches in reads
that consist of many (often low-quality) ``N`` bases. Use
``--match-read-wildcards`` to enable wildcards also in reads.

Use the option ``-N`` to disable interpretation of wildcard characters even in
the adapters. If wildcards are disabled entirely, that is, when you use ``-N``
and *do not* use ``--match-read-wildcards``, then Cutadapt compares characters
by their ASCII value. Thus, both the read and adapter can be arbitrary strings
(such as ``SEQUENCE`` or ``ADAPTER`` as used here in the examples).


.. versionadded:: 4.2
   Inosine ``I``

Repeated bases
--------------

If you have many repeated bases in the adapter sequence, such as many ``N`` s or
many ``A`` s, you do not have to spell them out. For example, instead of writing
ten ``A`` in a row (``AAAAAAAAAA``), write ``A{10}`` instead. The number within
the curly braces specifies how often the character that preceeds it will be
repeated. This works also for IUPAC wildcard characters, as in ``N{5}``.

It is recommended that you use quotation marks around your adapter sequence if
you use this feature. For poly-A trimming, for example, you would write::

    cutadapt -a "A{100}" -o output.fastq input.fastq


.. _modifying-reads:

Modifying reads
===============

This section describes in which ways reads can be modified other than adapter
removal.

.. seealso::

   :ref:`Read modification order <read-modification-order>`


.. _changing-what-is-done-when-an-adapter-is-found:
.. _action:

``--action`` changes what is done when an adapter is found
----------------------------------------------------------

The ``--action`` option can be used to change what is done when an adapter match
is found in a read.

The default is ``--action=trim``, which will remove the adapter and the
sequence before or after it from the read. For 5' adapters, the adapter and
the sequence preceding it is removed. For 3' adapters, the adapter and the
sequence following it is removed. Since linked adapters are a combination of
a 5' and 3' adapter, in effect only the sequence between the 5' and the 3'
adapter matches is kept.

With ``--action=retain``, the read is trimmed, but the adapter sequence itself
is not removed. Up- and downstream sequences are removed in the same way as
for the ``trim`` action. For linked adapters, both adapter sequences are kept.

.. note::
    Because it is somewhat unclear what should happen, ``--action=retain`` can
    at the moment not be combined with ``--times`` (multiple rounds of adapter
    removal).

Use ``--action=none`` to not change the read even if there is a match.
This is useful because the statistics will still be updated as before
and because the read will still be considered "trimmed" for the read
filtering options. Combining this with ``--untrimmed-output``, for
example, can be used to copy reads without adapters to a different
file. Other read modification options, if used, may still change
the read.

Use ``--action=mask`` to write ``N`` characters to those parts of the read
that would otherwise have been removed.

Use ``--action=lowercase`` to change to lowercase those parts of the read that
would otherwise have been removed. The rest is converted to uppercase.

.. versionadded:: 3.1
    The ``retain`` action.


.. _cut-bases:

Removing a fixed number of bases
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
the ``--cut`` is applied *before* any adapter trimming.


.. _quality-trimming:

Quality trimming
----------------

The ``-q`` (or ``--quality-cutoff``) parameter can be used to trim
low-quality ends from reads. If you specify a single cutoff value, the
3' end of each read is trimmed::

    cutadapt -q 10 -o output.fastq input.fastq

For Illumina reads, this is sufficient as their quality is high at the beginning,
but degrades towards the 3' end.

It is also possible to also trim from the 5' end by specifying two
comma-separated cutoffs as *5' cutoff,3' cutoff*. For example, ::

    cutadapt -q 15,10 -o output.fastq input.fastq

will quality-trim the 5' end with a cutoff of 15 and the 3' end with a cutoff
of 10. To only trim the 5' end, use a cutoff of 0 for the 3' end, as in
``-q 15,0``.

Quality trimming is done before any adapter trimming.

For paired-end data, quality trimming is by default applied to both reads using
the same cutoff(s). Use option ``-Q`` to specify different cutoffs for R2::

    cutadapt -q 5 -Q 15,20 -o out.1.fastq -p out.2.fastq in.1.fastq in.2.fastq

To disable quality-trimming of R2, use ``-Q 0``.

By default, quality values are assumed to be encoded as
ascii(phred quality + 33). Nowadays, this should always be the case.
Some old Illumina FASTQ files encode qualities as ascii(phred quality + 64).
For those, you must add ``--quality-base=64`` to the command line.

A :ref:`description of the quality-trimming algorithm is also
available <quality-trimming-algorithm>`. The algorithm is the same as used by BWA.


.. versionadded:: 3.5
    The ``-Q`` option


.. _nextseq-trim:

Quality trimming of reads using two-color chemistry (NextSeq)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some Illumina instruments use a two-color chemistry to encode the four bases.
This includes the NextSeq and the NovaSeq. In those instruments, a
'dark cycle' (with no detected color)
encodes a ``G``. However, dark cycles also occur when sequencing "falls
off" the end of the fragment. The read then `contains a run of high-quality, but
incorrect “G” calls <https://sequencing.qcfail.com/articles/illumina-2-colour-chemistry-can-overcall-high-confidence-g-bases/>`_
at its 3' end.

Since the regular quality-trimming algorithm cannot deal with this situation,
you need to use the ``--nextseq-trim`` option::

    cutadapt --nextseq-trim=20 -o out.fastq input.fastq

This works like regular quality trimming (where one would use ``-q 20``
instead), except that the qualities of ``G`` bases are ignored.

.. versionadded:: 1.10


Shortening reads to a fixed length
----------------------------------

To shorten each read down to a certain length, use the ``--length`` option or
the short version ``-l``::

    cutadapt -l 10 -o output.fastq.gz input.fastq.gz

This shortens all reads from ``input.fastq.gz`` down to 10 bases. The removed bases
are those on the 3' end.

If you want to remove a fixed number of bases from each read, use
:ref:`the --cut option instead <cut-bases>`.


.. _modifying-read-names:

Modifying read names
--------------------

If you feel the need to modify the names of processed reads, some of the
following options may be useful.

These options exist; they are explained in more detail in the following
sections:

- ``--rename`` changes a read name according to a template.
- ``--prefix`` (or ``-x``) adds a prefix to read names.
- ``--suffix`` (or ``-y``) adds a suffix to read names.
- ``--length-tag`` updates a “length tag” such as ``length=`` with the correct read length
- ``--strip-suffix`` removes a known suffix from read names

The ``--prefix`` and ``--suffix`` options are outdated as they do not ensure that paired-end
read names remain consistent, and you should prefer to use ``--rename``.
``--prefix`` and ``--suffix`` can currently not be used together with ``--rename``.

.. _rename:
.. _read-renaming:

``--rename`` renames reads
~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``--rename`` option can be used to rename both single-end and paired-end reads.
This section describes how it can be used to rename single-end reads.

We use the following terminology: The FASTQ or FASTA record header line consists of a
*read ID* and is optionally followed by a separator (whitespace) and a *comment*.

For example, in this FASTQ header, the read ID is ``read1234`` and the comment is ``value=17``
(sequence and qualities not shown)::

    @read1234 value=17


The ``--rename`` option expects a *template string* such as
``{id} extra_info {adapter_name}`` as a parameter. It can contain regular text
and placeholders that consist of a name enclosed in curly braces (``{placeholdername}``).
The character sequence ``\t`` will be replaced by a tab character (this is currently the only
allowed escape sequence).

The read name will be set to the template string in which the placeholders are
replaced with the actual values relevant for the current read.

The following placeholders are currently available for single-end reads:

* ``{header}`` -- the full, unchanged header
* ``{id}`` -- the read ID, that is, the part of the header before the first whitespace
* ``{comment}`` -- the part of the header after the whitespace following the ID
* ``{adapter_name}`` -- the name of adapter that was found in this read or
  ``no_adapter`` if there was no adapter match. If you use ``--times`` to do
  multiple rounds of adapter matching, this is the name of the *last* found adapter.
* ``{match_sequence}`` -- the sequence of the read that matched the adapter (including
  errors). If there was no adapter match, this is set to an empty string. If you use a
  linked adapter, this is to the two matching strings, separated by a comma.
* ``{cut_prefix}`` -- the prefix removed by the ``--cut`` (or ``-u``) option (that is, when
  used with a positive length argument)
* ``{cut_suffix}`` -- the suffix removed by the ``--cut`` (or ``-u``) option (that is, when
  used with a negative length argument)
* ``{rc}`` -- this is replaced with the string ``rc`` if the read was reverse complemented.
  This only applies when :ref:`reverse complementing <reverse-complement>` was requested.
* ``\t`` -- not a placeholder, but will be replaced with the tab character.

For example, assume you have this input read in ``in.fasta``::

    >myread extra info
    ACGTAAAATTTTCCCC

Running the command ::

    cutadapt -a myadapter=TTTT -u 4 --rename='{id} barcode={cut_prefix} adapter={adapter_name} {comment}' in.fasta

Will result in this modified read::

    >myread barcode=ACGT adapter=myadapter extra info
    AAAA


.. versionadded:: 3.2

    The ``{rn}`` placeholder.

.. versionadded:: 3.3
    The ``{rc}`` placeholder.

.. versionadded:: 3.6
    The ``{match_sequence}`` placeholder.

.. versionadded:: 4.3
    The ``\t`` escape sequence.

``--rename`` also renames paired-end reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the ``--rename`` option is used with paired-end data, the template is applied
separately to both R1 and R2. That is, for R1, the placeholders are replaced with values
from R1, and for R2, the placeholders are replaced with values from R2. For example,
``{comment}`` becomes R1’s comment in R1 and it becomes R2’s comment in R2.

As another example, using ``--rename='{id} please note: {comment}'``, the paired-end reads ::

    >myread important comment
    ...

    >myread also quite important
    ...

are renamed to ::

    >myread please note: important comment
    ...

    >myread please note: also quite important
    ...

For paired-end data, the placeholder ``{rn}`` is available (“read number”),
and it is replaced with ``1`` in R1 and with ``2`` in R2.

In addition, it is possible to write a placeholder as ``{r1.placeholdername}`` or
``{r2.placeholdername}``, which always takes the replacement value from R1 or R2,
respectively.

For example, assume R1 starts with a 4 nt barcode that you want to “move” from the
sequence into the ID of both reads. You can use
``--cut=4 --rename='{id}_{r1.cut_prefix} {comment}'`` and the read pair ::

    >myread this is R1
    ACGTAAAATTTT

    >myread this is R2
    GGGGCCCC

will be changed to ::

    >myread_ACGT this is R1
    AAAATTTT

    >myread_ACGT this is R2
    GGGGCCCC

The ``{r1.placeholder}`` and ``{r2.placeholder}`` notation is available for all
placeholders except ``{rn}`` and ``{id}`` because the read ID needs to be
identical for both reads.


In general, the read IDs of R1 and R2 need to be identical. Cutadapt
enforces this when reading paired-end FASTQ files, except that it allows a single trailing
"1" or "2" as the only difference between the read IDs. This allows for read IDs ending in
``/1`` and ``/2`` (some old formats are like this) or ``.1`` and ``.2`` (``fastq-dump``
produces this).
If you use ``--rename``, Cutadapt will also enforce this when *writing* paired-end reads.

.. versionadded:: 3.2
    The ``--rename`` option


Other read name modifications
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use ``-y`` (or its alias ``--suffix``) to append a text to read names. The given string can
contain the placeholder ``{name}``, which will be replaced with the name of the
adapter found in that read. For example, writing ::

    cutadapt -a adapter1=ACGT -y ' we found {name}' input.fastq

changes a read named ``read1`` to ``read1 we found adapter1`` if the adapter
``ACGT`` was found.

The option ``-x`` (and its alias ``--prefix``) work the same, except that the text
is added in front of the read name. For both options, spaces need to be
specified explicitly, as in the above example. If no adapter was found in a
read, the text ``no_adapter`` is inserted for ``{name}``.

We recommend that you no longer use the ``-x``/``--prefix``/``-y``/``--suffix``
options and use ``--rename`` instead, which is more general.

In order to remove a suffix of each read name, use ``--strip-suffix``.

Some old 454 read files contain the length of the read in the name::

    >read1 length=17
    ACGTACGTACAAAAAAA

If you want to update this to the correct length after trimming, use the option
``--length-tag``. In this example, this would be ``--length-tag 'length='``.
After trimming, the read would perhaps look like this::

    >read1 length=10
    ACGTACGTAC


.. _filtering:

Filtering reads
===============

By default, all processed reads, no matter whether they were trimmed or not,
are written to the output file specified by the ``-o`` option (or to standard
output if ``-o`` was not provided). For paired-end reads, the second read in a
pair is always written to the file specified by the ``-p`` option.

The options described here make it possible to filter reads by either discarding
them entirely or by redirecting them to other files. When redirecting reads,
the basic rule is that *each read is written to at most one file*. You cannot
write reads to more than one output file.

Filters are applied to *all* processed reads, no matter whether they have been
modified by adapter- or quality trimming.

``--minimum-length LENGTH`` or ``-m LENGTH``
    Discard processed reads that are shorter than LENGTH.

    If you do not use this option, reads that have a length of zero (empty
    reads) are kept in the output. Some downstream tools may have problems
    with zero-length sequences. In that case, specify at least ``-m 1``.

``--too-short-output FILE``
    Instead of discarding the reads that are too short according to ``-m``,
    write them to *FILE* (in FASTA/FASTQ format).

``--maximum-length LENGTH`` or ``-M LENGTH``
    Discard processed reads that are longer than LENGTH.

``--too-long-output FILE``
    Instead of discarding reads that are too long (according to ``-M``),
    write them to *FILE* (in FASTA/FASTQ format).

``--untrimmed-output FILE``
    Write all reads without adapters to *FILE* (in FASTA/FASTQ format) instead
    of writing them to the regular output file.

``--discard-trimmed``
    Discard reads in which an adapter was found.

``--discard-untrimmed``
    Discard reads in which *no* adapter was found. This has the same effect as
    specifying ``--untrimmed-output /dev/null``.

The options ``--too-short-output`` and ``--too-long-output`` are applied first.
This means, for example, that a read that is too long will never end up in the
``--untrimmed-output`` file when ``--too-long-output`` was given, no matter
whether it was trimmed or not.

The options ``--untrimmed-output``, ``--discard-trimmed`` and ``-discard-untrimmed``
are mutually exclusive.

The following filtering options do not have a corresponding option for redirecting
reads. They always discard those reads for which the filtering criterion applies.

``--max-n COUNT_or_FRACTION``
    Discard reads with more than COUNT ``N`` bases. If ``COUNT_or_FRACTION`` is
    a number between 0 and 1, it is interpreted as a fraction of the read length

``--max-expected-errors ERRORS`` or ``--max-ee ERRORS``
    Discard reads with more than ERRORS expected errors. The number of expected
    errors is computed as described in
    `Edgar et al. (2015) <https://academic.oup.com/bioinformatics/article/31/21/3476/194979>`_,
    (Section 2.2).

``--discard-casava``
    Discard reads that did not pass CASAVA filtering. Illumina’s CASAVA pipeline in
    version 1.8 adds an *is_filtered* header field to each read. Specifying this
    option, the reads that did not pass filtering (these are the reads that have
    a ``Y`` for *is_filtered*) will be discarded. Reads for which the header cannot
    be recognized are kept.


.. _paired-end:

Trimming paired-end reads
=========================

Cutadapt supports trimming of paired-end reads. To enable this, provide two
input files and a second output file with the ``-p`` option (this is the short
form of ``--paired-output``). This is the basic command line syntax::

    cutadapt -a ADAPTER_FWD -A ADAPTER_REV -o out.1.fastq -p out.2.fastq reads.1.fastq reads.2.fastq

Here, the input reads are in ``reads.1.fastq`` and ``reads.2.fastq``, and the
result will be written to ``out.1.fastq`` and ``out.2.fastq``.

In paired-end mode, the options ``-a``, ``-b``, ``-g`` and ``-u`` that also
exist in single-end mode are applied to the forward reads only. To modify
the reverse read, these options have uppercase versions ``-A``, ``-B``,
``-G`` and ``-U`` that work just like their counterparts.
In the example above, ``ADAPTER_FWD`` will therefore be trimmed from the
forward reads and ``ADAPTER_REV`` from the reverse reads.

====================== ===========================
Single-end/R1 option   Corresponding option for R2
====================== ===========================
``--adapter``, ``-a``  ``-A``
``--front``, ``-g``    ``-G``
``--anywhere``, ``-b`` ``-B``
``--cut``, ``-u``      ``-U``
``--output``, ``-o``   ``--paired-output``, ``-p``
====================== ===========================

In paired-end mode, Cutadapt checks whether the input files are
properly paired. An error is raised if one of the files contains more reads than
the other or if the read names in the two files do not match. The read name
comparison ignores a trailing ``/1`` or ``/2`` to allow processing some old
Illumina paired-end files.

In some cases, it works to run Cutadapt twice in single-end mode on the input
files, but we recommend against it as this skips the consistency checks that
Cutadapt can do otherwise.

Also, as soon as you start to use one of the filtering options that discard
reads, it is mandatory you process both files at the same time to make sure that the
output files are kept synchronized. If a read is removed from one of the files,
Cutadapt will always ensure that it is also removed from the other file.

The following command-line options are applied to *both* reads:

* ``-q`` (along with ``--quality-base``)
* ``--times`` applies to all the adapters given
* ``--trim-n``
* ``--action``
* ``--length``
* ``--length-tag``
* ``--prefix``, ``--suffix``

The following limitations still exist:

* The ``--info-file``, ``--rest-file`` and ``--wildcard-file`` options write out
  information only from the first read.


.. _filtering-paired:

Filtering paired-end reads
--------------------------

The :ref:`filtering options listed above <filtering>` can also be used when
trimming paired-end data.

Importantly, Cutadapt *always discards both reads of a pair* if it determines
that the pair should be discarded. This ensures that the reads in the output
files are in sync. (If you don’t want or need this, you can run Cutadapt
separately on the R1 and R2 files.)

The same applies also to the options that redirect reads to other files if they
fulfill a filtering criterion, such as
``--too-short-output``/``--too-short-paired-output``. That is, the reads are
always sent in pairs to these alternative output files.

The ``--pair-filter`` option determines how to combine the filters for
R1 and R2 into a single decision about the read pair.

The default is ``--pair-filter=any``, which means that a read pair is discarded
(or redirected) if at least *one of* the reads (R1 or R2) fulfills the filtering criterion.
As an example, if option ``--minimum-length=20`` is used and paired-end data is
processed, a read pair is discarded if at least one of the reads is shorter than
20 nt.

With ``--pair-filter=both``, you can require that filtering criteria must apply
to *both* reads in order for a read pair to be discarded.

Finally, ``--pair-filter=first`` will make a decision about the read pair
by inspecting whether the filtering criterion applies to the first read,
ignoring the second read.

The following table describes the effect for some filtering options.

+----------------------------+------------------------------------------------+-----------------------------------------+
| Filtering option           | With ``--pair-filter=any``, the pair           | With ``--pair-filter=both``, the pair   |
|                            | is discarded if ...                            | is discarded if ...                     |
+============================+================================================+=========================================+
| ``--minimum-length``       | one of the reads is too short                  | both reads are too short                |
+----------------------------+------------------------------------------------+-----------------------------------------+
| ``--maximum-length``       | one of the reads is too long                   | both reads are too long                 |
+----------------------------+------------------------------------------------+-----------------------------------------+
| ``--discard-trimmed``      | one of the reads contains an adapter           | both reads contain an adapter           |
+----------------------------+------------------------------------------------+-----------------------------------------+
| ``--discard-untrimmed``    | one of the reads does not contain an adapter   | both reads do not contain an adapter    |
+----------------------------+------------------------------------------------+-----------------------------------------+
| ``--max-n``                | one of the reads contains too many ``N`` bases | both reads contain too many ``N`` bases |
+----------------------------+------------------------------------------------+-----------------------------------------+

There is currently no way to change the pair-filter mode for each filter individually.

.. note::

    As an exception, when you specify adapters *only* for R1 (``-a``/``-g``/``-b``) or *only* for
    R2 (``-A``/``-G``/``-B``), then the ``--pair-filter`` mode for ``--discard-untrimmed`` is
    forced to be ``both`` (and accordingly, also for the ``--untrimmed-(paired-)output`` options).

    Otherwise, with the default ``--pair-filter=any`` setting, all pairs would be considered
    untrimmed because it would always be the case that one of the reads in the pair does not contain
    an adapter.

    The pair-filter mode for the other filtering options, such as ``--minimum-length``, is
    not overridden in the same way and remains ``any`` unless changed explicitly with the
    ``--pair-filter`` option.

These are the paired-end specific filtering and output options:

``--minimum-length LENGTH1:LENGTH2`` or ``-m LENGTH1:LENGTH2``
    When trimming paired-end reads, the minimum lengths for R1 and R2 can be specified
    separately by separating them with a colon (``:``). If the colon syntax is not used,
    the same minimum length applies to both reads, as discussed above. Also, one of the
    values can be omitted to impose no restrictions. For example, with ``-m 17:``,
    the length of R1 must be at least 17, but the length of R2 is ignored.

``--maximum-length LENGTH1:LENGTH2`` or ``-M LENGTH1:LENGTH2``
    Maximum lengths can also be specified separately, see the explanation of ``-m`` above.

``--paired-output FILE`` or ``-p FILE``
    Write the second read of each processed pair to *FILE* (in FASTA/FASTQ
    format).

``--untrimmed-paired-output FILE``
    Used together with ``--untrimmed-output``. The second read in a pair is
    written to this file when the processed pair was *not* trimmed.

``--too-short-paired-output FILE``
    Write the second read in a pair to this file if pair is too short. Use
    together with ``--too-short-output``.

``--too-long-paired-output FILE``
    Write the second read in a pair to this file if pair is too long. Use
    together with ``--too-long-output``.

``--pair-filter=(any|both|first)``
    Which of the reads in a paired-end read have to match the filtering
    criterion in order for it to be filtered.


Note that the option names can be abbreviated as long as it is clear which
option is meant (unique prefix). For example, instead of ``--untrimmed-output``
and ``--untrimmed-paired-output``, you can write ``--untrimmed-o`` and
``--untrimmed-p``.

.. versionadded:: 1.18
    ``--pair-filter=first``


.. _paired-adapters:
.. _pair-adapters:

Paired adapters
---------------

When processing paired-end data, Cutadapt has two sets of adapters to work with: The ones that
are to be found and removed in the forward read (R1), specified with ``-a``/``-g``/``-b``,
and the ones to be found and removed in the reverse read (R2), specified with ``-A``/``-G``/``-B``.

Normally, the program looks at the R1 and R2 reads independently. That is, the best matching R1
adapter is removed from R1 and the best matching R2 adapter is removed from R2.

To change this, the option ``--pair-adapters`` can be used. It causes each R1 adapter to be
paired up with its corresponding R2 adapters. The first R1 adapter will be paired up with the first
R2 adapter, and so on. The adapters are then always removed in pairs from a read pair. It is an
error if the number of provided adapters is not identical for the R1 and R2 sets.

Example::

    cutadapt --pair-adapters -a AAAAA -a GGGG -A CCCCC -a TTTT -o out.1.fastq -p out.2.fastq in.1.fastq in.2.fastq

Here, the adapter pairs are (``AAAAA``, ``CCCCC``) and (``GGGG``, ``TTTT``). That is, paired-end
reads will only be trimmed if either

* ``AAAAA`` is found in R1 *and* ``CCCCC`` is found in R2,
* or ``GGGG`` is found in R1 *and* ``TTTT`` is found in R2.

There is one limitation of the algorithm at the moment: The program looks for the best-matching R1 adapter
first and then checks whether the corresponding R2 adapter can be found. If not, the read pair
remains unchanged. However, it is in theory possible that a different R1 adapter that does not
fit as well would have a partner that *can* be found. Some read pairs may therefore remain untrimmed.

This option was added to help with
:ref:`demultiplexing Illumina unique dual indices (UDIs) <unique-dual-indices>`.

.. versionadded:: 2.1


Interleaved paired-end reads
----------------------------

Cutadapt supports reading and writing paired-end reads from a single FASTQ file
in which the entries for the first and second read from each pair alternate.
The first read in each pair comes before the second. This is called “interleaved”
format. Enable this file format by adding the ``--interleaved`` option to the
command-line. Then, if you provide only a single file where usually two would be
expected, reads are automatically read or written interleaved.

For example, to read interleaved from ``reads.fastq`` and to write interleaved to ``trimmed.fastq``::

    cutadapt --interleaved -q 20 -a ACGT -A TGCA -o trimmed.fastq reads.fastq

In the following example, the input ``reads.fastq`` is interleaved, but output is
written to two files ``trimmed.1.fastq`` and ``trimmed.2.fastq``::

    cutadapt --interleaved -q 20 -a ACGT -A TGCA -o trimmed.1.fastq -p trimmed.2.fastq reads.fastq

Reading two-file input and writing interleaved is also possible by providing
a second input file::

    cutadapt --interleaved -q 20 -a ACGT -A TGCA -o trimmed.1.fastq reads.1.fastq reads.2.fastq

The following options also supported interleaved output::

  * ``--untrimmed-output`` (omit ``--untrimmed-paired-output``)
  * ``--too-short-output`` (omit ``--too-short-paired-output``)
  * ``--too-long-output`` (omit ``--too-long-paired-output``)

If you omit ``--interleaved`` but trim paired-end files, the above options must be used in pairs.

Cutadapt will detect if an input file is not properly interleaved by checking
whether read names match and whether the file contains an even number of entries.


Trimming paired-end reads separately
------------------------------------

.. warning::

    Trimming paired-end data in this way is not recommended as it
    bypasses all paired-end error-checking, such as checking whether
    the number of reads is the same in both files. You should use
    the normal paired-end trimming mode with the ``-o``/``--p``
    options described above.

If you do not use any of the filtering options that discard reads, such
as ``--discard``, ``--minimum-length`` or ``--maximum-length``, you can run
Cutadapt on each file separately::

    cutadapt -a ADAPTER_FWD -o trimmed.1.fastq.gz reads1.fastq.gz
    cutadapt -a ADAPTER_REV -o trimmed.2.fastq.gz reads2.fastq.gz


You can use the options that are listed under 'Additional modifications'
in Cutadapt's help output without problems. For example, if you want to
quality-trim the first read in each pair with a threshold of 10, and the
second read in each pair with a threshold of 15, then the commands could
be::

    cutadapt -q 10 -a ADAPTER_FWD -o trimmed.1.fastq reads1.fastq
    cutadapt -q 15 -a ADAPTER_REV -o trimmed.2.fastq reads2.fastq


.. note::

    Previous Cutadapt versions (up to 1.18) had a “legacy mode” that was
    activated under certain conditions and in which the read-modifying
    options such as ``-q`` would only apply to the forward/R1 reads.
    This mode no longer exists.


.. _multiple-adapters:

Multiple adapters
=================

It is possible to specify more than one adapter sequence by using the options
``-a``, ``-b`` and ``-g`` more than once. Any combination is allowed, such as
five ``-a`` adapters and two ``-g`` adapters. Each read will be searched for
all given adapters, but **only the best matching adapter is removed**. (But it
is possible to :ref:`trim more than one adapter from each
read <more-than-one>`). This is how a command may look to trim one of two
possible 3' adapters::

    cutadapt -a TGAGACACGCA -a AGGCACACAGGG -o output.fastq input.fastq

The adapter sequences can also be read from a FASTA file. Instead of giving an
explicit adapter sequence, you need to write ``file:`` followed by the name of
the FASTA file::

    cutadapt -a file:adapters.fasta -o output.fastq input.fastq

All of the sequences in the file ``adapters.fasta`` will be used as 3'
adapters. The other adapter options ``-b`` and ``-g`` also support this.

With ``-g``, you can also write ``-g ^file:adapters.fasta`` to specify that
all adapters read from ``adapters.fasta`` should be anchored.

Similarly, with ``-a``, you can also write ``-a file$:adapters.fasta`` to
anchor all adapters to the 3' end.

The ``file:`` syntax can be combined with the regular way of specifying an
adapter. But no matter how you specify multiple adapter sequences, remember
that only the best matching adapter is trimmed from each read.

When Cutadapt has multiple adapter sequences to work with, either specified
explicitly on the command line or via a FASTA file, it decides in the
following way which adapter should be trimmed:

* All given adapter sequences are matched to the read.
* Adapter matches where the overlap length (see the ``-O`` parameter) is too
  small or where the error rate is too high (``-e``) are removed from further
  consideration.
* Among the remaining matches, the one with the largest alignment score
  is chosen.
* If there is a tie, the first adapter wins. The order of adapters is the order
  in which they are given on the command line or in which they are found in the
  FASTA file.

If your adapter sequences are all similar and differ only by a variable barcode
sequence, you can use a single adapter sequence instead that
:ref:`contains wildcard characters <wildcards>`.

If you want to search for a combination of a 5' and a 3' adapter, you may want
to provide them as a single so-called :ref:`"linked adapter" <linked-adapters>`
instead.


.. versionadded:: 4.1
   Ability to anchor 5’ adapters from an external file with ``-g ^file:``

.. versionadded:: 4.3
   Ability to anchor 3' adapters from an external file with ``-a file$:``


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

    ADAPTERADAPTERADAPTERmysequence

To be on the safe side, you assume that there are at most five copies of the
adapter sequence. This command can be used to trim the reads correctly::

    cutadapt -g ^ADAPTER -n 5 -o output.fastq.gz input.fastq.gz

To search for a combination of a 5' and a 3' adapter, have a look
at the :ref:`support for "linked adapters" <linked-adapters>` instead, which
works better for that particular case because it is allows you to require that
the 3' adapter is trimmed only when the 5' adapter also occurs, and it cannot
happen that the same adapter is trimmed twice.

Before Cutadapt supported linked adapters, the ``--times`` option was the
recommended way to search for 5'/3' linked adapters. For completeness, we
describe how it was done. For example, when the 5' adapter is *FIRST* and the
3' adapter is *SECOND*, then the read could look like this::

    FIRSTmysequenceSECOND

That is, the sequence of interest is framed by the 5' and the 3' adapter. The
following command would be used to trim such a read::

    cutadapt -g ^FIRST -a SECOND -n 2 ...


.. _demultiplexing:

Demultiplexing
==============

Cutadapt supports demultiplexing, which means that reads are written to different
output files depending on which adapter was found in them. To use this, include
the string ``{name}`` in the name of the output file and :ref:`give each adapter
a name <named-adapters>`.
The path is then interpreted as a template and each trimmed read is written
to the path in which ``{name}`` is replaced with the name of the adapter that
was found in the read. Reads in which no adapter was found will be written to a
file in which ``{name}`` is replaced with ``unknown``.

Example::

    cutadapt -a one=TATA -a two=GCGC -o trimmed-{name}.fastq.gz input.fastq.gz

This command will create the three files ``demulti-one.fastq.gz``,
``demulti-two.fastq.gz`` and ``demulti-unknown.fastq.gz``.

More realistically, your “adapters” would actually be barcode sequences that you
will want to :ref:`provide in a FASTA file <multiple-adapters>`. Here is a
made-up example for such a ``barcodes.fasta`` file::

    >barcode01
    TTAAGGCC
    >barcode02
    TAGCTAGC
    >barcode03
    ATGATGAT

Since our barcodes are located at the 5’ end of the R1 read, we use the ``-g``
option to provide Cutadapt with the adapter sequences, as in ``-g ^file:barcodes.fasta``.
Also, we prefix the ``^file:`` with the ``^`` character to specify that we want to
:ref:`anchor the 5’ adapters <anchored-5adapters>`.

Since these barcode sequences have a length of 8 and the default maximum error
rate is 10%, Cutadapt would by default not allow any errors when matching them
(a single error would result in an error rate of 1/8=12.5%). We therefore use
``-e 1`` to allow one error.

Here is the final command::

    cutadapt -e 1 -g ^file:barcodes.fasta -o "trimmed-{name}.fastq.gz" input.fastq.gz

Demultiplexing is also supported for paired-end data if you provide the ``{name}`` template
in both output file names (``-o`` and ``-p``). Example::

    cutadapt -e 1 -g ^file:barcodes.fasta -o trimmed-{name}.1.fastq.gz -p trimmed-{name}.2.fastq.gz input.1.fastq.gz input.2.fastq.gz

Paired-end demultiplexing always uses the adapter matches of the *first* read to decide where a
read should be written. If adapters for read 2 are given (``-A``/``-G``), they are detected and
removed as normal, but these matches do not influence where the read pair is written. This is
to ensure that read 1 and read 2 are always synchronized.

To demultiplex using a barcode that is located on read 2,
you can "cheat" and swap the roles of R1 and R2 for
both the input and output files ::

    cutadapt -e 1 -g ^file:barcodes.fasta -o trimmed-{name}.2.fastq.gz -p trimmed-{name}.1.fastq.gz input.2.fastq.gz input.1.fastq.gz

If you do this in a script or pipeline, it may be a good idea to add a comment to clarify that
this reversal of R1 and R2 is intended.

More advice on demultiplexing:

* You can use ``--untrimmed-output`` to change the name of the output file that receives the
  untrimmed reads (those in which no barcode could be found).
* Similarly, you can use ``--untrimmed-paired-output`` to change the name of the output file that
  receives the untrimmed R2 reads.
* If you want to demultiplex, but keep the barcode in the reads, use the option ``--action=none``.


.. _combinatorial-demultiplexing:

Demultiplexing paired-end reads with combinatorial dual indexes
---------------------------------------------------------------

`Illumina’s combinatorial dual indexing strategy <https://support.illumina.com/bulletins/2018/08/understanding-unique-dual-indexes--udi--and-associated-library-p.html>`_
uses a set of indexed adapters on R1 and another one on R2. Unlike unique dual indexes (UDI)
(described on the same page) all combinations of indexes are possible.

For demultiplexing this type of data ("combinatorial demultiplexing"), it is necessary to write each
read pair to an output file depending on the adapters found on R1 *and* R2.

Doing this is similar to doing normal demultiplexing as described above, but you need
to use ``{name1}`` and ``{name2}`` in both output file name templates. For example::

    cutadapt \
        -e 0.15 --no-indels \
        -g ^file:barcodes_fwd.fasta \
        -G ^file:barcodes_rev.fasta \
        -o {name1}-{name2}.1.fastq.gz -p {name1}-{name2}.2.fastq.gz \
        input.1.fastq.gz input.2.fastq.gz

The ``{name1}`` will be replaced with the name of the best-matching R1 adapter and ``{name2}`` will
be replaced with the name of the best-matching R2 adapter.

If there was no match of an R1 adapter, ``{name1}`` is set to "unknown", and if there is no match of
an R2 adapter, ``{name2}`` is set to "unknown". To discard read pairs for which one or both adapters
could not be found, use ``--discard-untrimmed``.

The ``--untrimmed-output`` and ``--untrimmed-paired-output`` options cannot be used.

Read the :ref:`demultiplexing <demultiplexing>` section for how to choose the error rate etc.
Also, the tips below about how to speed up demultiplexing apply even with combinatorial
demultiplexing.

When doing the above, you will end up with lots of files named ``first-second.x.fastq.gz``, where
*first* is the name of the first indexed adapter and *second* is the name of the second indexed
adapter, and *x* is 1 or 2. Each indexed adapter combination may correspond to a sample name and
you may want to name your files according to the sample name, not the name of the adapters.
Cutadapt does not have built-in functionality to achieve this, but you can use an external
tool such as ``mmv`` (“multiple move”). First, create a list of patterns in ``patterns.txt``::

    fwdindex1-revindex1.[12].fastq.gz sampleA.#1.fastq.gz
    fwdindex1-revindex2.[12].fastq.gz sampleB.#1.fastq.gz
    fwdindex1-revindex3.[12].fastq.gz sampleC.#1.fastq.gz
    fwdindex2-revindex1.[12].fastq.gz sampleD.#1.fastq.gz
    fwdindex2-revindex2.[12].fastq.gz sampleE.#1.fastq.gz
    ...

Here, *fwdindex1*/*revindex1* etc. are the names of indexes, and *sampleA* etc.
are your sample names. Then rename all files at once with ::

    mmv < patterns.txt


.. versionadded:: 2.4


.. _paired-adapters-dual-indices:
.. _unique-dual-indices:

Demultiplexing unique dual indices
----------------------------------

`Illumina’s unique dual indexing (UDI) scheme <https://support.illumina.com/bulletins/2018/08/understanding-unique-dual-indexes--udi--and-associated-library-p.html>`_
(“non-redundant indexing”) uses 96 unique i5 indices and 96 unique i7
indices, which are only used in pairs. That is, the first i5 index is always
used with the first i7 index and so on.

To demultiplex this type of data, the
:ref:`--pair-adapters option <pair-adapters>` needs to be used. Example::

    cutadapt -j 8 -e 1 --no-indels --pair-adapters -g ^file:i5indices.fasta -G ^file:i7indices.fasta -o 'demultiplexed-{name}_R1.fastq.gz' -p 'demultiplexed-{name}_R2.fastq.gz' input.R1.fastq.gz input.R2.fastq.gz


.. note::
    If the adapters do not come in pairs, but all combinations are possible, use
    :ref:`combinatorial demultiplexing <combinatorial-demultiplexing>`.


.. _speed-up-demultiplexing:

Speeding up demultiplexing
--------------------------

Finding many adapters/barcodes simultaneously (which is what demultiplexing in Cutadapt is about),
can be sped up tremendously by using the right options since Cutadapt will then be able to create an
index of the barcode sequences instead of checking for each barcode separately. Currently, the
following conditions need to be met in order for index creation to be enabled:

* The barcodes/adapters must be anchored:
  For 5’ adapters, use ``-g ^ADAPTER`` or ``-g ^file:adapters.fasta``.
  For 3’ adapters, use ``-a ADAPTER$`` or ``-a file$:adapters.fasta``.
* The maximum error rate (``-e``) must be set such that at most 2 errors are allowed:
  Use ``-e 0``, ``-e 1`` or ``-e 2``.
* No IUPAC wildcards must be used in the barcode/adapter. Also, you cannot use the option
  ``--match-read-wildcards``.

An index will be built for all the adapters that fulfill these criteria if there are at least two
of them. You can provide additional adapters/barcodes, and they will just not be included in the
index. Whether an index is created or not should not affect the results, only how fast you get them.

To see whether an index is created, look for a message like this in the first few lines of
Cutadapt’s output::

    Building index of 23 adapters ...

Hopefully some of the above restrictions will be lifted in the future.

.. versionadded:: 1.15
   Demultiplexing of paired-end data.

.. versionadded:: 2.0
   Added ability to use an index of adapters for speeding up demultiplexing

.. versionadded:: 3.0
   An index can be built even when indels are allowed (that is, ``--no-indels``
   is no longer required).


Demultiplexing paired-end reads in mixed orientation
----------------------------------------------------

For some protocols, the barcode will be located either on R1 or on R2
depending on the orientation in which the DNA fragment was sequenced.

For example, the read layout could be either this ::

    R1: barcode-forwardprimer-sequence  R2: reverseprimer-sequence

or this ::

    R1: reverseprimer-sequence  R2: barcode-forwardprimer-sequence

To demultiplex such data with Cutadapt, choose one of the orientations first and
demultiplex the reads as if only that existed in the data, using a command like this ::

    cutadapt -g ^file:barcodes.fasta \
        -o round1-{name}.R1.fastq.gz \
        -p round1-{name}.R2.fastq.gz \
        R1.fastq.gz R2.fastq.gz

Then all the read pairs in which no barcode could be found will end up in
``round1-unknown.R1.fastq.gz`` and ``round1-unknown.R2.fastq.gz``. This will
also include the pairs in which the barcode was not actually in R1, but in R2. To
demultiplex these reads as well, run Cutadapt a second time with those “unknown”
files as input, but also reverse the roles of R1 and R2 ::

    cutadapt -g ^file:barcodes.fasta \
        -o round2-{name}.R2.fastq.gz \
        -p round2-{name}.R1.fastq.gz \
        round1-unknown.R2.fastq.gz round1-unknown.R1.fastq.gz


.. _truseq:

Illumina TruSeq
===============

Illumina makes their adapter sequences available in the
`Illumina Adapter Sequences Document <https://support.illumina.com/downloads/illumina-adapter-sequences-document-1000000002694.html>`_.

As an example for how to use that information with Cutadapt, we show
how to trim TruSeq adapters. The document gives the adapter sequence
for read 1 as ``AGATCGGAAGAGCACACGTCTGAACTCCAGTCA`` and for read 2
as ``AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT``. When using Cutadapt, this
means you should trim your paired-end data as follows::

    cutadapt \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -o trimmed.R1.fastq.gz -p trimmed.R2.fastq.gz \
        reads.R1.fastq.gz reads.R2.fastq.gz

See also the :ref:`section about paired-end adapter trimming above <paired-end>`.

Keep in mind that Cutadapt removes the adapter that it finds and also the sequence
following it, so even if the actual adapter sequence that is used in a protocol
is longer than that (and possibly contains a variable index), it is sufficient to
specify a prefix of the sequence(s).

.. note::
   Previous versions of this document also recommended using ``AGATCGGAAGAGC``
   as adapter sequence for both read 1 and read 2, but you should avoid doing so
   as that sequence occurs multiple times in the human genome.

To understand the structure of Illumina libraries and what the i5, i7, P5, P7
sequences are, see
`this overview <https://teichlab.github.io/scg_lib_structs/methods_html/Illumina.html>`_.

Some older information is also available in the document `Illumina TruSeq Adapters
De-Mystified <http://tucf-genomics.tufts.edu/documents/protocols/TUCF_Understanding_Illumina_TruSeq_Adapters.pdf>`_,
but keep in mind that it does not cover newer protocols.


Under some circumstances, you may want to consider not trimming adapters at all.
For example, a good library prepared for exome, genome or transcriptome
sequencing should contain very few reads with adapters anyway. Also, some read
mapping programs including BWA-MEM and STAR will soft-clip bases at the 3' ends
of reads that do not match the reference, which will take care of adapters
implicitly.


.. _warnbase:

Warning about incomplete adapter sequences
------------------------------------------

Sometimes Cutadapt’s report ends with these lines::

    WARNING:
        One or more of your adapter sequences may be incomplete.
        Please see the detailed output above.

Further up, you’ll see a message like this::

    Bases preceding removed adapters:
      A: 95.5%
      C: 1.0%
      G: 1.6%
      T: 1.6%
      none/other: 0.3%
    WARNING:
        The adapter is preceded by "A" extremely often.
        The provided adapter sequence may be incomplete.
        To fix the problem, add "A" to the beginning of the adapter sequence.

This means that in 95.5% of the cases in which an adapter was removed from a
read, the base coming *before* that was an ``A``. If your DNA fragments are
not random, such as in amplicon sequencing, then this is to be expected and
the warning can be ignored. If the DNA fragments are supposed to be random,
then the message may be genuine: The adapter sequence may be incomplete and
should include an additional ``A`` in the beginning.

This warning exists because some documents list the Illumina TruSeq adapters
as starting with ``GATCGGA...``. While that is technically correct, the
library preparation actually results in an additional ``A`` before that
sequence, which also needs to be removed. See the :ref:`previous
section <truseq>` for the correct sequence.

.. _n-bases:
.. _dealing-with-ns:

Dealing with ``N`` bases
========================

Cutadapt supports the following options to deal with ``N`` bases in your reads:

``--max-n COUNT``
    Discard reads containing more than *COUNT* ``N`` bases. A fractional *COUNT*
    between 0 and 1 can also be given and will be treated as the proportion of
    maximally allowed ``N`` bases in the read. For example, ``--max-n 0``
    removes all reads that contain any ``N`` bases.

``--trim-n``
    Remove flanking ``N`` bases from each read. That is, a read such as this::

        NNACGTACGTNNNN

    Is trimmed to just ``ACGTACGT``. This option is applied *after* adapter
    trimming. If you want to get rid of ``N`` bases before adapter removal, use
    quality trimming: ``N`` bases typically also have a low quality value
    associated with them.

.. _cutadapt-s-output:

Cutadapt's output
=================

Reporting
---------

Cutadapt will by default print a full report after it has finished processing
the reads. To suppress all output except error messages, use the option
``--quiet``.

The report type can be changed to a one-line summary with the option
``--report=minimal``. The output will be a tab-separated table (tsv) with one
header row and one row of content. Here is an example::

    $ cutadapt --report=minimal -a ... -m 20 -q 10 -o ... -p ... in.[12].fastq.gz
    status in_reads in_bp     too_short too_long too_many_n out_reads w/adapters qualtrim_bp out_bp w/adapters2 qualtrim2_bp out2_bp
    OK     1000000  202000000 24827     0        0          975173    28968      1674222     97441426 0 0 98492473

This is the meaning of each column:

=============== ==========================================================
Column heading  Explanation
=============== ==========================================================
status          Incomplete adapter warning (``OK`` or ``WARN``)
in_reads        Number of processed reads (read pairs for paired-end)
in_bp           Number of processed basepairs
too_short       Number of reads/read pairs that were too short
too_long        Number of reads/read pairs that were too long
too_many_n      Number of reads/read pairs that contained too many ``N``
out_reads       Number of reads written
w/adapters      Number of reads containing at least one adapter
qualtrim_bp     Number of bases removed from R1 reads by quality trimming
out_bp          Number of bases written to R1 reads
w/adapters2     Number of R2 reads containing at least one adapter
qualtrim2_bp    Number of bases removed from R2 reads by quality trimming
out2_bp         Number of bases written
=============== ==========================================================

The last three fields are omitted for single-end data.

.. versionadded:: 1.18


How to read the report
----------------------

After every run, Cutadapt prints out per-adapter statistics. The output
starts with something like this::

    Sequence: 'ACGTACGTACGTTAGCTAGC'; Length: 20; Trimmed: 2402 times.

If option ``--revcomp`` was used,
this line will additionally contain something like ``Reverse-complemented:
984 times``. This describes how many times of the 2402 total times the
adapter was found on the reverse complement of the read.

The next piece of information is this::

    No. of allowed errors:
    0-7 bp: 0; 8-15 bp: 1; 16-20 bp: 2

The adapter, as was shown above, has a length of 20
characters. We are using a custom error rate of 0.12. What this
implies is shown above: Matches up to a length of 7 bp are allowed to
have no errors. Matches of lengths 8-15 bp are allowd to have 1 error
and matches of length 16 or more can have 2 errors. See also :ref:`the section about
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
    7       13      0.3     0       1 12
    8       31      0.1     1       0 31
    ...
    100     397     0.0     3       358 36 3

The first row tells us the following: Three bases were removed in 140
reads; randomly, one would expect this to occur 156.2 times; the maximum
number of errors at that match length is 0 (this is actually redundant
since we know already that no errors are allowed at lengths 0-7 bp).

The last column shows the number of reads that had 0, 1, 2 ... errors.
In the last row, for example, 358 reads matched the adapter with zero
errors, 36 with 1 error, and 3 matched with 2 errors.

In the row for length 7 is an apparent anomaly, where the max.err column
is 0 and yet we have 31 reads matching with 1 error. This is because the
matches are actually contributed by alignments to the first 8 bases of
the adapter with one deletion, so 7 bases are removed but the error
cut-off applied is for length 8.

The "expect" column gives only a rough estimate of the number of
sequences that is expected to match randomly, but it can help to
estimate whether the matches that were found are true adapter matches
or if they are due to chance. At lengths 6, for example, only 2.4
reads are expected, but 35 do match, which hints that most of these
matches are due to actual adapters.
For slightly more accurate estimates, you can provide the correct
GC content (as a percentage) of your reads with the option
``--gc-content``. The default is ``--gc-content=50``.

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


JSON report
-----------

With ``--json=filename.cutadapt.json``, a report in JSON format is written to the given file.

We strongly recommend that you use the ``.cutadapt.json`` file name extension for this file for
easier discoverability by log-parsing tools such as `MultiQC <https://multiqc.info>`_.

See the :ref:`description of the JSON report file format <json-report-format>`.

.. versionadded:: 3.5


.. _info-file:

Info file
---------

When the ``--info-file=info.tsv`` command-line parameter is given, detailed
information about where adapters were found in each read are written
to the given text file as tab-separated values.

See the :ref:`description of the info file format <info-file-format>`.
