=======
Recipes
=======

This section contains short how-to guides for doing certain tasks.


Remove more than one adapter
----------------------------

If you want to remove a 5' and 3' adapter at the same time, :ref:`use the
support for linked adapters <linked-adapters>`.

If your situation is different, for example, when you have many 5' adapters
but only one 3' adapter, then you have two options.

First, you can specify the adapters and also ``--times=2`` (or the short
version ``-n 2``). For example::

    cutadapt -g ^TTAAGGCC -g ^AAGCTTA -a TACGGACT -n 2 -o output.fastq input.fastq

This instructs Cutadapt to run two rounds of adapter finding and removal. That
means that, after the first round and only when an adapter was actually found,
another round is performed. In both rounds, all given adapters are searched and
removed. The problem is that it could happen that one adapter is found twice (so
the 3' adapter, for example, could be removed twice).

The second option is to not use the ``-n`` option, but to run Cutadapt twice,
first removing one adapter and then the other. It is easiest if you use a pipe
as in this example::

    cutadapt -g ^TTAAGGCC -g ^AAGCTTA input.fastq | cutadapt -a TACGGACT - > output.fastq


Separating trimmed and untrimmed reads
--------------------------------------

To send trimmed and untrimmed reads to separate output files, use the the
``--untrimmed-output`` option::

    cutadapt -a TTAAGGCC --untrimmed-output=untrimmed.fastq.gz -o trimmed.fastq.gz input.fastq.gz

For paired-end data, use also the ``--untrimmed-paired-output`` option::

    cutadapt \
      -g ^AAGGCC \
      -G ^TTGGAA \
      --untrimmed-output=untrimmed.1.fastq.gz \
      --untrimmed-paired-output=untrimmed.2.fastq.gz
      -o trimmed.2.fastq.gz \
      -p trimmed.2.fastq.gz \
      input.1.fastq.gz \
      input.2.fastq.gz


Trim poly-A tails
-----------------

If you want to trim a poly-A tail from the 3' end of your reads, use the 3'
adapter type (``-a``) with an adapter sequence of many repeated ``A``
nucleotides. Use the
following notation to specify a sequence that consists of 100 ``A``::

    cutadapt -a "A{100}" -o output.fastq input.fastq

This also works when there are sequencing errors in the poly-A tail. So this
read ::

    TACGTACGTACGTACGAAATAAAAAAAAAAA

will be trimmed to::

    TACGTACGTACGTACG

If for some reason you would like to use a shorter sequence of ``A``, you can
do so: The matching algorithm always picks the leftmost match that it can find,
so Cutadapt will do the right thing even when the tail has more ``A`` than you
used in the adapter sequence. However, sequencing errors may result in shorter
matches than desired. For example, using ``-a "A{10}"``, the read above (where
the ``AAAT`` is followed by eleven ``A``) would be trimmed to::

    TACGTACGTACGTACGAAAT

Depending on your application, perhaps a variant of ``-a A{10}N{90}`` is an
alternative, forcing the match to be located as much to the left as possible,
while still allowing for non-``A`` bases towards the end of the read.


Trim a fixed number of bases preceding each adapter
---------------------------------------------------

If the adapters you want to remove are preceded by some unknown sequence (such
as a random tag/molecular identifier), you can specify this as part of the
adapter sequence in order to remove both in one go.

For example, assume you want to trim Illumina adapters preceded by 10 bases
that you want to trim as well. Instead of this command::

    cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC ...

Use this command::

    cutadapt -O 13 -a N{10}AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC ...

The ``-O 13`` is the minimum overlap for an adapter match, where the 13 is
computed as 3 plus 10 (where 3 is the default minimum overlap and 10 is the
length of the unknown section). If you do not specify it, the adapter sequence
would match the end of every read (because ``N`` matches anything), and ten
bases would then be removed from every read.


Trimming (amplicon-) primers from paired-end reads
--------------------------------------------------

If the reads are shorter than the amplicon, use ::

    cutadapt -g ^FWDPRIMER -G ^REVPRIMER --discard-untrimmed -o out.1.fastq.gz -p out.2.fastq.gz in.1.fastq.gz in.2.fastq.gz

If the reads can be longer than the amplicon, use a :ref:`"linked adapter" <linked-adapters>`::

    cutadapt -a ^FWDPRIMER...RCREVPRIMER -A ^REVPRIMER...RCFWDPRIMER --discard-untrimmed -o out.1.fastq.gz -p out.2.fastq.gz in.1.fastq.gz in.2.fastq.gz

You need to insert your own sequences as follows. The three dots ``...`` need to be written as they
are.

FWDPRIMER
    Sequence of the forward primer

REVPRIMER
    Sequence of the reverse primer

RCFWDPRIMER
    Reverse-complemented sequence of the forward primer

RCREVPRIMER
    Reverse-complemented sequence of the reverse primer

Explanation
~~~~~~~~~~~

The full DNA fragment that is put on the sequencer looks like this
(looking only at the forward strand):

   5' sequencing primer -- forward primer -- sequence of interest -- reverse complement of reverse primer -- reverse complement of 3' sequencing primer

Since sequencing of R1 starts after the 5' sequencing primer, R1 will
start with the forward primer and then continue into the sequence of
interest and possibly into the two primers to the right of it, depending on
the read length and how long the sequence of interest is.

If the reads are sufficiently short, R1 will not extend into the reverse primer, and
R2 will not extend into the forward primer. In that case, only the forward primer on R1
and the reverse primer on R2 need to be removed::

    -g ^FWDPRIMER -G ^REVPRIMER --discard-untrimmed

If the reads are so long that they can possibly extend into the primer on the respective other side,
:ref:`linked adapters <linked-adapters>` for both R1 and R2 can be used. For R1::

    -a ^FWDPRIMER...RCREVPRIMER

Sequencing of R2 starts before the 3' sequencing primer and
proceeds along the reverse-complementary strand. For the correct
linked adapter, the sequences from above therefore need to be
swapped and reverse-complemented::

    -A ^REVPRIMER...RCFWDPRIMER

The uppercase ``-A`` specifies that this option is meant to work on R2.
Cutadapt does not reverse-complement any sequences of its own; you
will have to do that yourself.

Finally, you may want to filter the trimmed read pairs.
Option ``--discard-untrimmed`` throws away all read pairs in
which R1 doesn’t start with ``FWDPRIMER`` or in which R2
does not start with ``REVPRIMER``.

A note on how the filtering works: In linked adapters, by default
the first part (before the ``...``) is anchored. Anchored
sequences *must* occur. If they don’t, then the other sequence
(after the ``...``) is not even searched for and the entire
read is internally marked as “untrimmed”. This is done for both
R1 and R2 and as soon as any of them is marked as “untrimmed”,
the entire pair is considered to be “untrimmed”. If
``--discard-untrimmed`` is used, this means that the entire
pair is discarded if R1 or R2 are untrimmed. (Option
``--pair-filter=both`` can be used to change this to require
that *both* were marked as untrimmed.)


Piping paired-end data
----------------------

Sometimes it is necessary to run Cutadapt twice on your data. For example, when
you want to change the order in which read modification or filtering options are
applied. To simplify this, you can use Unix pipes (``|``), but this is more
difficult with paired-end data since then input and output consists of two files
each.

The solution is to interleave the paired-end data, send it over the pipe
and then de-interleave it in the other process. Here is how this looks in
principle::

    cutadapt [options] --interleaved in.1.fastq.gz in.2.fastq.gz | \
      cutadapt [options] --interleaved -o out.1.fastq.gz -p out.2.fastq.gz -

Note the ``-`` character in the second invocation to Cutadapt.


Support for concatenated compressed files
-----------------------------------------

Cutadapt supports concatenated gzip and bzip2 input files.


Check whether a FASTQ file is properly formatted
------------------------------------------------

    cutadapt -o /dev/null input.fastq

Any problems with the FASTQ file will be detected and reported.


Check whether FASTQ files are properly paired
---------------------------------------------

    cutadapt -o /dev/null -p /dev/null input.R1.fastq input.R2.fastq

Any problems with the individual FASTQ files or improperly paired reads
(mismatching read ids) will be detected and reported.

See :ref:`the requirements for properly paired reads <properly-paired-reads>`.


Rescuing single reads from paired-end reads that were filtered
--------------------------------------------------------------

When trimming and filtering paired-end reads, Cutadapt always discards entire read pairs. If you
want to keep one of the reads, you need to write the filtered read pairs to an output file and
postprocess it.

For example, assume you are using ``-m 30`` to discard too short reads. Cutadapt discards all
read pairs in which just one of the reads is too short (but see the ``--pair-filter`` option).
To recover those (individual) reads that are long enough, you can first use the
``--too-short-(paired)-output`` options to write the filtered pairs to a file, and then postprocess
those files to keep only the long enough reads.


    cutadapt -m 30 -q 20 -o out.1.fastq.gz -p out.2.fastq.gz --too-short-output=tooshort.1.fastq.gz --too-short-paired-output=tooshort.2.fastq.gz in.1.fastq.gz in.2.fastq.gz
    cutadapt -m 30 -o rescued.a.fastq.gz tooshort.1.fastq.gz
    cutadapt -m 30 -o rescued.b.fastq.gz tooshort.2.fastq.gz

The two output files ``rescued.a.fastq.gz`` and ``rescued.b.fastq.gz`` contain those individual
reads that are long enough. Note that the file names do not end in ``.1.fastq.gz`` and
``.2.fastq.gz`` to make it very clear that these files no longer contain synchronized paired-end
reads.


.. _bisulfite:

Bisulfite sequencing (RRBS)
---------------------------

When trimming reads that come from a library prepared with the RRBS (reduced
representation bisulfite sequencing) protocol, the last two 3' bases must be
removed in addition to the adapter itself. This can be achieved by using not
the adapter sequence itself, but by adding two wildcard characters to its
beginning. If the adapter sequence is ``ADAPTER``, the command for trimming
should be::

    cutadapt -a NNADAPTER -o output.fastq input.fastq

Details can be found in `Babraham bioinformatics' "Brief guide to
RRBS" <http://www.bioinformatics.babraham.ac.uk/projects/bismark/RRBS_Guide.pdf>`_.
A summary follows.

During RRBS library preparation, DNA is digested with the restriction enzyme
MspI, generating a two-base overhang on the 5' end (``CG``). MspI recognizes
the sequence ``CCGG`` and cuts
between ``C`` and ``CGG``. A double-stranded DNA fragment is cut in this way::

    5'-NNNC|CGGNNN-3'
    3'-NNNGGC|CNNN-5'

The fragment between two MspI restriction sites looks like this::

    5'-CGGNNN...NNNC-3'
      3'-CNNN...NNNGGC-5'

Before sequencing (or PCR) adapters can be ligated, the missing base positions
must be filled in with GTP and CTP::

    5'-ADAPTER-CGGNNN...NNNCcg-ADAPTER-3'
    3'-ADAPTER-gcCNNN...NNNGGC-ADAPTER-5'

The filled-in bases, marked in lowercase above, do not contain any original
methylation information, and must therefore not be used for methylation calling.
By prefixing the adapter sequence with ``NN``, the bases will be automatically
stripped during adapter trimming.


.. _file-format-conversion:

Convert FASTQ to FASTA
----------------------

Cutadapt detects the output format from the output file name extension.
Convert FASTQ to FASTA format::

    cutadapt -o output.fasta.gz input.fastq.gz

Cutadapt detects FASTA output and omits the qualities.

If output is written to standard output, no output file name is available,
so the same format as the input is used.

To force FASTA output even in this case, use the ``--fasta`` option::

    cutadapt --fasta input.fastq.gz > out.fasta


Trim qualities
--------------

Quality-trim 3' ends::

    cutadapt -q 20 -o output.fastq.gz input.fastq.gz


.. _json-jq:

Extract information from the JSON report with ``jq``
----------------------------------------------------

The :ref:`JSON report <json-report-format>` that is written when using the ``--json``
option can be read by `jq <https://stedolan.github.io/jq/>`_.

Get the number of reads (or read pairs) written::

   jq '.read_counts.output' mysample.cutadapt.json

Get the percentage of reads that contain an adapter::

   jq '.read_counts.read1_with_adapter / .read_counts.input * 100' mysample.cutadapt.json

Get how often the first adapter was found::

   jq '.adapters_read1[0].total_matches' mysample.cutadapt.json


Quickly test how Cutadapt trims a single sequence
-------------------------------------------------

Use ``echo`` to write the sequence in FASTA format, and run Cutadapt with ``--quiet``::

    echo -e ">r\nAACCGGTT" | cutadapt --quiet -a CCGGTTGGAA -

Output::

    >r
    AA
