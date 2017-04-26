=============
Recipes (FAQ)
=============

This section gives answers to frequently asked questions. It shows you how to
get cutadapt to do what you want it to do!


Avoid internal adapter matches
------------------------------

To force matches to be at the end of the read and thus avoiding internal
adapter matches, append a few ``X`` characters to the adapter sequence, like
this: ``-a TACGGCATXXX``. The ``X`` is counted as a mismatch and will force the
match to be at the end. Just make sure that there are more ``X`` characters than
the length of the adapter times the error rate. This is not the same as an
anchored 3' adapter since partial matches are still allowed.


Remove more than one adapter
----------------------------

If you want to remove a 5' and 3' adapter at the same time, :ref:`use the
support for linked adapters <linked-adapters>`.

If your situation is different, for example, when you have many 5' adapters
but only one 3' adapter, then you have two options.

First, you can specify the adapters and also ``--times=2`` (or the short
version ``-n 2``). For example::

	cutadapt -g ^TTAAGGCC -g ^AAGCTTA -a TACGGACT -n 2 -o output.fastq input.fastq

This instructs cutadapt to run two rounds of adapter finding and removal. That
means that, after the first round and only when an adapter was actually found,
another round is performed. In both rounds, all given adapters are searched and
removed. The problem is that it could happen that one adapter is found twice (so
the 3' adapter, for example, could be removed twice).

The second option is to not use the ``-n`` option, but to run cutadapt twice,
first removing one adapter and then the other. It is easiest if you use a pipe
as in this example::

	cutadapt -g ^TTAAGGCC -g ^AAGCTTA input.fastq | cutadapt -a TACGGACT - > output.fastq


Trim poly-A tails
-----------------

If you want to trim a poly-A tail from the 3' end of your reads, use the 3'
adapter type (``-a``) with an adapter sequence of many repeated ``A``
nucleotides. Starting with version 1.8 of cutadapt, you can use the
following notation to specify a sequence that consists of 100 ``A``::

	cutadapt -a "A{100}" -o output.fastq input.fastq

This also works when there are sequencing errors in the poly-A tail. So this
read ::

	TACGTACGTACGTACGAAATAAAAAAAAAAA

will be trimmed to::

	TACGTACGTACGTACG

If for some reason you would like to use a shorter sequence of ``A``, you can
do so: The matching algorithm always picks the leftmost match that it can find,
so cutadapt will do the right thing even when the tail has more ``A`` than you
used in the adapter sequence. However, sequencing errors may result in shorter
matches than desired. For example, using ``-a "A{10}"``, the read above (where
the ``AAAT`` is followed by eleven ``A``) would be trimmed to::

	TACGTACGTACGTACGAAAT

Depending on your application, perhaps a variant of ``-a A{10}N{90}`` is an
alternative, forcing the match to be located as much to the left as possible,
while still allowing for non-``A`` bases towards the end of the read.


Trim a fixed number of bases after adapter trimming
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


Trimming (amplicon-) primers from both ends of paired-end reads
---------------------------------------------------------------

If you want to remove primer sequences that flank your sequence of
interest, you should use a :ref:`"linked adapter" <linked-adapters>`
to remove them. If you have paired-end data (with R1 and R2), you
can correctly trim both R1 and R2 by using linked adapters for both
R1 and R2. Here is how to do this.

The full DNA fragment that is put on the sequencer looks like this
(looking only at the forward strand):

   5' sequencing primer -- forward primer -- sequence of interest -- reverse complement of reverse primer -- reverse complement of 3' sequencing primer

Since sequencing of R1 starts after the 5' sequencing primer, R1 will
start with the forward primer and then continue into the sequence of
interest and into the two primers to the right of it, depending on
the read length and how long the sequence of interest is. For R1,
the linked adapter option that needs to be used is therefore ::

    -a FWDPRIMER...RCREVPRIMER

where ``FWDPRIMER`` needs to be replaced with the sequence of your
forward primer and ``RCREVPRIMER`` with the reverse complement of
the reverse primer. The three dots ``...`` need to be entered
as they are -- they tell cutadapt that this is a linked adapter
with a 5' and a 3' part.

Sequencing of R2 starts before the 3' sequencing primer and
proceeds along the reverse-complementary strand. For the correct
linked adapter, the sequences from above therefore need to be
swapped and reverse-complemented::

    -A REVPRIMER...RCFWDPRIMER

The uppercase ``-A`` specifies that this option is
meant to work on R2. Similar to above, ``REVPRIMER`` is
the sequence of the reverse primer and ``RCFWDPRIMER`` is the
reverse-complement of the forward primer. Note that cutadapt
does not reverse-complement any sequences of its own; you
will have to do that yourself.

Finally, you may want to filter the trimmed read pairs.
Use ``--discard-untrimmed`` to throw away all read pairs in
which R1 doesn’t start with ``FWDPRIMER`` or in which R2
does not start with ``REVPRIMER``.

A note on how the filtering works: In linked adapters, by default
the first part (before the ``...``) is anchored. Anchored
sequences *must* occur. If they don’t, then the other sequence
(after the ``...``) is not even searched for and the entire
read is internally marked as “untrimmed”. This is done for both
R1 and R2 and as soon as *any* of them is marked as “untrimmed”,
the entire pair is considered to be “untrimmed”. If
``--discard-untrimmed`` is used, this means that the entire
pair is discarded if R1 or R2 are untrimmed. (Option
``--pair-filter=both`` can be used to change this to require
that *both* were marked as untrimmed.)

In summary, this is how to trim your data and discard all
read pairs that do not contain the primer sequences that
you know must be there::

    cutadapt -a FWDPRIMER...RCREVPRIMER -A REVPRIMER...RCFWDPRIMER --discard-untrimmed -o out.1.fastq.gz -p out.2.fastq.gz in.1.fastq.gz in.2.fastq.gz


Other things (unfinished)
-------------------------

* How to detect adapters
* Use cutadapt for quality-trimming only
* Use it for minimum/maximum length filtering
* Use it for conversion to FASTQ
