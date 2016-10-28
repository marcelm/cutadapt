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


Other things (unfinished)
-------------------------

* How to detect adapters
* Use cutadapt for quality-trimming only
* Use it for minimum/maximum length filtering
* Use it for conversion to FASTQ
