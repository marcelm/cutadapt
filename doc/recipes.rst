=======
Recipes
=======

For some trimming applications, the pre-defined adapter types behave differently
from what you would like to have. In this section, we show some ways in which
cutadapt can be made to behave in the desired way.

.. note:: This section is still being written.


Avoiding internal adapter matches
---------------------------------

To force matches to be at the end of the read and thus avoiding internal
adapter matches, append a few ``X`` characters to the adapter sequence, like
this: ``-a TACGGCATXXX``. The ``X`` is counted as a mismatch and will force the
match to be at the end. Just make sure that there are more ``X`` characters than
the length of the adapter times the error rate. This is not the same as an
anchored 3' adapter since partial matches are still allowed.


Removing more than one adapter
------------------------------

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


Trimming poly-A tails
---------------------

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


Other things (unfinished)
-------------------------

* How to detect adapters
* Use cutadapt for quality-trimming only
* Use it for minimum/maximum length filtering
* Use it for conversion to FASTQ
