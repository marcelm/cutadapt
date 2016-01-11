Ideas/To Do
===========

This is a rather unsorted list of features that would be nice to have, of
things that could be improved in the source code, and of possible algorithmic
improvements.

- show average error rate
- In colorspace and probably also for Illumina data, gapped alignment
  is not necessary
- ``--progress``
- run pylint, pychecker
- length histogram
- check whether input is FASTQ although -f fasta is given
- search for adapters in the order in which they are given on the
  command line
- more tests for the alignment algorithm
- deprecate ``--rest-file``
- ``--detect`` prints out best guess which of the given adapters is the correct one
- alignment algorithm: make a 'banded' version
- it seems the str.find optimization isn't very helpful. In any case, it should be
  moved into the Aligner class.
- allow to remove not the adapter itself, but the sequence before or after it
- instead of trimming, convert adapter to lowercase
- warn when given adapter sequence contains non-IUPAC characters
- try multithreading again, this time use os.pipe() or 0mq
- extensible file type detection
- the --times setting should be an attribute of Adapter

Specifying adapters
-------------------

The idea is to deprecate the ``-b``,  ``-g`` and ``-u`` parameters. Only ``-a``
is used with a special syntax for each adapter type. This makes it a bit easier
to add new adapter types in the feature.

.. csv-table::

    back,``-a ADAPTER``,``-a ADAPTER`` or ``-a ...ADAPTER``
    suffix,``-a ADAPTER$``,``-a ...ADAPTER$``
    front,``-g ADAPTER``,``-a ADAPTER...``
    prefix,``-g ^ADAPTER``,``-a ^ADAPTER...`` (or have anchoring by default?)
    anywhere,``-b ADAPTER``, ``-a ...ADAPTER...`` ???
    unconditional,``-u +10``,``-a 10...`` (collides with colorspace)
    unconditional,``-u -10``,``-a ...10$``
    linked,(not implemented),``-a ADAPTER...ADAPTER`` or ``-a ^ADAPTER...ADAPTER``

Or add only ``-a ADAPTER...`` as an alias for ``-g ^ADAPTER`` and
``-a ...ADAPTER`` as an alias for ``-a ADAPTER``.

The ``...`` would be equivalent to ``N*`` as in regular expressions.

Another idea: Allow something such as ``-a ADAP$TER`` or ``-a ADAPTER$NNN``.
This would be a way to specify less strict anchoring.

Make it possible to specify that the rightmost or leftmost match should be
picked. Default right now: Leftmost, even for -g adapters.

Allow ``N{3,10}`` as in regular expressions (for a variable-length sequence).

Use parentheses to specify the part of the sequence that should be kept:

* ``-a (...)ADAPTER`` (default)
* ``-a (...ADAPTER)`` (default)
* ``-a ADAPTER(...)`` (default)
* ``-a (ADAPTER...)`` (??)

Or, specify the part that should be removed:

    ``-a ...(ADAPTER...)``
    ``-a ...ADAPTER(...)``
    ``-a (ADAPTER)...``

Model somehow all the flags that exist for semiglobal alignment. For start of the adapter:

* Start of adapter can be degraded or not
* Bases are allowed to be before adapter or not

Not degraded and no bases before allowed = anchored.
Degraded and bases before allowed = regular 5'

By default, the 5' end should be anchored, the 3' end not.

* ``-a ADAPTER...`` → not degraded, no bases before allowed
* ``-a N*ADAPTER...`` → not degraded, bases before allowed
* ``-a ADAPTER^...`` → degraded, no bases before allowed
* ``-a N*ADAPTER^...`` → degraded, bases before allowed
* ``-a ...ADAPTER`` → degraded, bases after allowed
* ``-a ...ADAPTER$`` → not degraded, no bases after allowed



Paired-end trimming
-------------------

* Could also use a paired-end read merger, then remove adapters with -a and -g

Available/used letters for command-line options
-----------------------------------------------

* Remaining characters: All uppercase letters except A, B, G, M, N, O, U
* Lowercase letters: i, j, k, l, s, w
* Planned/reserved: Q (paired-end quality trimming), j (multithreading)
