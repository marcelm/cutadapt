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
- convert adapter to lowercase
- warn when given adapter sequence contains non-IUPAC characters


Specifying adapters
-------------------

The idea is to deprecate the ``-b`` and ``-g`` parameters. Only ``-a`` is used
with a special syntax for each adapter type. This makes it a bit easier to add
new adapter types in the feature.

.. csv-table::

    back,``-a ADAPTER``,``-a ADAPTER`` or ``-a ...ADAPTER``
    suffix,``-a ADAPTER$``,``-a ...ADAPTER$``
    front,``-g ADAPTER``,``-a ADAPTER...``
    prefix,``-g ^ADAPTER``,``-a ^ADAPTER...``
    anywhere,``-b ADAPTER``, ``-a ...ADAPTER...`` ???
    paired,(not implemented),``-a ADAPTER...ADAPTER`` or ``-a ^ADAPTER...ADAPTER``

Or add only ``-a ADAPTER...`` as an alias for ``-g ^ADAPTER`` and
``-a ...ADAPTER`` as an alias for ``-a ADAPTER``.

Another idea: Allow something such as ``-a ADAP$TER`` or ``-a ADAPTER$NNN``.
This would be a way to specify less strict anchoring.


Paired-end trimming
-------------------

* Could also use a paired-end read merger, then remove adapters with -a and -g
* Should minimum overlap be sum of the two overlaps in each read?


Single-letter command-line options
----------------------------------

Remaining characters: All uppercase letters except A, B, G, M, N, O
Lowercase letters: i, j, k, l, s, w
