Ideas/To Do
===========

This is a rather unsorted list of features that would be nice to have, of
things that could be improved in the source code, and of possible algorithmic
improvements.


- show average error rate
- In colorspace and probably also for Illumina data, gapped alignment
  is not necessary
- use ``str.format`` instead of ``%``
- allow to change scores at runtime (using command-line parameters)
- ``--progress``
- run pylint, pychecker
- length histogram
- check whether input is FASTQ although -f fasta is given
- close on StopIteration
- search for adapters in the order in which they are given on the
  command line
- number reads as given on command line
- more tests for the alignment algorithm
- allow adapters in ``file:`` to be anchored
- deprecate ``--rest-file``
- ``--detect`` prints out best guess which of the given adapters is the correct one
- alignment algorithm: make a 'banded' version


Specifying adapters
-------------------

The idea is to do away with the ``-b`` and ``-g`` parameters. Only ``-a`` is used
with a special syntax for each adapter type. This makes it a bit easier to add
new adapter types in the feature.

.. csv-table::

    back,``-a ADAPTER``,``-a ADAPTER`` or ``-a ...ADAPTER``
    front,``-g ADAPTER``,``-a ADAPTER...``
    prefix,``-g ^ADAPTER``,``-a ^ADAPTER...``
    anywhere,``-b ADAPTER``, ``-a ...ADAPTER...`` ???
    paired,(not implemented),``-a ADAPTER...ADAPTER`` or ``-a ^ADAPTER...ADAPTER``

Or add only ``-a ADAPTER...`` as an alias for ``-g ^ADAPTER`` and
``-a ...ADAPTER`` as an alias for ``-a ADAPTER``.
