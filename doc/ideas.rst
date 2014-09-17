Ideas/To Do
===========

This is a rather unsorted list of features that would be nice to have, of
things that could be improved in the source code, and of possible algorithmic
improvements.


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
