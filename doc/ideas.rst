Ideas/To Do
===========

This is a rather unsorted list of features that would be nice to have, of
things that could be improved in the source code, and of possible algorithmic
improvements.

- show average error rate
- run pylint, pychecker
- length histogram
- ``--detect`` prints out best guess which of the given adapters is the correct one
- allow to remove not the adapter itself, but the sequence before or after it
- warn when given adapter sequence contains non-IUPAC characters
- extensible file type detection


Backwards-incompatible changes
------------------------------

- Drop ``--rest-file`` support
- Possibly drop wildcard-file support, extend info-file instead
- For non-anchored 5' adapters, find rightmost match


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
    linked,``-a ADAPTER...ADAPTER``,``-a ADAPTER...ADAPTER`` or ``-a ^ADAPTER...ADAPTER``

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

API
---

* https://github.com/marcelm/cutadapt/labels/API

Paired-end trimming
-------------------

* Could also use a paired-end read merger, then remove adapters with -a and -g

Available letters for command-line options
------------------------------------------

* Lowercase letters: i, k, s, w
* Uppercase letters: C, D, E, F, H, I, J, K, L, P, R, S, T, V, W
* Deprecated, could be re-used: c, d, t
* Planned/reserved: Q (paired-end quality trimming), V (alias for --version)
