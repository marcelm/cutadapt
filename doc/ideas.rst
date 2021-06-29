Ideas/To Do
-----------

This is a rather unsorted list of features that would be nice to have, of
things that could be improved in the source code, and of possible algorithmic
improvements.

- show average error rate
- length histogram
- ``--detect`` prints out best guess which of the given adapters is the correct one
- warn when given adapter sequence contains non-IUPAC characters


Specifying adapters
~~~~~~~~~~~~~~~~~~~

Allow something such as ``-a ADAP$TER`` or ``-a ADAPTER$NNN``.
This would be a way to specify less strict anchoring.

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


Available letters for command-line options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Lowercase letters: i, k, s, w
* Uppercase letters: C, D, E, F, H, I, J, K, L, P, R, S, T, V, W
* Deprecated, could be re-used: c, d, t
* Planned/reserved: Q (paired-end quality trimming), V (alias for --version)
