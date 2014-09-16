(Frequently) Asked Questions
============================

How does the `--overlap` parameter work?
----------------------------------------

The minimum overlap length specified by the ``--overlap`` parameter helps to
reduce trimming of randomly matching adapters. The process is as follows: First,
the adapter is matched to the read. The trailing bases of the read that match
the initial bases of the adapter are the *overlap*. The assumption is that this
is the part of the read that should be discarded. However, the last few bases of
the read and the first few bases of the adapter will often match purely by
chance. For example, in every fourth read (on average) the last base will match
the first base of the adapter. As a consequence, many reads would be trimmed by
just a few bases, simply because the adapter matches by chance. In order to
reduce this effect, there's the ``--overlap`` option. If the length of the
overlap is below the threshold provided with that option, the read will not be
trimmed.

The default value for ``--overlap`` is 3.


