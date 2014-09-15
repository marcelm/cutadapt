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

How does cutadapt decide which adapter to trim when multiple adapters are provided?
-----------------------------------------------------------------------------------

When multiple adapters are provided on the command line via the ``-a``, ``-b``
or ``-g`` parameters, all adapters are first matched to the read.

Adapter matches where the overlap length is too small or where the error rate is
too high are removed from further consideration. Among the remaining matches,
the criterion for deciding which match is best is the *number of matching
bases*. If there is a tie, the first adapter wins. The order of adapters is the
order in which they are given on the command line.

Percentage identity would be another possible criterion, but the idea was to
prefer long over short matches. For that, the absolute number of matching bases
is more appropriate.

Why does the -g option delete adapters even if they occur at the end or within the read?
----------------------------------------------------------------------------------------

The only difference between the ``-a`` and ``-g`` options is that ``-g`` finds
the adapter anywhere within the read and removes everything *before* it. If you
expect the read to begin with the adapter, then add the character ``^`` before
the adapter sequence on the command line. For example::

    cutadapt -g ^ADAPTER -o output.fastq input.fastq
