=================
Algorithm details
=================


.. _adapter-alignment-algorithm:

Adapter alignment algorithm
===========================

Since the publication of the `EMBnet journal application note about
Cutadapt <http://dx.doi.org/10.14806/ej.17.1.200>`_, the alignment algorithm
used for finding adapters has changed significantly. An overview of this new
algorithm is given in this section. An even more detailed description is
available in Chapter 2 of my PhD thesis `Algorithms and tools for the analysis
of high-throughput DNA sequencing data <http://hdl.handle.net/2003/31824>`_.

The algorithm is based on *semiglobal alignment*, also called *free-shift*,
*ends-free* or *overlap* alignment. In a regular (global) alignment, the
two sequences are compared from end to end and all differences occuring over
that length are counted. In semiglobal alignment, the sequences are allowed to
freely shift relative to each other and differences are only penalized in the
overlapping region between them::

      FANTASTIC
   ELEFANT

The prefix ``ELE`` and the suffix ``ASTIC`` do not have a counterpart in the
respective other row, but this is not counted as an error. The overlap ``FANT``
has a length of four characters.

Traditionally, *alignment scores* are used to find an optimal overlap aligment:
This means that the scoring function assigns a positive value to matches,
while mismatches, insertions and deletions get negative values. The optimal
alignment is then the one that has the maximal total score. Usage of scores
has the disadvantage that they are not at all intuitive: What does a total score
of *x* mean? Is that good or bad? How should a threshold be chosen in order to
avoid finding alignments with too many errors?

For Cutadapt, the adapter alignment algorithm primarily uses *unit costs* instead.
This means that mismatches, insertions and deletions are counted as one error, which
is easier to understand and allows to specify a single parameter for the
algorithm (the maximum error rate) in order to describe how many errors are
acceptable.

There is a problem with this: When using costs instead of scores, we would like
to minimize the total costs in order to find an optimal alignment. But then the
best alignment would always be the one in which the two sequences do not overlap
at all! This would be correct, but meaningless for the purpose of finding an
adapter sequence.

The optimization criteria are therefore a bit different. The basic idea is to
consider the alignment optimal that maximizes the overlap between the two
sequences, as long as the allowed error rate is not exceeded.

Conceptually, the procedure is as follows:

1. Consider all possible overlaps between the two sequences and compute an
   alignment for each, minimizing the total number of errors in each one.
2. Keep only those alignments that do not exceed the specified maximum error
   rate.
3. Then, keep only those alignments that have a maximal number of matches
   (that is, there is no alignment with more matches). (Note: This has been
   changed, see the section below for an update.)
4. If there are multiple alignments with the same number of matches, then keep
   only those that have the smallest error rate.
5. If there are still multiple candidates left, choose the alignment that starts
   at the leftmost position within the read.

In Step 1, the different adapter types are taken into account: Only those
overlaps that are actually allowed by the adapter type are actually considered.


.. _algorithm-indel-scores:

Alignment algorithm changes in Cutadapt 4
=========================================

The above algorithm has been tweaked slightly in Cutadapt 4.
The main problem was that the idea of maximizing the number of matches
(criterion 3 in the section above) sometimes leads to unintuitive results.

For example, the previous algorithm would prefer an alignment such as this one::

    CCAGTCCTTTCCTGAGAGT   Read
    ||||||||   ||
    CCAGTCCT---CT         5' adapter

This alignment was considered to be the best one because it contains 10 matches,
which is the maximum possible.
The three consecutive deletions are ignored when making that decision.
To the user, the unexpected result is visible because the read would end up as
``GAGAGT`` after trimming.

With the tuned algorithm, the alignment is more sensible::

    CCAGTCCTTTCCTGAGAGT   Read
    ||||||||X|
    CCAGTCCTCT            5' adapter

The trimmed read is now ``CCTGAGAGT``, which is what one would likely expect.

The alignment algorithm in Cutadapt can perhaps now be described as
a *hybrid* algorithm that uses both edit distance and score:

- Edit distance is used to fill out the dynamic programming matrix.
  Conceptually, this can be seen as computing the edit distance for all
  possible overlaps between the read and the adapter.
  We need to use the edit distance as optimization criterion at this
  stage because we want to be able to let the user provide a maximum
  error rate (``-e``).
  Also, using edit distance (that is, unit costs) allows using some
  optimizations while filling in the matrix (Ukkonen’s trick).
- A second matrix with scores is filled in simultaneously.
  The value in a cell is the score of the edit-distance-based alignment,
  the score is not used as optimization criterion.
- Finally, the score is used to decide which of the overlaps between
  read and adapter is the best one.
  (This means looking into the last row and column of the score matrix.)

The score function is currently: match: +1, mismatch: -1, indel: -2

A second change in the alignment algorithm is relevant if there are
multiple adapter occurrences in a read (such as adapter dimers).
With the new algorithm, leftmost (earlier) adapter occurrences are now more
reliably preferred even if a later match has fewer errors.

Here are two examples from the SRR452441 dataset (R1 only),
trimmed with the standard Illumina adapter.
The top row shows the alignment as found by the previous algorithm,
the middle row shows the sequencing read,
and the last row shows the alignment as found by the updated algorithm. ::

    @SRR452441.2151945
                                             AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  Previous alignment
                                             ||||||||||||||||||||||||||||||||||
    -GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCACACGAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGCACACGAATCTCGTATGCCGTCTTCT
    X|||||||||||||||||||||||||||||||||
    AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  New alignment

Previously the read was trimmed to the first 40 bases,
now the earlier, nearly full-length occurrence is taken into account,
and the read is empty after trimming. ::

    @SRR452441.2157038
                                    AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  Previous alignment
                                    ||||||||||||||||||||||||||||||||||
    -GATCGGAAGAGCACACGTCTGAACTCCAGTCAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGCACACGAATCTCGTATGCCGTCTTCTGCTTGAAAA
    X||||||||||||||||||||||||||||||||X
    AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  New alignment

Only very few reads should be affected by the above changes (in SRR452441, which has 2.2 million
reads, only four reads were trimmed differently).
In those cases where it matters, however, there should now be fewer surprises.


.. _quality-trimming-algorithm:

Quality trimming algorithm
==========================

The trimming algorithm implemented in Cutadapt is the same as the one used by
BWA, but applied to both
ends of the read in turn (if requested). That is: Subtract the given cutoff
from all qualities; compute partial sums from all indices to the end of the
sequence; cut the sequence at the index at which the sum is minimal. If both
ends are to be trimmed, repeat this for the other end.

The basic idea is to remove all bases starting from the end of the read whose
quality is smaller than the given threshold. This is refined a bit by allowing
some good-quality bases among the bad-quality ones. In the following example,
we assume that the 3' end is to be quality-trimmed.

Assume you use a threshold of 10 and have these quality values:

42, 40, 26, 27, 8, 7, 11, 4, 2, 3

Subtracting the threshold gives:

32, 30, 16, 17, -2, -3, 1, -6, -8, -7

Then sum up the numbers, starting from the end (partial sums). Stop early if
the sum is greater than zero:

(70), (38), 8, -8, -25, -23, -20, -21, -15, -7

The numbers in parentheses are not computed (because 8 is greater than zero),
but shown here for completeness. The position of the minimum (-25) is used as
the trimming position. Therefore, the read is trimmed to the first four bases,
which have quality values 42, 40, 26, 27.
