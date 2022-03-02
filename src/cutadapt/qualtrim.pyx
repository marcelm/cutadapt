# cython: profile=False, emit_code_comments=False, language_level=3
"""
Quality trimming.
"""

cdef extern from *:
    unsigned char * PyUnicode_1BYTE_DATA(object o)
    int PyUnicode_KIND(object o)
    int PyUnicode_1BYTE_KIND


cdef class HasNoQualities(Exception):
    pass


def quality_trim_index(str qualities, int cutoff_front, int cutoff_back, int base=33):
    """
    Find the positions at which to trim low-quality ends from a nucleotide sequence.
    Return tuple (start, stop) that indicates the good-quality segment.

    Qualities are assumed to be ASCII-encoded as chr(qual + base).

    The algorithm is the same as the one used by BWA within the function
    'bwa_trim_read':
    - Subtract the cutoff value from all qualities.
    - Compute partial sums from all indices to the end of the sequence.
    - Trim sequence at the index at which the sum is minimal.
    """
    if qualities is None:
        raise HasNoQualities("Cannot do quality trimming when no qualities are available")
    cdef:
        int s
        int max_qual
        int n = len(qualities)
        int stop = n
        int start = 0
        int i
        char* qual

    if not PyUnicode_KIND(qualities) == PyUnicode_1BYTE_KIND:
        raise ValueError("Quality data is not ASCII")
    qual = <char *>PyUnicode_1BYTE_DATA(qualities)

    # find trim position for 5' end
    s = 0
    max_qual = 0
    for i in range(n):
        s += cutoff_front - (qual[i] - base)
        if s < 0:
            break
        if s > max_qual:
            max_qual = s
            start = i + 1

    # same for 3' end
    max_qual = 0
    s = 0
    for i in reversed(range(n)):
        s += cutoff_back - (qual[i] - base)
        if s < 0:
            break
        if s > max_qual:
            max_qual = s
            stop = i
    if start >= stop:
        start, stop = 0, 0
    return (start, stop)


def nextseq_trim_index(sequence, int cutoff, int base=33):
    """
    Variant of the above quality trimming routine that works on NextSeq data.
    With Illumina NextSeq, bases are encoded with two colors. 'No color' (a
    dark cycle) usually means that a 'G' was sequenced, but that also occurs
    when sequencing falls off the end of the fragment. The read then contains
    a run of high-quality G bases in the end.

    This routine works as the one above, but counts qualities belonging to 'G'
    bases as being equal to cutoff - 1.
    """
    bases = sequence.sequence
    qualities = sequence.qualities
    if qualities is None:
        raise HasNoQualities()
    cdef:
        int s = 0
        int max_qual = 0
        int max_i
        int i, q
        char* qual

    if not PyUnicode_KIND(qualities) == PyUnicode_1BYTE_KIND:
        raise ValueError("Quality data is not ASCII")
    qual = <char *>PyUnicode_1BYTE_DATA(qualities)

    s = 0
    max_qual = 0
    max_i = len(qualities)
    for i in reversed(range(max_i)):
        q = qual[i] - base
        if bases[i] == 'G':
            q = cutoff - 1
        s += cutoff - q
        if s < 0:
            break
        if s > max_qual:
            max_qual = s
            max_i = i
    return max_i


def expected_errors(str qualities, int base=33):
    """
    Return the number of expected errors (as double) from a read’s
    qualities.

    This uses the formula in Edgar et al. (2015),
    see Section 2.2 in <https://academic.oup.com/bioinformatics/article/31/21/3476/194979>.

    qualities -- ASCII-encoded qualities (chr(qual + base))
    """
    cdef:
        int i, q
        bytes quals = qualities.encode()
        char* cq = quals
        double e = 0.0

    for i in range(len(qualities)):
        q = cq[i] - base
        e += 10 ** (-q / 10)
    return e
