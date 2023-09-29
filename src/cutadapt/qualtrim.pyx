# cython: profile=False, emit_code_comments=False, language_level=3
"""
Quality trimming.
"""
from cpython.unicode cimport PyUnicode_GET_LENGTH
from libc.stdint cimport uint8_t

cdef extern from *:
    unsigned char * PyUnicode_1BYTE_DATA(object o)
    void *PyUnicode_DATA(object o)
    bint PyUnicode_IS_COMPACT_ASCII(object o)
    int PyUnicode_KIND(object o)
    int PyUnicode_1BYTE_KIND

cdef extern from "expected_errors.h":
    float expected_errors_from_phreds(const uint8_t *phreds, size_t phreds_length, uint8_t base)

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


def poly_a_trim_index(str s, bint revcomp = False):
    """
    Return start index of poly-A tail

    If revcomp is True, return end of poly-T head instead.

    Poly-A tails shorter than 3 are ignored.
    """
    if not PyUnicode_KIND(s) == PyUnicode_1BYTE_KIND:
        raise ValueError("Sequence is not ASCII")
    cdef:
        char* s_ptr = <char *>PyUnicode_1BYTE_DATA(s)
        int n = len(s)
        int best_score = 0
        int score = 0
        int i
        char c
        int errors = 0
        int best_index

    if revcomp:
        best_index = 0
        for i in range(n):
            if s_ptr[i] == b"T":
                score += 1
            else:
                score -= 2
                errors += 1

            if score > best_score and errors * 5 <= i + 1:  # max error rate 0.2
                best_score = score
                best_index = i + 1
        if best_index < 3:
            best_index = 0
    else:
        best_index = n
        for i in reversed(range(n)):
            if s_ptr[i] == b"A":
                score += 1
            else:
                score -= 2
                errors += 1

            if score > best_score and errors * 5 <= n - i:  # max error rate 0.2
                best_score = score
                best_index = i
        if best_index > n - 3:
            best_index = n
    return best_index


def expected_errors(str qualities, uint8_t base=33):
    """
    Return the number of expected errors (as double) from a read’s
    qualities.

    This uses the formula in Edgar et al. (2015),
    see Section 2.2 in <https://academic.oup.com/bioinformatics/article/31/21/3476/194979>.

    qualities -- ASCII-encoded qualities (chr(qual + base))
    """
    if not PyUnicode_IS_COMPACT_ASCII(qualities):
        raise ValueError(f"Quality string contains non-ASCII values: {qualities}")
    cdef:
        uint8_t *quals = <uint8_t *>PyUnicode_DATA(qualities)
        size_t qual_length = PyUnicode_GET_LENGTH(qualities)
        double e = expected_errors_from_phreds(quals, qual_length, base)

    if e < 0.0:
        for q in qualities:
            if ord(q) < base or ord(q) > 126:
                raise ValueError(f"Not a valid phred value {ord(q)} for character {q}")
    return e
