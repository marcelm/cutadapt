__all__ = [
    'Aligner',
    'PrefixComparer',
    'SuffixComparer',
    'hamming_sphere',
    'hamming_environment',
]

from cutadapt._align import Aligner, PrefixComparer, SuffixComparer

# flags for global alignment

# The interpretation of the first flag is:
# An initial portion of seq1 may be skipped at no cost.
# This is equivalent to saying that in the alignment,
# gaps in the beginning of seq2 are free.
#
# The other flags have an equivalent meaning.
START_WITHIN_SEQ1 = 1
START_WITHIN_SEQ2 = 2
STOP_WITHIN_SEQ1 = 4
STOP_WITHIN_SEQ2 = 8

# Use this to get regular semiglobal alignment
# (all gaps in the beginning or end are free)
SEMIGLOBAL = START_WITHIN_SEQ1 | START_WITHIN_SEQ2 | STOP_WITHIN_SEQ1 | STOP_WITHIN_SEQ2


def hamming_sphere(s, k):
    """
    Yield all strings t for which the hamming distance between s and t is exactly k,
    assuming the alphabet is A, C, G, T.
    """
    assert k >= 0
    if k == 0:
        yield s
        return
    n = len(s)

    # i is the first position that is varied
    for i in range(n - k + 1):
        prefix = s[:i]
        c = s[i]
        suffix = s[i+1:]
        for ch in 'ACGT':
            if ch == c:
                continue
            for t in hamming_sphere(suffix, k - 1):
                y = prefix + ch + t
                assert len(y) == n
                yield y


def hamming_environment(s, k):
    """
    Find all strings t for which the hamming distance between s and t is at most k,
    assuming the alphabet is A, C, G, T.

    Yield tuples (t, e, m), where e is the hamming distance between s and t and
    m is the number of matches (equal to len(t) - e).
    """
    n = len(s)
    for e in range(k + 1):
        for t in hamming_sphere(s, e):
            yield t, e, n - e
