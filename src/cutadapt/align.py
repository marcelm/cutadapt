__all__ = [
    'Aligner',
    'PrefixComparer',
    'SuffixComparer',
    'hamming_sphere',
    'hamming_environment',
    'edit_environment',
    'edit_distance',
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


def edit_distance(s: str, t: str):
    """
    Return the edit distance between the strings s and t.
    The edit distance is the sum of the numbers of insertions, deletions,
    and mismatches that is minimally necessary to transform one string
    into the other.
    """
    m = len(s)  # index i
    n = len(t)  # index j
    costs = list(range(m + 1))

    for j in range(1, n + 1):
        prev = costs[0]
        costs[0] += 1
        for i in range(1, m + 1):
            match = int(s[i - 1] == t[j - 1])
            c = min(
                prev + 1 - match,
                costs[i] + 1,
                costs[i - 1] + 1,
            )
            prev = costs[i]
            costs[i] = c
    return costs[-1]


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


def naive_edit_environment(s: str, k: int):
    """
    Apply all possible edits up to edit distance k to string s.
    A string may be returned more than once.
    """
    yield s
    if k == 0:
        return
    for s in naive_edit_environment(s, k - 1):
        n = len(s)
        for ch in 'ACGT':
            # all insertions
            for i in range(n+1):
                x = s[:i] + ch + s[i:]
                assert len(x) == len(s) + 1
                yield x
            # all substitutions
            for i in range(n):
                x = s[:i] + ch + s[i+1:]
                assert len(x) == len(s)
                yield x
        # all deletions
        for i in range(n):
            x = s[:i] + s[i+1:]
            assert len(x) == len(s) - 1
            yield x


def edit_environment(s: str, k: int):
    """
    Find all strings t for which the edit distance between s and t is at most k,
    assuming the alphabet is A, C, G, T.

    Yield tuples (t, e, m), where e is the edit distance between s and t and
    m is the number of matches in the optimal alignment.
    """
    n = len(s)
    alphabet = "TGCA"
    work_stack = [(
        "",
        list(range(n + 1)),
        [0] * (n + 1),
    )]
    while work_stack:
        # t is the current prefix
        # costs is a row at index len(t) in the DP matrix
        # matches is a row in the corresponding matrix of the no. of matches
        t, costs, matches = work_stack.pop()

        # The row is the last row of the DP matrix for aligning t against s
        i = len(t)
        if costs[-1] <= k:
            # The costs of an optimal alignment of t against s are at most k,
            # so t is within the edit environment.
            yield t, costs[-1], matches[-1]

        if i == n + k:
            # Last row reached
            continue

        # Runtime heuristic: The entries in the DP matrix cannot get lower
        # in subsequent rows, so don’t try longer suffixs if all entries are
        # greater than k.
        if min(costs) > k:
            continue

        # compute next row in DP matrix for all characters of the alphabet
        for ch in alphabet:
            # create a new DP matrix row for each character of the alphabet
            next_costs = [0] * (n + 1)
            next_costs[0] = len(t) + 1
            next_matches = [0] * (n + 1)
            for j in range(1, n + 1):
                match = 0 if s[j - 1] == ch else 1
                assert j > 0

                diag = costs[j-1] + match
                left = next_costs[j-1] + 1
                up = costs[j] + 1
                if diag <= left and diag <= up:
                    c, m = diag, matches[j-1] + (1 - match)
                elif left <= up:
                    c, m = left, next_matches[j-1]
                else:
                    c, m = up, matches[j]
                next_costs[j] = c
                next_matches[j] = m
            work_stack.append((t + ch, next_costs, next_matches))
