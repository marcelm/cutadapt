__all__ = [
    "EndSkip",
    "Aligner",
    "PrefixComparer",
    "SuffixComparer",
    "hamming_sphere",
    "hamming_environment",
    "edit_environment",
    "edit_distance",
]

from enum import IntFlag
from typing import Iterator, Tuple

from cutadapt._align import (
    Aligner,
    PrefixComparer,
    SuffixComparer,
    hamming_sphere,
    edit_environment,
)


class EndSkip(IntFlag):
    """
    Flags for the Aligner that indicate which ends of reference or query may be skipped at
    no cost. Setting all four flags at the same time results in semiglobal alignment.
    """

    REFERENCE_START = 1  # a prefix of the reference may be skipped at no cost
    QUERY_START = 2  # a prefix of the query may be skipped at no cost
    REFERENCE_END = 4  # a suffix of the reference may be skipeed at no cost
    QUERY_STOP = 8  # a suffix of the query may be skipeed at no cost
    SEMIGLOBAL = 15  # all of the above


def edit_distance(s: str, t: str) -> int:
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


def hamming_environment(s: str, k: int) -> Iterator[Tuple[str, int, int]]:
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


def naive_edit_environment(s: str, k: int) -> Iterator[str]:
    """
    Apply all possible edits up to edit distance k to string s
    and yield the resulting strings.
    A string may be returned more than once.
    """
    yield s
    if k == 0:
        return
    for s in naive_edit_environment(s, k - 1):
        n = len(s)
        for ch in "ACGT":
            for i in range(n):
                prefix = s[:i] + ch
                yield prefix + s[i:]  # insertion
                yield prefix + s[i + 1 :]  # substitution
            yield s + ch  # insertion into final position
        # all deletions
        for i in range(n):
            yield s[:i] + s[i + 1 :]


def py_edit_environment(s: str, k: int) -> Iterator[Tuple[str, int, int]]:
    """
    Find all strings t for which the edit distance between s and t is at most k,
    assuming the alphabet is A, C, G, T.

    Yield tuples (t, e, score), where e is the edit distance between s and t and
    score is the score of the optimal alignment.
    """
    rate = k / len(s) if s else 0
    aligner = Aligner(s, max_error_rate=rate, flags=0, min_overlap=len(s))
    seen = set()
    for t in naive_edit_environment(s, k):
        if t in seen:
            continue
        seen.add(t)
        result = aligner.locate(t)
        score, errors = result[-2:]  # type: ignore
        yield t, errors, score


def slow_edit_environment(s: str, k: int) -> Iterator[Tuple[str, int, int]]:
    """
    Find all strings t for which the edit distance between s and t is at most k,
    assuming the alphabet is A, C, G, T.

    Yield tuples (t, e, m), where e is the edit distance between s and t and
    m is the number of matches in the optimal alignment.

    This is slow and only used in testing.
    """
    n = len(s)
    alphabet = "TGCA"
    work_stack = [
        (
            "",
            list(range(n + 1)),
            [0] * (n + 1),
        )
    ]
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

                diag = costs[j - 1] + match
                left = next_costs[j - 1] + 1
                up = costs[j] + 1
                if diag <= left and diag <= up:
                    c, m = diag, matches[j - 1] + (1 - match)
                elif left <= up:
                    c, m = left, next_matches[j - 1]
                else:
                    c, m = up, matches[j]
                next_costs[j] = c
                next_matches[j] = m
            work_stack.append((t + ch, next_costs, next_matches))
