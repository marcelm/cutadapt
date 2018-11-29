__all__ = [
    'Aligner',
    'compare_prefixes',
    'compare_suffixes',
]

from cutadapt._align import Aligner, compare_prefixes


def compare_suffixes(s1, s2, wildcard_ref=False, wildcard_query=False):
    """
    Find out whether one string is the suffix of the other one, allowing
    mismatches. Used to find an anchored 3' adapter when no indels are allowed.
    """
    s1 = s1[::-1]
    s2 = s2[::-1]
    _, length, _, _, matches, errors = compare_prefixes(s1, s2, wildcard_ref, wildcard_query)
    return (len(s1) - length, len(s1), len(s2) - length, len(s2), matches, errors)


# convenience function (to avoid having to instantiate an Aligner manually)
def locate(reference, query, max_error_rate, flags=15, wildcard_ref=False,
        wildcard_query=False, min_overlap=1):
    aligner = Aligner(
        reference,
        max_error_rate,
        min_overlap=min_overlap,
        start_in_reference=flags & 1,
        start_in_query=flags & 2,
        stop_in_reference=flags & 4,
        stop_in_query=flags & 8,
        wildcard_ref=wildcard_ref,
        wildcard_query=wildcard_query,
    )
    return aligner.locate(query)
