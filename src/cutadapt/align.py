__all__ = [
    'Aligner',
    'PrefixComparer',
    'SuffixComparer',
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


# convenience function (to avoid having to instantiate an Aligner manually)
def locate(reference, query, max_error_rate, flags=SEMIGLOBAL, wildcard_ref=False,
        wildcard_query=False, min_overlap=1):
    aligner = Aligner(reference, max_error_rate, flags, wildcard_ref, wildcard_query)
    aligner.min_overlap = min_overlap
    return aligner.locate(query)
