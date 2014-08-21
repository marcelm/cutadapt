"""
Alignment module.
"""
from __future__ import print_function, division, absolute_import

from array import array
import sys

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
ALLOW_WILDCARD_SEQ1 = 1
ALLOW_WILDCARD_SEQ2 = 2

GLOBAL = 0

# Use this to get regular semiglobal alignment
# (all gaps in the beginning or end are free)
SEMIGLOBAL = START_WITHIN_SEQ1 | START_WITHIN_SEQ2 | STOP_WITHIN_SEQ1 | STOP_WITHIN_SEQ2


from cutadapt._align import Aligner, compare_prefixes, locate
