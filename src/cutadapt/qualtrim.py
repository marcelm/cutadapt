"""
Quality trimming.
"""
import sys

if sys.version > '3':
	xrange = range


from cutadapt._qualtrim import quality_trim_index, nextseq_trim_index
