# coding: utf-8
"""
Quality trimming.
"""
from __future__ import print_function, division, absolute_import

import sys

if sys.version > '3':
	xrange = range


from cutadapt._qualtrim import quality_trim_index, nextseq_trim_index
