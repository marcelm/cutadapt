from __future__ import print_function, division

from cutadapt.seqio import ColorspaceSequence
from cutadapt.adapters import ColorspaceAdapter, PREFIX
from cutadapt.scripts.cutadapt import RepeatedAdapterCutter

def test_cs_5p():
	read = ColorspaceSequence("name", "0123", "DEFG", "T")
	adapter = ColorspaceAdapter("CG", PREFIX, 0.1)
	cutter = RepeatedAdapterCutter([adapter])
	matches = cutter.find_match(read)
	# no assertion here, just make sure the above code runs without
	# an exception
