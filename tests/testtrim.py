from __future__ import print_function, division

from cutadapt.seqio import ColorspaceSequence
from cutadapt.adapters import ColorspaceAdapter, PREFIX
from cutadapt.scripts.cutadapt import AdapterCutter

def test_cs_5p():
	read = ColorspaceSequence("name", "0123", "DEFG", "T")
	adapter = ColorspaceAdapter("CG", PREFIX, 0.1)
	cutter = AdapterCutter([adapter])
	matches = cutter.find_matches(read)
	# no assertion here, just make sure the above code runs without
	# an exception
