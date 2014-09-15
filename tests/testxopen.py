from __future__ import print_function, division, absolute_import

from cutadapt.xopen import xopen

base = "tests/data/small.fastq"
files = [ base + ext for ext in ['', '.gz', '.bz2' ] ]

def test_xopen():
	for name in files:
		f = xopen(name, 'r')
		lines = list(f)
		assert len(lines) == 12
		assert lines[5] == 'AGCCGCTANGACGGGTTGGCCCTTAGACGTATCT\n', name
		f.close()
