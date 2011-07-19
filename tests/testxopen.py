from cutadapt.xopen import xopen
from cutadapt import seqio
import sys

uncompressed = "tests/data/simple.fasta"
compressed = "tests/data/compressed.fastq.gz"

def test_xopen():
	f = xopen(uncompressed)
	assert f.readline() == '# a comment\n'
	f.close()

	f = xopen(compressed)
	assert f.readline() == '@first_sequence\n'
	f.close()
