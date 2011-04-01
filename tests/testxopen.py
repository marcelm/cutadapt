from cutadapt.xopen import xopen
from cutadapt import seqio
import sys

uncompressed = "tests/data/simple.fasta"
compressed = "tests/data/compressed.fastq.gz"

if sys.version_info[0] < 3:
	def test_xopen():
		for name in [compressed, uncompressed]:
			f = xopen(name)
			assert not f.closed
			f.close()
			assert f.closed

			with xopen(name) as f:
				assert not f.closed
				data = f.read()
			assert f.closed
