from __future__ import print_function, division, absolute_import
from nose.tools import raises
from cutadapt.xopen import xopen
import gzip
import os
import random
from utils import temporary_path

base = "tests/data/small.fastq"
files = [ base + ext for ext in ['', '.gz', '.bz2' ] ]

def test_xopen():
	for name in files:
		f = xopen(name, 'r')
		lines = list(f)
		assert len(lines) == 12
		assert lines[5] == 'AGCCGCTANGACGGGTTGGCCCTTAGACGTATCT\n', name
		f.close()


def create_truncated_file(path):
	# Random text
	text = ''.join(random.choice('ABCDEFGHIJKLMNOPQRSTUVWXYZ') for _ in range(200))
	f = xopen(path, 'w')
	f.write(text)
	f.close()
	f = open(path, 'a')
	f.truncate(os.stat(path).st_size - 10)
	f.close()


@raises(EOFError)
def test_truncated_gz():
	with temporary_path('truncated.gz') as path:
		create_truncated_file(path)

		f = xopen(path, 'r')
		f.read()
		f.close()


@raises(EOFError)
def test_truncated_gz_iter():
	with temporary_path('truncated.gz') as path:
		create_truncated_file(path)

		f = xopen(path, 'r')
		for line in f:
			pass
		f.close()
