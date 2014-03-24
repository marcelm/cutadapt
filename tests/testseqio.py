from __future__ import print_function, division, absolute_import

from nose.tools import raises
from cutadapt import seqio
import sys

from io import BytesIO

# files tests/data/simple.fast{q,a}
simple_fastq = [
	seqio.Sequence("first_sequence", b"SEQUENCE1", b":6;;8<=:<"),
	seqio.Sequence("second_sequence", b"SEQUENCE2", b"83<??:(61")
	]

simple_fasta = [ seqio.Sequence(x.name, x.sequence, None) for x in simple_fastq ]


def test_fastareader():
	with seqio.FastaReader("tests/data/simple.fasta") as f:
		reads = list(f)
	assert reads == simple_fasta


def test_fastqreader():
	with seqio.FastqReader("tests/data/simple.fastq") as f:
		reads = list(f)
	assert reads == simple_fastq


@raises(seqio.FormatError)
def test_fastq_wrongformat():
	with seqio.FastqReader("tests/data/withplus.fastq") as f:
		reads = list(f)


def test_sequence_reader():
	# test the autodetection
	with seqio.SequenceReader("tests/data/simple.fastq") as f:
		reads = list(f)
	assert reads == simple_fastq

	with seqio.SequenceReader("tests/data/simple.fasta") as f:
		reads = list(f)
	assert reads == simple_fasta

	with open("tests/data/simple.fastq", 'rb') as f:
		reads = list(seqio.SequenceReader(f))
	assert reads == simple_fastq

	# make the name attribute unavailable
	f = BytesIO(open("tests/data/simple.fastq", 'rb').read())
	reads = list(seqio.SequenceReader(f))
	assert reads == simple_fastq

	f = BytesIO(open("tests/data/simple.fasta", 'rb').read())
	reads = list(seqio.SequenceReader(f))
	assert reads == simple_fasta


def test_context_manager():
	filename = "tests/data/simple.fasta"
	with open(filename, 'rb') as f:
		assert not f.closed
		reads = list(seqio.SequenceReader(f))
		assert not f.closed
	assert f.closed

	with seqio.FastaReader(filename) as sr:
		tmp_sr = sr
		assert not sr.fp.closed
		reads = list(sr)
		assert not sr.fp.closed
	assert tmp_sr.fp.closed
