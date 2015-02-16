# coding: utf-8
from __future__ import print_function, division, absolute_import

from nose.tools import raises
from cutadapt import seqio
import sys
from cutadapt.compat import StringIO


# files tests/data/simple.fast{q,a}
simple_fastq = [
	seqio.Sequence("first_sequence", "SEQUENCE1", ":6;;8<=:<"),
	seqio.Sequence("second_sequence", "SEQUENCE2", "83<??:(61")
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


def test_fastareader_keeplinebreaks():
	with seqio.FastaReader("tests/data/simple.fasta", keep_linebreaks=True) as f:
		reads = list(f)
	assert reads[0] == simple_fasta[0]
	assert reads[1].sequence == 'SEQUEN\nCE2'


@raises(seqio.FormatError)
def test_fastq_wrongformat():
	with seqio.FastqReader("tests/data/withplus.fastq") as f:
		reads = list(f)


@raises(seqio.FormatError)
def test_too_many_qualities():
	seqio.Sequence(name="name", sequence="ACGT", qualities="#####")


@raises(seqio.FormatError)
def test_too_many_qualities_colorspace():
	seqio.ColorspaceSequence(name="name", sequence="T0123", qualities="#####")


@raises(seqio.FormatError)
def test_invalid_primer():
	seqio.ColorspaceSequence(name="name", sequence="K0123", qualities="####")


@raises(seqio.FormatError)
def test_mismatching_read_names():
	fasta = StringIO(">name\nACG")
	qual = StringIO(">nome\n3 5 7")
	list(seqio.FastaQualReader(fasta, qual))


@raises(seqio.FormatError)
def test_invalid_quality_value():
	fasta = StringIO(">name\nACG")
	qual = StringIO(">name\n3 xx 7")
	list(seqio.FastaQualReader(fasta, qual))


def test_sequence_reader():
	# test the autodetection
	with seqio.SequenceReader("tests/data/simple.fastq") as f:
		reads = list(f)
	assert reads == simple_fastq

	with seqio.SequenceReader("tests/data/simple.fasta") as f:
		reads = list(f)
	assert reads == simple_fasta

	with open("tests/data/simple.fastq") as f:
		reads = list(seqio.SequenceReader(f))
	assert reads == simple_fastq

	# make the name attribute unavailable
	f = StringIO(open("tests/data/simple.fastq").read())
	reads = list(seqio.SequenceReader(f))
	assert reads == simple_fastq

	f = StringIO(open("tests/data/simple.fasta").read())
	reads = list(seqio.SequenceReader(f))
	assert reads == simple_fasta


def test_context_manager():
	filename = "tests/data/simple.fasta"
	with open(filename) as f:
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
