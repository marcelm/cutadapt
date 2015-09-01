# coding: utf-8
from __future__ import print_function, division, absolute_import

import sys
from textwrap import dedent
from nose.tools import raises
from cutadapt import seqio
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

	fasta = StringIO(">first_sequence\nSEQUENCE1\n>second_sequence\nSEQUENCE2\n")
	reads = list(seqio.FastaReader(fasta))
	assert reads == simple_fasta


def test_fastareader_with_comments():
	fasta = StringIO(dedent(
		"""
		# a comment
		# another one
		>first_sequence
		SEQUENCE1
		>second_sequence
		SEQUENCE2
		"""))
	reads = list(seqio.FastaReader(fasta))
	assert reads == simple_fasta


@raises(seqio.FormatError)
def test_wrong_fasta_format():
	fasta = StringIO(dedent(
		"""
		# a comment
		# another one
		unexpected
		>first_sequence
		SEQUENCE1
		>second_sequence
		SEQUENCE2
		"""))
	reads = list(seqio.FastaReader(fasta))


def test_fastqreader():
	with seqio.FastqReader("tests/data/simple.fastq") as f:
		reads = list(f)
	assert reads == simple_fastq


def test_fastqreader_dos():
	with seqio.FastqReader("tests/data/dos.fastq") as f:
		dos_reads = list(f)
	with seqio.FastqReader("tests/data/small.fastq") as f:
		unix_reads = list(f)
	assert dos_reads == unix_reads


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
def test_fastq_incomplete():
	fastq = StringIO("@name\nACGT+\n")
	with seqio.FastqReader(fastq) as fq:
		list(fq)


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
	with seqio.open("tests/data/simple.fastq") as f:
		reads = list(f)
	assert reads == simple_fastq

	with seqio.open("tests/data/simple.fasta") as f:
		reads = list(f)
	assert reads == simple_fasta

	with open("tests/data/simple.fastq") as f:
		reads = list(seqio.open(f))
	assert reads == simple_fastq

	# make the name attribute unavailable
	f = StringIO(open("tests/data/simple.fastq").read())
	reads = list(seqio.open(f))
	assert reads == simple_fastq

	f = StringIO(open("tests/data/simple.fasta").read())
	reads = list(seqio.open(f))
	assert reads == simple_fasta


def test_fasta_context_manager():
	filename = "tests/data/simple.fasta"
	with open(filename) as f:
		assert not f.closed
		reads = list(seqio.open(f))
		assert not f.closed
	assert f.closed

	with seqio.FastaReader(filename) as sr:
		tmp_sr = sr
		assert not sr.fp.closed
		reads = list(sr)
		assert not sr.fp.closed
	assert tmp_sr.fp is None


def test_fastq_context_manager():
	filename = "tests/data/simple.fastq"
	with open(filename) as f:
		assert not f.closed
		reads = list(seqio.open(f))
		assert not f.closed
	assert f.closed

	with seqio.FastqReader(filename) as sr:
		tmp_sr = sr
		assert not sr.fp.closed
		reads = list(sr)
		assert not sr.fp.closed
	assert tmp_sr.fp is None


def test_interleaved_reader():
	expected = [
		(seqio.Sequence('read1/1 some text', 'TTATTTGTCTCCAGC', '##HHHHHHHHHHHHH'),
		seqio.Sequence('read1/2 other text', 'GCTGGAGACAAATAA', 'HHHHHHHHHHHHHHH')),
		(seqio.Sequence('read3/1', 'CCAACTTGATATTAATAACA', 'HHHHHHHHHHHHHHHHHHHH'),
		seqio.Sequence('read3/2', 'TGTTATTAATATCAAGTTGG', '#HHHHHHHHHHHHHHHHHHH'))
	]
	reads = list(seqio.InterleavedSequenceReader("tests/cut/interleaved.fastq"))
	for (r1, r2), (e1, e2) in zip(reads, expected):
		print(r1, r2, e1, e2)

	assert reads == expected
	with seqio.open("tests/cut/interleaved.fastq", interleaved=True) as f:
		reads = list(f)
	assert reads == expected
