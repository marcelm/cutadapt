# coding: utf-8
from __future__ import print_function, division, absolute_import

import sys
import os
import shutil
from textwrap import dedent
from nose.tools import raises
from tempfile import mkdtemp
from cutadapt.seqio import (Sequence, ColorspaceSequence, FormatError,
	FastaReader, FastqReader, FastaQualReader, InterleavedSequenceReader,
	FastaWriter, FastqWriter, InterleavedSequenceWriter, open as openseq,
	sequence_names_match)
from cutadapt.compat import StringIO


# files tests/data/simple.fast{q,a}
simple_fastq = [
	Sequence("first_sequence", "SEQUENCE1", ":6;;8<=:<"),
	Sequence("second_sequence", "SEQUENCE2", "83<??:(61")
	]

simple_fasta = [ Sequence(x.name, x.sequence, None) for x in simple_fastq ]


class TestSequence:
	@raises(FormatError)
	def test_too_many_qualities(self):
		Sequence(name="name", sequence="ACGT", qualities="#####")

	@raises(FormatError)
	def test_too_many_qualities_colorspace(self):
		ColorspaceSequence(name="name", sequence="T0123", qualities="#####")

	@raises(FormatError)
	def test_invalid_primer(self):
		ColorspaceSequence(name="name", sequence="K0123", qualities="####")


class TestFastaReader:
	def test(self):
		with FastaReader("tests/data/simple.fasta") as f:
			reads = list(f)
		assert reads == simple_fasta

		fasta = StringIO(">first_sequence\nSEQUENCE1\n>second_sequence\nSEQUENCE2\n")
		reads = list(FastaReader(fasta))
		assert reads == simple_fasta

	def test_with_comments(self):
		fasta = StringIO(dedent(
			"""
			# a comment
			# another one
			>first_sequence
			SEQUENCE1
			>second_sequence
			SEQUENCE2
			"""))
		reads = list(FastaReader(fasta))
		assert reads == simple_fasta

	@raises(FormatError)
	def test_wrong_format(self):
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
		reads = list(FastaReader(fasta))

	def test_fastareader_keeplinebreaks(self):
		with FastaReader("tests/data/simple.fasta", keep_linebreaks=True) as f:
			reads = list(f)
		assert reads[0] == simple_fasta[0]
		assert reads[1].sequence == 'SEQUEN\nCE2'

	def test_context_manager(self):
		filename = "tests/data/simple.fasta"
		with open(filename) as f:
			assert not f.closed
			reads = list(openseq(f))
			assert not f.closed
		assert f.closed

		with FastaReader(filename) as sr:
			tmp_sr = sr
			assert not sr._file.closed
			reads = list(sr)
			assert not sr._file.closed
		assert tmp_sr._file is None
		# Open it a second time
		with FastaReader(filename) as sr:
			pass


class TestFastqReader:
	def test_fastqreader(self):
		with FastqReader("tests/data/simple.fastq") as f:
			reads = list(f)
		assert reads == simple_fastq

	def test_fastqreader_dos(self):
		with FastqReader("tests/data/dos.fastq") as f:
			dos_reads = list(f)
		with FastqReader("tests/data/small.fastq") as f:
			unix_reads = list(f)
		assert dos_reads == unix_reads

	@raises(FormatError)
	def test_fastq_wrongformat(self):
		with FastqReader("tests/data/withplus.fastq") as f:
			reads = list(f)

	@raises(FormatError)
	def test_fastq_incomplete(self):
		fastq = StringIO("@name\nACGT+\n")
		with FastqReader(fastq) as fq:
			list(fq)

	def test_context_manager(self):
		filename = "tests/data/simple.fastq"
		with open(filename) as f:
			assert not f.closed
			reads = list(openseq(f))
			assert not f.closed
		assert f.closed

		with FastqReader(filename) as sr:
			tmp_sr = sr
			assert not sr._file.closed
			reads = list(sr)
			assert not sr._file.closed
		assert tmp_sr._file is None


class TestFastaQualReader:
	@raises(FormatError)
	def test_mismatching_read_names(self):
		fasta = StringIO(">name\nACG")
		qual = StringIO(">nome\n3 5 7")
		list(FastaQualReader(fasta, qual))

	@raises(FormatError)
	def test_invalid_quality_value(self):
		fasta = StringIO(">name\nACG")
		qual = StringIO(">name\n3 xx 7")
		list(FastaQualReader(fasta, qual))


class TestSeqioOpen:
	def setup(self):
		self._tmpdir = mkdtemp()

	def teardown(self):
		shutil.rmtree(self._tmpdir)

	def test_sequence_reader(self):
		# test the autodetection
		with openseq("tests/data/simple.fastq") as f:
			reads = list(f)
		assert reads == simple_fastq

		with openseq("tests/data/simple.fasta") as f:
			reads = list(f)
		assert reads == simple_fasta

		with open("tests/data/simple.fastq") as f:
			reads = list(openseq(f))
		assert reads == simple_fastq

		# make the name attribute unavailable
		f = StringIO(open("tests/data/simple.fastq").read())
		reads = list(openseq(f))
		assert reads == simple_fastq

		f = StringIO(open("tests/data/simple.fasta").read())
		reads = list(openseq(f))
		assert reads == simple_fasta

	def test_autodetect_fasta_format(self):
		path = os.path.join(self._tmpdir, 'tmp.fasta')
		with openseq(path, mode='w') as f:
			assert isinstance(f, FastaWriter)
			for seq in simple_fastq:
				f.write(seq)
		assert list(openseq(path)) == simple_fasta

	def test_write_qualities_to_fasta(self):
		path = os.path.join(self._tmpdir, 'tmp.fasta')
		with openseq(path, mode='w', qualities=True) as f:
			assert isinstance(f, FastaWriter)
			for seq in simple_fastq:
				f.write(seq)
		assert list(openseq(path)) == simple_fasta

	def test_autodetect_fastq_format(self):
		path = os.path.join(self._tmpdir, 'tmp.fastq')
		with openseq(path, mode='w') as f:
			assert isinstance(f, FastqWriter)
			for seq in simple_fastq:
				f.write(seq)
		assert list(openseq(path)) == simple_fastq

	@raises(ValueError)
	def test_fastq_qualities_missing(self):
		path = os.path.join(self._tmpdir, 'tmp.fastq')
		openseq(path, mode='w', qualities=False)


class TestInterleavedReader:
	def test(self):
		expected = [
			(Sequence('read1/1 some text', 'TTATTTGTCTCCAGC', '##HHHHHHHHHHHHH'),
			Sequence('read1/2 other text', 'GCTGGAGACAAATAA', 'HHHHHHHHHHHHHHH')),
			(Sequence('read3/1', 'CCAACTTGATATTAATAACA', 'HHHHHHHHHHHHHHHHHHHH'),
			Sequence('read3/2', 'TGTTATTAATATCAAGTTGG', '#HHHHHHHHHHHHHHHHHHH'))
		]
		reads = list(InterleavedSequenceReader("tests/cut/interleaved.fastq"))
		for (r1, r2), (e1, e2) in zip(reads, expected):
			print(r1, r2, e1, e2)

		assert reads == expected
		with openseq("tests/cut/interleaved.fastq", interleaved=True) as f:
			reads = list(f)
		assert reads == expected

	@raises(FormatError)
	def test_missing_partner(self):
		s = StringIO('@r1\nACG\n+\nHHH')
		list(InterleavedSequenceReader(s))

	@raises(FormatError)
	def test_incorrectly_paired(self):
		s = StringIO('@r1/1\nACG\n+\nHHH\n@wrong_name\nTTT\n+\nHHH')
		list(InterleavedSequenceReader(s))


class TestFastaWriter:
	def setup(self):
		self._tmpdir = mkdtemp()
		self.path = os.path.join(self._tmpdir, 'tmp.fasta')

	def teardown(self):
		shutil.rmtree(self._tmpdir)

	def test(self):
		with FastaWriter(self.path) as fw:
			fw.write("name", "CCATA")
			fw.write("name2", "HELLO")
		assert fw._file.closed
		with open(self.path) as t:
			assert t.read() == '>name\nCCATA\n>name2\nHELLO\n'

	def test_linelength(self):
		with FastaWriter(self.path, line_length=3) as fw:
			fw.write("r1", "ACG")
			fw.write("r2", "CCAT")
			fw.write("r3", "TACCAG")
		assert fw._file.closed
		with open(self.path) as t:
			d = t.read()
			assert d == '>r1\nACG\n>r2\nCCA\nT\n>r3\nTAC\nCAG\n'

	def test_write_sequence_object(self):
		with FastaWriter(self.path) as fw:
			fw.write(Sequence("name", "CCATA"))
			fw.write(Sequence("name2", "HELLO"))
		assert fw._file.closed
		with open(self.path) as t:
			assert t.read() == '>name\nCCATA\n>name2\nHELLO\n'

	def test_write_to_file_like_object(self):
		sio = StringIO()
		with FastaWriter(sio) as fw:
			fw.write(Sequence("name", "CCATA"))
			fw.write(Sequence("name2", "HELLO"))
			assert sio.getvalue() == '>name\nCCATA\n>name2\nHELLO\n'
		assert not fw._file.closed

	def test_write_zero_length_sequence(self):
		sio = StringIO()
		with FastaWriter(sio) as fw:
			fw.write(Sequence("name", ""))
			assert sio.getvalue() == '>name\n\n', '{0!r}'.format(sio.getvalue())


class TestFastqWriter:
	def setup(self):
		self._tmpdir = mkdtemp()
		self.path = os.path.join(self._tmpdir, 'tmp.fastq')

	def teardown(self):
		shutil.rmtree(self._tmpdir)

	def test(self):
		with FastqWriter(self.path) as fq:
			fq.writeseq("name", "CCATA", "!#!#!")
			fq.writeseq("name2", "HELLO", "&&&!&&")
		assert fq._file.closed
		with open(self.path) as t:
			assert t.read() == '@name\nCCATA\n+\n!#!#!\n@name2\nHELLO\n+\n&&&!&&\n'

	def test_twoheaders(self):
		with FastqWriter(self.path) as fq:
			fq.write(Sequence("name", "CCATA", "!#!#!", second_header=True))
			fq.write(Sequence("name2", "HELLO", "&&&!&", second_header=True))
		assert fq._file.closed
		with open(self.path) as t:
			assert t.read() == '@name\nCCATA\n+name\n!#!#!\n@name2\nHELLO\n+name2\n&&&!&\n'

	def test_write_to_file_like_object(self):
		sio = StringIO()
		with FastqWriter(sio) as fq:
			fq.writeseq("name", "CCATA", "!#!#!")
			fq.writeseq("name2", "HELLO", "&&&!&&")
		assert sio.getvalue() == '@name\nCCATA\n+\n!#!#!\n@name2\nHELLO\n+\n&&&!&&\n'


class TestInterleavedWriter:
	def test(self):
		reads = [
			(Sequence('A/1 comment', 'TTA', '##H'),
			Sequence('A/2 comment', 'GCT', 'HH#')),
			(Sequence('B/1', 'CC', 'HH'),
			Sequence('B/2', 'TG', '#H'))
		]
		sio = StringIO()
		with InterleavedSequenceWriter(sio) as writer:
			for read1, read2 in reads:
				writer.write(read1, read2)
		assert sio.getvalue() == '@A/1 comment\nTTA\n+\n##H\n@A/2 comment\nGCT\n+\nHH#\n@B/1\nCC\n+\nHH\n@B/2\nTG\n+\n#H\n'


class TestPairedSequenceReader:
	def test_sequence_names_match(self):
		def match(name1, name2):
			seq1 = Sequence(name1, 'ACGT')
			seq2 = Sequence(name2, 'AACC')
			return sequence_names_match(seq1, seq2)

		assert match('abc', 'abc')
		assert match('abc/1', 'abc/2')
		assert match('abc.1', 'abc.2')
		assert match('abc1', 'abc2')
		assert not match('abc', 'xyz')

