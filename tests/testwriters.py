# coding: utf-8
import os
from cutadapt import writers
from cutadapt.seqio import Sequence
from .utils import temporary_directory, temporary_path, files_equal, cutpath, datapath

def test_single_end():
	w = writers.Writers()
	with temporary_path("result.fq") as path:
		w.add_writer(None, path)
		assert w.has_writer(None)
		assert isinstance(w.get_writer(None), writers.SingleEndWriter)
		read = Sequence(
			name="prefix:1_13_1440/1", 
			sequence="CAAGATCTNCCCTGCCACATTGCCCTAGTTAAAC",
			qualities="<=A:A=57!7<';<6?5;;6:+:=)71>70<,=:"
		)
		w.write(None, read)
		w.close()
		assert files_equal(cutpath('no-trim.fastq'), path)
		summary = w.summary()
		assert summary[0] == 1
		assert summary[1][0] == 34
		assert summary[1][1] == 0

def test_paired_end():
	w = writers.Writers()
	with temporary_path("result1.fq") as path1:
		with temporary_path("result2.fq") as path2:
			w.add_writer(None, path1, path2)
			assert w.has_writer(None)
			assert isinstance(w.get_writer(None), writers.PairedEndWriter)
			read1 = Sequence(
				name="read1/1 some text", 
				sequence="TTATTTGTCTCCAGCTTAGACATATCGCCT",
				qualities="##HHHHHHHHHHHHHHHHHHHHHHHHHHHH"
			)
			read2 = Sequence(
				name="read1/2 other text", 
				sequence="GCTGGAGACAAATAACAGTGGAGTAGTTTT",
				qualities="HHHHHHHHHHHHHHHHHHHHHHHHHHHHHH"
			)
			w.write(None, read1, read2)
			w.close()
			assert files_equal(datapath('small.paired.1.fastq'), path1)
			assert files_equal(datapath('small.paired.2.fastq'), path2)
			summary = w.summary()
			assert summary[0] == 1
			assert summary[1][0] == 30
			assert summary[1][1] == 30

def test_multiplexed():
	w = writers.Writers(multiplexed=True, name_pattern=os.path.join(temporary_directory(), "foo.{name}.fa"))
	with temporary_path("foo.alpha.fa") as path1:
		with temporary_path("foo.beta.fa") as path2:
			assert isinstance(w.get_multiplex_writer("alpha"), writers.SingleEndWriter)
			assert isinstance(w.get_multiplex_writer("beta"), writers.SingleEndWriter)
		
			class DummyAdapter(object):
				def __init__(self, name):
					self.name = name
			class DummyMatch(object):
				def __init__(self, name):
					self.adapter = DummyAdapter(name)

			read1 = Sequence(
				name="read1",
				sequence="GATCCTCCTGGAGCTGGCTGATACCAGTATACCAGTGCTGATTGTTG",
				match = DummyMatch("alpha")
			)
			w.write(None, read1)
			read2 = Sequence(
				name="read2",
				sequence="CTCGAGAATTCTGGATCCTCTCTTCTGCTACCTTTGGGATTTGCTTGCTCTTG",
				match=DummyMatch("beta")
			)
			w.write(None, read2)
			w.close()
			assert files_equal(cutpath('twoadapters.first.fasta'), path1)
			assert files_equal(cutpath('twoadapters.second.fasta'), path2)
			summary = w.summary()
			assert summary[0] == 2
			assert summary[1][0] == 100
			assert summary[1][1] == 0

def test_force_create():
	w = writers.Writers()
	with temporary_path("empty.fastq") as path:
		w.add_writer(None, path, force_create=True)
		w.close()
		assert files_equal(datapath("empty.fastq"), path)
