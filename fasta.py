#!/usr/bin/env python
import gzip
"""
This module provides routines for reading and writing FASTA and FASTQ files.

TODO
  * some functions accept files, some accept file names
  * make all functions work with gzipped data
"""

def readfastq(f, colorspace=False):
	"""
	Reads a FASTQ file and returns a generator of tuples: (description, sequence, qualities)

	colorspace -- Usually (when this is False), there must be n characters in the sequence and
	n quality values. When this is True, there must be n+1 characters in the sequence and n quality values.
	"""
	assert type(f) is not str

	lengthdiff = 1 if colorspace else 0
	for i, line in enumerate(f):
		if i % 4 == 0:
			assert line[0] == '@'
			description = line.strip()
			description = description[1:]
		elif i % 4 ==  1:
			sequence = line.strip()
		elif i % 4 == 2:
			assert line.startswith('+')
		elif i % 4 == 3:
			qualities = line[:-1]
			if len(qualities) + lengthdiff != len(sequence):
				raise ValueError, "Length of quality sequence and length of read do not match (%d+%d!=%d)" % (len(qualities), lengthdiff, len(sequence))
			yield (description, sequence, qualities)

def readfasta(f):
	"""
	Read a FASTA file and return an iterator over tuples: (comment, sequence).

	f is a filename or a file-like object. If f is a filename, then .gz files are supported.

	NOTE: This is (in theory) quadratic for many FASTA files since += is used for the
	sequence string. This is not a problem in practice.
	"""
	close = False
	if type(f) is str:
		if f.endswith('.gz'):
			f = gzip.open(f)
			close = True
		else:
			f = open(f)
			close = True
	comment = None
	seq = ""

	for line in f:
		# stripping also deals with DOS line break issues
		line = line.strip()
		if line and line[0] == ">":
			if comment is not None:
				assert seq.find('\n') == -1
				yield (comment, seq)
			comment = line[1:]
			seq = ""
		else:
			seq += line
	if comment is not None:
		assert seq.find('\n') == -1
		yield (comment, seq)
	# TODO better use a context manager
	if close:
		f.close()

def writefastq(f, sequences):
	"""
	Write sequences to the file f.
	seqlist -- a list of (description, sequence, qualities) tuples
	"""
	for description, sequence, qualities in sequences:
		f.write('@%s\n%s\n+\n%s\n' % (description, sequence, qualities))

class UnknownFileType(Exception):
	pass

def fastafiletype(fname):
	"""
	Determine file type of fname. Return the string FASTQ or FASTA or
	raise an UnknownFileType exception.
	"""
	with open(fname) as f:
		for line in f:
			if line.startswith('#'):
				continue
			if line.startswith('@'):
				return 'FASTQ'
			if line.startswith('>'):
				return 'FASTA'
			raise UnknownFileType("neither FASTQ nor FASTA")

def readfastaq(name, colorspace=False):
	"""
	Read a FASTA or FASTQ file. The file type is recognized automatically.
	Return tuples (description, sequence, qualities).
	All elements of the tuple are strings (the qualities are not decoded).
	If a FASTA file was detected, qualities is always None.
	"""
	if fastafiletype(name) == 'FASTQ':
		with open(name) as fastqfile:
			for desc, seq, qual in readfastq(fastqfile, colorspace):
				yield desc, seq, qual
	else:
		with open(name) as fastafile:
			for desc, seq in readfasta(fastafile):
				yield desc, seq, None
