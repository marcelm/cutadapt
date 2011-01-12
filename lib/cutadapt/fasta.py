#!/usr/bin/env python
"""
This module provides routines for reading and writing FASTA and FASTQ files.
"""

def readfastq(infile, colorspace=False):
	"""
	Reads a FASTQ file and returns a generator of tuples: (description, sequence, qualities)

	infile -- file-like object (for example, a file or a GzipFile)
	colorspace -- Usually (when this is False), there must be n characters in the sequence and
	n quality values. When this is True, there must be n+1 characters in the sequence and n quality values.
	"""
	lengthdiff = 1 if colorspace else 0
	for i, line in enumerate(infile):
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
				raise ValueError("Length of quality sequence and length of read do not match (%d+%d!=%d)" % (len(qualities), lengthdiff, len(sequence)))
			yield (description, sequence, qualities)


def readfasta(infile):
	"""
	Read a FASTA file and return an iterator over tuples: (comment, sequence).

	f -- file-like object (for example, a file or a GzipFile)

	NOTE: This is (in theory) quadratic for many FASTA files since += is used for the
	sequence string. This is not a problem in practice.
	"""
	comment = None
	seq = ""

	for line in infile:
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


def writefastq(outfile, sequences):
	"""
	Write sequences to the file outfile in FASTQ format.
	outfile -- a file-like object that has a 'write' method
	seqlist -- a list of (description, sequence, qualities) tuples
	"""
	for description, sequence, qualities in sequences:
		outfile.write('@%s\n%s\n+\n%s\n' % (description, sequence, qualities))
