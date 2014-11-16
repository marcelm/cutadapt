from __future__ import print_function, division, absolute_import

from nose.tools import raises
from cutadapt.scripts import cutadapt
from utils import run, files_equal, datapath, cutpath, redirect_stderr, temporary_path

def test_paired_separate():
	'''test separate trimming of paired-end reads'''
	run('-a TTAGACATAT', 'paired-separate.1.fastq', 'paired.1.fastq')
	run('-a CAGTGGAGTA', 'paired-separate.2.fastq', 'paired.2.fastq')


def test_paired_end_legacy():
	'''--paired-output, no -A/-B/-G'''
	with temporary_path("paired-tmp.fastq") as pairedtmp:
		# the -m 14 filters out one read, which should then also be filtered out in the second output file
		run(['-a', 'TTAGACATAT', '-m', '14', '--paired-output', pairedtmp], 'paired.m14.1.fastq', 'paired.1.fastq', 'paired.2.fastq')
		assert files_equal(cutpath('paired.m14.2.fastq'), pairedtmp)


def test_paired_end():
	'''-p, -m and -A'''
	with temporary_path("paired.1.fastq") as p1:
		with temporary_path("paired.2.fastq") as p2:
			params = [
				'-a', 'TTAGACATAT',
				'-A', 'CAGTGGAGTA',
				'-m', '14',
				'-o', p1, '-p', p2,
				datapath('paired.1.fastq'), datapath('paired.2.fastq')
			]
			assert cutadapt.main(params) is None
			assert files_equal(cutpath('paired.1.fastq'), p1)
			assert files_equal(cutpath('paired.2.fastq'), p2)


def test_untrimmed_paired_output():
	with temporary_path("tmp-paired.1.fastq") as tmp1:
		with temporary_path("tmp-paired.2.fastq") as tmp2:
			with temporary_path("tmp-untrimmed.1.fastq") as untrimmed1:
				with temporary_path("tmp-untrimmed.2.fastq") as untrimmed2:
					params = [
						'-a', 'TTAGACATAT',
						'-o', tmp1, '-p', tmp2,
						'--untrimmed-output', untrimmed1,
						'--untrimmed-paired-output', untrimmed2,
						datapath('paired.1.fastq'), datapath('paired.2.fastq')
					]
					assert cutadapt.main(params) is None
					assert files_equal(cutpath('paired-untrimmed.1.fastq'), untrimmed1)
					assert files_equal(cutpath('paired-untrimmed.2.fastq'), untrimmed2)
					assert files_equal(cutpath('paired-trimmed.1.fastq'), tmp1)
					assert files_equal(cutpath('paired-trimmed.2.fastq'), tmp2)


def test_explicit_format_with_paired():
	with temporary_path("paired-tmp.fastq") as pairedtmp:
		run(['--format=fastq', '-a', 'TTAGACATAT', '-m', '14', '-p', pairedtmp], 'paired.m14.1.fastq', 'paired.1.txt', 'paired.2.txt')
		assert files_equal(cutpath('paired.m14.2.fastq'), pairedtmp)


def test_no_trimming():
	# make sure that this doesn't divide by zero
	cutadapt.main(['-a', 'XXXXX', '-o', '/dev/null', '-p', '/dev/null', datapath('paired.1.fastq'), datapath('paired.2.fastq')])


@raises(SystemExit)
def test_paired_end_missing_file():
	with redirect_stderr():
		cutadapt.main(['-a', 'XX', '--paired-output', 'out.fastq', datapath('paired.1.fastq')])


@raises(SystemExit)
def test_first_too_short():
	# paired-truncated.1.fastq is paired.1.fastq without the last read
	with redirect_stderr():
		cutadapt.main('-a XX --paired-output out.fastq'.split() + [datapath('paired-truncated.1.fastq'), datapath('paired.2.fastq')])


@raises(SystemExit)
def test_second_too_short():
	# paired-truncated.2.fastq is paired.2.fastq without the last read
	with redirect_stderr():
		cutadapt.main('-a XX --paired-output out.fastq'.split() + [datapath('paired.1.fastq'), datapath('paired-truncated.2.fastq')])


@raises(SystemExit)
def test_unmatched_read_names():
	# paired-swapped.1.fastq: paired.1.fastq with reads 2 and 3 swapped
	with redirect_stderr():
		cutadapt.main('-a XX --paired-output out.fastq'.split() + [datapath('paired-swapped.1.fastq'), datapath('paired.2.fastq')])
