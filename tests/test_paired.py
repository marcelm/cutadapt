# coding: utf-8
from __future__ import print_function, division, absolute_import

import os.path
import shutil
import tempfile

import pytest
from nose.tools import raises

from cutadapt.compat import PY3
from cutadapt.__main__ import main
from utils import run, assert_files_equal, datapath, cutpath, redirect_stderr, temporary_path


if PY3:
	@pytest.fixture(params=[1, 2])
	def cores(request):
		return request.param
else:
	@pytest.fixture
	def cores():
		return 1


def run_paired(params, in1, in2, expected1, expected2, cores):
	if type(params) is str:
		params = params.split()
	params += ['--cores', str(cores), '--buffer-size=512']
	with temporary_path('tmp1-' + expected1) as p1:
		with temporary_path('tmp2-' + expected2) as p2:
			params += ['-o', p1, '-p', p2]
			params += [datapath(in1), datapath(in2)]
			assert main(params) is None
			assert_files_equal(cutpath(expected1), p1)
			assert_files_equal(cutpath(expected2), p2)


def run_interleaved(params, inpath1, inpath2=None, expected1=None, expected2=None, cores=1):
	"""
	Interleaved input or output (or both)
	"""
	assert not (inpath1 and inpath2 and expected1 and expected2)
	assert not (expected2 and not expected1)
	assert not (inpath2 and not inpath1)
	if type(params) is str:
		params = params.split()
	params += ['--interleaved', '--cores', str(cores), '--buffer-size=512']
	with temporary_path('tmp1-' + expected1) as tmp1:
		params += ['-o', tmp1]
		paths = [datapath(inpath1)]
		if inpath2:
			paths += [datapath(inpath2)]
		if expected2:
			with temporary_path('tmp2-' + expected2) as tmp2:
				params += ['-p', tmp2]
				assert main(params + paths) is None
				assert_files_equal(cutpath(expected2), tmp2)
		else:
			assert main(params + paths) is None
		assert_files_equal(cutpath(expected1), tmp1)


def test_paired_separate():
	"""test separate trimming of paired-end reads"""
	run('-a TTAGACATAT', 'paired-separate.1.fastq', 'paired.1.fastq')
	run('-a CAGTGGAGTA', 'paired-separate.2.fastq', 'paired.2.fastq')


def test_paired_end_legacy(cores):
	"""--paired-output, not using -A/-B/-G"""
	# the -m 14 filters out one read, which should then also be filtered out in the second output file
	# -q 10 should not change anything: qualities in file 1 are high enough,
	# qualities in file 2 should not be inspected.
	run_paired(
		'-a TTAGACATAT -m 14 -q 10',
		in1='paired.1.fastq', in2='paired.2.fastq',
		expected1='paired.m14.1.fastq', expected2='paired.m14.2.fastq',
		cores=cores
	)


def test_untrimmed_paired_output():
	with temporary_path("tmp-untrimmed.1.fastq") as untrimmed1:
		with temporary_path("tmp-untrimmed.2.fastq") as untrimmed2:
			run_paired(
				['-a', 'TTAGACATAT',
					'--untrimmed-output', untrimmed1,
					'--untrimmed-paired-output', untrimmed2],
				in1='paired.1.fastq', in2='paired.2.fastq',
				expected1='paired-trimmed.1.fastq', expected2='paired-trimmed.2.fastq',
				cores=1
			)
			assert_files_equal(cutpath('paired-untrimmed.1.fastq'), untrimmed1)
			assert_files_equal(cutpath('paired-untrimmed.2.fastq'), untrimmed2)


def test_explicit_format_with_paired():
	# Use --format=fastq with input files whose extension is .txt
	with temporary_path("paired.1.txt") as txt1:
		with temporary_path("paired.2.txt") as txt2:
			shutil.copyfile(datapath("paired.1.fastq"), txt1)
			shutil.copyfile(datapath("paired.2.fastq"), txt2)
			run_paired(
				'--format=fastq -a TTAGACATAT -m 14',
				in1=txt1, in2=txt2,
				expected1='paired.m14.1.fastq',
				expected2='paired.m14.2.fastq',
				cores=1
			)


def test_no_trimming_legacy():
	# make sure that this doesn't divide by zero
	main([
		'-a', 'XXXXX', '-o', '/dev/null', '-p', '/dev/null',
		datapath('paired.1.fastq'), datapath('paired.2.fastq')])


def test_no_trimming():
	# make sure that this doesn't divide by zero
	main([
		'-a', 'XXXXX', '-A', 'XXXXX', '-o', '/dev/null', '-p', '/dev/null',
		datapath('paired.1.fastq'), datapath('paired.2.fastq')])


@raises(SystemExit)
def test_missing_file():
	with redirect_stderr():
		main(['-a', 'XX', '--paired-output', 'out.fastq', datapath('paired.1.fastq')])


def test_first_too_short(cores):
	with pytest.raises(SystemExit):
		with temporary_path("truncated.1.fastq") as trunc1:
			# Create a truncated file in which the last read is missing
			with open(datapath('paired.1.fastq')) as f:
				lines = f.readlines()
				lines = lines[:-4]
			with open(trunc1, 'w') as f:
				f.writelines(lines)
			with redirect_stderr():
				main(
					'-a XX -o /dev/null --paired-output out.fastq'.split()
					+ ['--cores', str(cores)]
					+ [trunc1, datapath('paired.2.fastq')]
				)


def test_second_too_short(cores):
	with pytest.raises(SystemExit):
		with temporary_path("truncated.2.fastq") as trunc2:
			# Create a truncated file in which the last read is missing
			with open(datapath('paired.2.fastq')) as f:
				lines = f.readlines()
				lines = lines[:-4]
			with open(trunc2, 'w') as f:
				f.writelines(lines)
			with redirect_stderr():
				main('-a XX -o /dev/null --paired-output out.fastq'.split()
					+ ['--cores', str(cores)]
					+ [datapath('paired.1.fastq'), trunc2])


def test_unmatched_read_names(cores):
	with pytest.raises(SystemExit):
		with temporary_path("swapped.1.fastq") as swapped:
			# Create a file in which reads 2 and 1 are swapped
			with open(datapath('paired.1.fastq')) as f:
				lines = f.readlines()
				lines = lines[0:4] + lines[8:12] + lines[4:8] + lines[12:]
			with open(swapped, 'w') as f:
				f.writelines(lines)
			with redirect_stderr():
				main('-a XX -o out1.fastq --paired-output out2.fastq'.split()
					+ ['--cores', str(cores)]
					+ [swapped, datapath('paired.2.fastq')])


def test_p_without_o(cores):
	"""Option -p given but -o missing"""
	with pytest.raises(SystemExit):
		main('-a XX -p /dev/null'.split()
			+ ['--cores', str(cores)]
			+ [datapath('paired.1.fastq'), datapath('paired.2.fastq')])


def test_paired_but_only_one_input_file(cores):
	"""Option -p given but only one input file"""
	with pytest.raises(SystemExit):
		main('-a XX -o /dev/null -p /dev/null'.split()
			+ ['--cores', str(cores)]
			+ [datapath('paired.1.fastq')])


def test_legacy_minlength(cores):
	"""Ensure -m is not applied to second read in a pair in legacy mode"""
	run_paired(
		'-a XXX -m 27',
		in1='paired.1.fastq', in2='paired.2.fastq',
		expected1='paired-m27.1.fastq', expected2='paired-m27.2.fastq',
		cores=cores
	)


def test_paired_end(cores):
	"""single-pass paired-end with -m"""
	run_paired(
		'-a TTAGACATAT -A CAGTGGAGTA -m 14',
		in1='paired.1.fastq', in2='paired.2.fastq',
		expected1='paired.1.fastq', expected2='paired.2.fastq',
		cores=cores
	)


def test_paired_anchored_back_no_indels():
	run_paired(
		'-a BACKADAPTER$ -A BACKADAPTER$ -N --no-indels',
		in1='anchored-back.fasta', in2='anchored-back.fasta',
		expected1='anchored-back.fasta', expected2="anchored-back.fasta",
		cores=1
	)


def test_paired_end_qualtrim(cores):
	"""single-pass paired-end with -q and -m"""
	run_paired(
		'-q 20 -a TTAGACATAT -A CAGTGGAGTA -m 14 -M 90',
		in1='paired.1.fastq', in2='paired.2.fastq',
		expected1='pairedq.1.fastq', expected2='pairedq.2.fastq',
		cores=cores
	)


def test_paired_end_qualtrim_swapped(cores):
	"""single-pass paired-end with -q and -m, but files swapped"""
	run_paired(
		'-q 20 -a CAGTGGAGTA -A TTAGACATAT -m 14',
		in1='paired.2.fastq', in2='paired.1.fastq',
		expected1='pairedq.2.fastq', expected2='pairedq.1.fastq',
		cores=cores
	)


def test_paired_end_cut(cores):
	run_paired(
		'-u 3 -u -1 -U 4 -U -2',
		in1='paired.1.fastq', in2='paired.2.fastq',
		expected1='pairedu.1.fastq', expected2='pairedu.2.fastq',
		cores=cores
	)


def test_paired_end_upper_a_only(cores):
	run_paired(
		'-A CAGTGGAGTA',
		in1='paired.1.fastq', in2='paired.2.fastq',
		expected1='paired-onlyA.1.fastq', expected2='paired-onlyA.2.fastq',
		cores=cores
	)


def test_discard_untrimmed(cores):
	# issue #146
	# the first adapter is a sequence cut out from the first read
	run_paired(
		'-a CTCCAGCTTAGACATATC -A XXXXXXXX --discard-untrimmed',
		in1='paired.1.fastq', in2='paired.2.fastq',
		expected1='empty.fastq', expected2='empty.fastq',
		cores=cores
	)


def test_discard_trimmed(cores):
	run_paired(
		'-A C -O 1 --discard-trimmed',  # applies everywhere
		in1='paired.1.fastq', in2='paired.2.fastq',
		expected1='empty.fastq', expected2='empty.fastq',
		cores=cores
	)


def test_interleaved_in_and_out(cores):
	"""Single-pass interleaved paired-end with -q and -m"""
	run_interleaved(
		'-q 20 -a TTAGACATAT -A CAGTGGAGTA -m 14 -M 90',
		inpath1='interleaved.fastq', expected1='interleaved.fastq',
		cores=cores
	)


def test_interleaved_in(cores):
	"""Interleaved input, two files output"""
	run_interleaved(
		'-q 20 -a TTAGACATAT -A CAGTGGAGTA -m 14 -M 90',
		inpath1='interleaved.fastq',
		expected1='pairedq.1.fastq', expected2='pairedq.2.fastq',
		cores=cores
	)


def test_interleaved_out(cores):
	"""Two files input, interleaved output"""
	run_interleaved(
		'-q 20 -a TTAGACATAT -A CAGTGGAGTA -m 14 -M 90',
		inpath1='paired.1.fastq', inpath2='paired.2.fastq',
		expected1='interleaved.fastq',
		cores=cores
	)


@raises(SystemExit)
def test_interleaved_neither_nor():
	"""Option --interleaved used, but pairs of files given for input and output"""
	with temporary_path("temp-paired.1.fastq") as p1:
		with temporary_path("temp-paired.2.fastq") as p2:
			params = '-a XX --interleaved'.split()
			with redirect_stderr():
				params += ['-o', p1, '-p1', p2, 'paired.1.fastq', 'paired.2.fastq']
				main(params)


def test_pair_filter(cores):
	run_paired(
		'--pair-filter=both -a TTAGACATAT -A GGAGTA -m 14',
		in1='paired.1.fastq', in2='paired.2.fastq',
		expected1='paired-filterboth.1.fastq', expected2='paired-filterboth.2.fastq',
		cores=cores
	)


def test_too_short_paired_output():
	with temporary_path("temp-too-short.1.fastq") as p1:
		with temporary_path("temp-too-short.2.fastq") as p2:
			run_paired(
				'-a TTAGACATAT -A CAGTGGAGTA -m 14 --too-short-output '
				'{0} --too-short-paired-output {1}'.format(p1, p2),
				in1='paired.1.fastq', in2='paired.2.fastq',
				expected1='paired.1.fastq', expected2='paired.2.fastq',
				cores=1
			)
			assert_files_equal(cutpath('paired-too-short.1.fastq'), p1)
			assert_files_equal(cutpath('paired-too-short.2.fastq'), p2)


def test_too_long_output():
	with temporary_path('temp-too-long.1.fastq') as p1:
		with temporary_path('temp-too-long.2.fastq') as p2:
			run_paired(
				'-a TTAGACATAT -A CAGTGGAGTA -M 14 --too-long-output '
				'{0} --too-long-paired-output {1}'.format(p1, p2),
				in1='paired.1.fastq', in2='paired.2.fastq',
				expected1='paired-too-short.1.fastq', expected2='paired-too-short.2.fastq',
				cores=1
			)
			assert_files_equal(cutpath('paired.1.fastq'), p1)
			assert_files_equal(cutpath('paired.2.fastq'), p2)


@raises(SystemExit)
def test_too_short_output_paired_option_missing():
	with temporary_path('temp-too-short.1.fastq') as p1:
		run_paired(
			'-a TTAGACATAT -A CAGTGGAGTA -m 14 --too-short-output '
			'{0}'.format(p1),
			in1='paired.1.fastq', in2='paired.2.fastq',
			expected1='paired.1.fastq', expected2='paired.2.fastq',
			cores=1
		)


def test_nextseq_paired(cores):
	run_paired('--nextseq-trim 22', in1='nextseq.fastq', in2='nextseq.fastq',
		expected1='nextseq.fastq', expected2='nextseq.fastq',
		cores=cores)


def test_paired_demultiplex():
	tempdir = tempfile.mkdtemp(prefix='cutadapt-tests.')
	multiout1 = os.path.join(tempdir, 'demultiplexed.{name}.1.fastq')
	multiout2 = os.path.join(tempdir, 'demultiplexed.{name}.2.fastq')
	params = [
		'-a', 'first=AACATTAGACA', '-a', 'second=CATTAGACATATCGG',
		'-A', 'ignored=CAGTGGAGTA', '-A', 'alsoignored=AATAACAGTGGAGTA',
		'-o', multiout1, '-p', multiout2,
		datapath('paired.1.fastq'), datapath('paired.2.fastq')]
	assert main(params) is None
	assert_files_equal(cutpath('demultiplexed.first.1.fastq'), multiout1.format(name='first'))
	assert_files_equal(cutpath('demultiplexed.second.1.fastq'), multiout1.format(name='second'))
	assert_files_equal(cutpath('demultiplexed.unknown.1.fastq'), multiout1.format(name='unknown'))
	assert_files_equal(cutpath('demultiplexed.first.2.fastq'), multiout2.format(name='first'))
	assert_files_equal(cutpath('demultiplexed.second.2.fastq'), multiout2.format(name='second'))
	assert_files_equal(cutpath('demultiplexed.unknown.2.fastq'), multiout2.format(name='unknown'))
	shutil.rmtree(tempdir)
