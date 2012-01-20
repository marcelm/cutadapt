from __future__ import print_function, division
import os
import imp

# import cutadapt although it does not have a .py extension.
# better solution:
# - symlinks
# - move code into modules, make cutadapt a simple wrapper
ca = imp.load_source('ca', 'cutadapt')


def dpath(path):
	"""get path to a data file (relative to the directory this
	test lives in)"""
	return os.path.join(os.path.dirname(__file__), path)

def datapath(path):
	return dpath(os.path.join('data', path))


def diff(path1, path2):
	assert os.system("diff -u {0} {1}".format(path1, path2)) == 0


def run(params, expected, inpath, inpath2=None):
	if type(params) is str:
		params = params.split()
	params += ['-o', dpath('tmp.fastaq') ] # TODO not parallelizable
	params += [ datapath(inpath) ]
	if inpath2:
		params += [ datapath(inpath2) ]

	ca.main(params)
	# TODO redirect standard output
	ret = os.system('diff -u ' + dpath(os.path.join('cut', expected)) + ' ' + dpath('tmp.fastaq'))
	assert ret == 0
	# TODO diff log files
	#echo "Running $CA $1 data/$3 ${second}"
	#if ! $CA $1 "data/$3" -o tmp.fastaq ${second} > tmp.log; then
		#cat tmp.log
		#exit 1
	#fi
	#sed -i '/Total time/d;/Time per read/d;/cutadapt version/d;/^Command line /d' tmp.log
	#diff -u cut/$2 tmp.fastaq
	#diff -u tmp.log log/$2.log
	#rm tmp.fastaq tmp.log


def test_example():
	run(["-b", "ADAPTER"], 'example.fa', 'example.fa')

def test_small():
	run(["-b", "TTAGACATATCTCCGTCG"], 'small.fastq', 'small.fastq')

def test_empty():
	'''empty input'''
	run(["-a", "TTAGACATATCTCCGTCG"], 'empty.fastq', 'empty.fastq')

def test_newlines():
	'''DOS/Windows newlines'''
	run("-e 0.12 -b TTAGACATATCTCCGTCG", "dos.fastq", "dos.fastq")

def test_lowercase():
	'''lower case adapter'''
	run("-b ttagacatatctccgtcg", "lowercase.fastq", "small.fastq")


def test_rest():
	'''-r/--rest-file'''
	run(['-b', 'ADAPTER', '-r', dpath('rest.tmp')], "rest.fa", "rest.fa")
	diff(dpath('rest.tmp'), datapath('rest.txt'))
	os.remove(dpath('rest.tmp'))


def test_restfront():
	run(['-g', 'ADAPTER', '-r', dpath('rest.tmp')], "restfront.fa", "rest.fa")
	diff(dpath('rest.tmp'), datapath('restfront.txt'))
	os.remove(dpath('rest.tmp'))


def test_discard():
	'''--discard'''
	run("-b TTAGACATATCTCCGTCG --discard", "discard.fastq", "small.fastq")


def test_plus():
	'''test if sequence name after the "+" is retained'''
	run("-e 0.12 -b TTAGACATATCTCCGTCG", "plus.fastq", "plus.fastq")


def test_extensiontxtgz():
	'''automatic recognition of "_sequence.txt.gz" extension'''
	run("-b TTAGACATATCTCCGTCG", "s_1_sequence.txt", "s_1_sequence.txt.gz")


def test_format():
	'''the -f/--format parameter'''
	run("-f fastq -b TTAGACATATCTCCGTCG", "small.fastq", "small.myownextension")


def test_minimum_length():
	'''-m/--minimum-length'''
	run("-c -m 5 -a 330201030313112312", "minlen.fa", "minlen.fa")


def test_too_short():
	'''--too-short-output'''
	run("-c -m 5 -a 330201030313112312 --too-short-output tooshort.tmp.fa", "minlen.fa", "minlen.fa")
	diff(datapath('tooshort.fa'), "tooshort.tmp.fa")
	os.remove('tooshort.tmp.fa')


def test_maximum_length():
	'''-M/--maximum-length'''
	run("-c -M 5 -a 330201030313112312", "maxlen.fa", "maxlen.fa")


def test_length_tag():
	'''454 data; -n and --length-tag'''
	run("-n 3 -e 0.1 --length-tag length= " \
		"-b TGAGACACGCAACAGGGGAAAGGCAAGGCACACAGGGGATAGG "\
		"-b TCCATCTCATCCCTGCGTGTCCCATCTGTTCCCTCCCTGTCTCA", '454.fa', '454.fa')


def test_overlap_a():
	'''-O/--overlap with -a (-c omitted on purpose)'''
	run("-O 10 -a 330201030313112312", "overlapa.fa", "overlapa.fa")

def test_overlap_b():
	'''-O/--overlap with -b'''
	run("-O 10 -b TTAGACATATCTCCGTCG", "overlapb.fa", "overlapb.fa")

def test_qualtrim():
	'''-q with low qualities'''
	run("-q 10 -a XXXXXX", "lowqual.fastq", "lowqual.fastq")

def test_qualbase():
	'''-q with low qualities, using ascii(quality+64) encoding'''
	run("-q 10 --quality-base 64 -a XXXXXX", "illumina64.fastq", "illumina64.fastq")

def test_twoadapters():
	'''two adapters'''
	run("-b CTCGAGAATTCTGGATCCTC -b GAGGATCCAGAATTCTCGAGTT", "twoadapters.fasta", "twoadapters.fasta")

def test_bwa():
	'''MAQ-/BWA-compatible output'''
	run("-c -e 0.12 -a 330201030313112312 -x 552: --maq", "solidmaq.fastq", "solid.csfasta", 'solid.qual')

def test_bfast():
	'''BFAST-compatible output'''
	run("-c -e 0.12 -a 330201030313112312 -x abc: --strip-f3", "solidbfast.fastq", "solid.csfasta", 'solid.qual')

def test_polya():
	'''poly-A tails'''
	run("-m 24 -O 10 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "polya.fasta", "polya.fasta")

def test_trim_095():
	'''some reads properly trimmed since cutadapt 0.9.5'''
	run("-c -e 0.122 -a 330201030313112312", "solid.fasta", "solid.fasta")

def test_gz_multiblock():
	'''compressed gz file with multiple blocks (created by concatenating two .gz files)'''
	run("-b TTAGACATATCTCCGTCG", "small.fastq", "multiblock.fastq.gz")

def test_suffix():
	'''-y/--suffix parameter, combined with _F3'''
	run("-c -e 0.12 -a 330201030313112312 -y _my_suffix --strip-f3", "suffix.fastq", "solid.csfasta", 'solid.qual')

def test_read_wildcard():
	'''test wildcards in reads'''
	run("--match-read-wildcards -b ACGTACGT", "wildcard.fa", "wildcard.fa")

def test_adapter_wildcard():
	'''test wildcards in adapter'''
	run("--wildcard-file - -b ACGTNNNACGT", "wildcard_adapter.fa", "wildcard_adapter.fa")

def test_wildcard_N():
	'''test 'N' wildcard matching with no allowed errors='''
	run("-e 0 -a GGGGGGG --match-read-wildcards", "wildcardN.fa", "wildcardN.fa")

def test_adapter_front():
	'''test adapter in front'''
	run("--front ADAPTER", "examplefront.fa", "example.fa")

def test_literal_N():
	'''test matching literal 'N's'''
	run("-N -e 0.2 -a NNNNNNNNNNNNNN", "trimN3.fasta", "trimN3.fasta")

def test_literal_N2():
	run("-N -O 1 -g NNNNNNNNNNNNNN", "trimN5.fasta", "trimN5.fasta")
