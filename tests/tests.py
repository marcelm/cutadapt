# TODO
# test --untrimmed-output
# test with the --output option
# test reading from standard input
from __future__ import print_function, division
import os
from cutadapt.scripts import cutadapt


def dpath(path):
	"""
	get path to a data file (relative to the directory this test lives in)
	"""
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

	assert cutadapt.main(params) == 0
	# TODO redirect standard output
	diff(dpath(os.path.join('cut', expected)), dpath('tmp.fastaq'))
	os.remove(dpath('tmp.fastaq'))
	# TODO diff log files
	#echo "Running $CA $1 data/$3 ${second}"
	#if ! $CA $1 "data/$3" -o tmp.fastaq ${second} > tmp.log; then
		#cat tmp.log
		#exit 1
	#fi
	#sed -i '/Total time/d;/Time per read/d;/cutadapt version/d;/^Command line /d' tmp.log
	#diff -u cut/$2 tmp.fastaq
	#diff -u tmp.log log/$2.log


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
	diff(datapath('rest.txt'), dpath('rest.tmp'))
	os.remove(dpath('rest.tmp'))


def test_restfront():
	run(['-g', 'ADAPTER', '-r', dpath('rest.tmp')], "restfront.fa", "rest.fa")
	diff(datapath('restfront.txt'), dpath('rest.tmp'))
	os.remove(dpath('rest.tmp'))


def test_discard():
	'''--discard'''
	run("-b TTAGACATATCTCCGTCG --discard", "discard.fastq", "small.fastq")


def test_discard_untrimmed():
	'''--discard-untrimmed'''
	run('-b CAAGAT --discard-untrimmed', 'discard-untrimmed.fastq', 'small.fastq')


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


def test_too_short_no_primer():
	'''--too-short-output and --trim-primer'''
	run("-c -m 5 -a 330201030313112312 --trim-primer --too-short-output tooshort.tmp.fa", "minlen.noprimer.fa", "minlen.fa")
	diff(datapath('tooshort.noprimer.fa'), "tooshort.tmp.fa")
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

def test_qualtrim_csfastaqual():
	'''-q with csfasta/qual files'''
	run("-c -q 10", "solidqual.fastq", "solid.csfasta", 'solid.qual')

def test_qualbase():
	'''-q with low qualities, using ascii(quality+64) encoding'''
	run("-q 10 --quality-base 64 -a XXXXXX", "illumina64.fastq", "illumina64.fastq")

def test_quality_trim_only():
	'''only trim qualities, do not remove adapters'''
	run("-q 10 --quality-base 64", "illumina64.fastq", "illumina64.fastq")

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
	'''wildcards in adapter'''
	wildcardtmp = dpath("wildcardtmp.txt")
	for adapter_type, expected in (("-a", "wildcard_adapter.fa"),
		("-b", "wildcard_adapter_anywhere.fa")):
		run("--wildcard-file {0} {1} ACGTNNNACGT".format(wildcardtmp, adapter_type),
			expected, "wildcard_adapter.fa")
		lines = open(wildcardtmp).readlines()
		lines = [ line.strip() for line in lines ]
		assert lines == ['AAA 1', 'GGG 2', 'CCC 3b', 'TTT 4b']
		os.remove(wildcardtmp)

def test_wildcard_N():
	'''test 'N' wildcard matching with no allowed errors'''
	run("-e 0 -a GGGGGGG --match-read-wildcards", "wildcardN.fa", "wildcardN.fa")

def test_adapter_front():
	'''test adapter in front'''
	run("--front ADAPTER", "examplefront.fa", "example.fa")

def test_literal_N():
	'''test matching literal 'N's'''
	run("-N -e 0.2 -a NNNNNNNNNNNNNN", "trimN3.fasta", "trimN3.fasta")

def test_literal_N2():
	run("-N -O 1 -g NNNNNNNNNNNNNN", "trimN5.fasta", "trimN5.fasta")

def test_anchored_front():
	run("-g ^FRONTADAPT", "anchored.fasta", "anchored.fasta")

def test_solid():
	run("-c -e 0.122 -a 330201030313112312", "solid.fastq", "solid.fastq")

def test_solid_basespace_adapter():
	'''colorspace adapter given in basespace'''
	run("-c -e 0.122 -a CGCCTTGGCCGTACAGCAG", "solid.fastq", "solid.fastq")

def test_solid5p():
	'''test 5' colorspace adapter'''
	# this is not a real adapter, just a random string
	# in colorspace: C0302201212322332333
	run("-c -e 0.1 --trim-primer -g CCGGAGGTCAGCTCGCTATA", "solid5p.fasta", "solid5p.fasta")

def test_solid5p_prefix_notrim():
	'''test anchored 5' colorspace adapter, no primer trimming'''
	run("-c -e 0.1 -g ^CCGGAGGTCAGCTCGCTATA", "solid5p-anchored.notrim.fasta", "solid5p.fasta")

def test_solid5p_prefix():
	'''test anchored 5' colorspace adapter'''
	run("-c -e 0.1 --trim-primer -g ^CCGGAGGTCAGCTCGCTATA", "solid5p-anchored.fasta", "solid5p.fasta")

def test_solid5p_fastq():
	'''test 5' colorspace adapter'''
	# this is not a real adapter, just a random string
	# in colorspace: C0302201212322332333
	run("-c -e 0.1 --trim-primer -g CCGGAGGTCAGCTCGCTATA", "solid5p.fastq", "solid5p.fastq")

def test_solid5p_prefix_notrim_fastq():
	'''test anchored 5' colorspace adapter, no primer trimming'''
	run("-c -e 0.1 -g ^CCGGAGGTCAGCTCGCTATA", "solid5p-anchored.notrim.fastq", "solid5p.fastq")

def test_solid5p_prefix_fastq():
	'''test anchored 5' colorspace adapter'''
	run("-c -e 0.1 --trim-primer -g ^CCGGAGGTCAGCTCGCTATA", "solid5p-anchored.fastq", "solid5p.fastq")

def test_sra_fastq():
	'''test SRA-formatted colorspace FASTQ'''
	run("-c -e 0.1 --format sra-fastq -a CGCCTTGGCCGTACAGCAG", "sra.fastq", "sra.fastq")

def test_issue_46():
	'''issue 46 - IndexError with --wildcard-file'''
	wildcardtmp = dpath("wildcardtmp.txt")
	run("--anywhere=AACGTN --wildcard-file={0}".format(wildcardtmp), "issue46.fasta", "issue46.fasta")
	os.remove(wildcardtmp)

def test_strip_suffix():
	run("--strip-suffix _sequence -a XXXXXXX", "stripped.fasta", "simple.fasta")

# note: the actual adapter sequence in the illumina.fastq.gz data set is
# GCCTAACTTCTTAGACTGCCTTAAGGACGT (fourth base is different)
def test_info_file():
	infotmp = dpath("infotmp.txt")
	run("--info-file {0} -a GCCGAACTTCTTAGACTGCCTTAAGGACGT".format(infotmp), "illumina.fastq", "illumina.fastq.gz")
	os.remove(infotmp)

def test_named_adapter():
	run("-a MY_ADAPTER=GCCGAACTTCTTAGACTGCCTTAAGGACGT", "illumina.fastq", "illumina.fastq.gz")
