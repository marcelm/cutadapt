# coding: utf-8
from __future__ import print_function, division, absolute_import

from cutadapt.colorspace import encode, decode
from cutadapt.scripts.cutadapt import main
from .utils import run, datapath

# If there are any unknown characters in the test sequence,
# round tripping will only work if all characters after the
# first unknown character are also unknown:
# encode("TNGN") == "T444", but
# decode("T444") == "TNNN".

sequences = [
	"",
	"C",
	"ACGGTC",
	"TN",
	"TN.",
	"TNN.N",
	"CCGGCAGCATTCATTACGACAACGTGGCACCGTGTTTTCTCGGTGGTA",
	"TGCAGTTGATGATCGAAGAAAACGACATCATCAGCCAGCAAGTGC",
	"CAGGGTTTGATGAGTGGCTGTGGGTGCTGGCGTATCCGGG"
	]


def test_encode():
	assert encode("AA") == "A0"
	assert encode("AC") == "A1"
	assert encode("AG") == "A2"
	assert encode("AT") == "A3"
	assert encode("CA") == "C1"
	assert encode("CC") == "C0"
	assert encode("CG") == "C3"
	assert encode("CT") == "C2"
	assert encode("GA") == "G2"
	assert encode("GC") == "G3"
	assert encode("GG") == "G0"
	assert encode("GT") == "G1"
	assert encode("TA") == "T3"
	assert encode("TC") == "T2"
	assert encode("TG") == "T1"
	assert encode("TT") == "T0"

	assert encode("TN") == "T4"
	assert encode("NT") == "N4"
	assert encode("NN") == "N4"

	assert encode("ACGGTC") == "A13012"
	assert encode("TTT.N") == "T0044"
	assert encode("TTNT.N") == "T04444"


def test_decode():
	for s in sequences:
		expected = s.replace('.', 'N')
		encoded = encode(s)
		assert decode(encoded) == expected
	assert decode('A.') == 'AN'
	assert decode('C.') == 'CN'
	assert decode('G.') == 'GN'
	assert decode('T.') == 'TN'


def test_qualtrim_csfastaqual():
	"""-q with csfasta/qual files"""
	run("-c -q 10", "solidqual.fastq", "solid.csfasta", 'solid.qual')


def test_E3M():
	"""Read the E3M dataset"""
	# not really colorspace, but a fasta/qual file pair
	main(['-o', '/dev/null', datapath("E3M.fasta"), datapath("E3M.qual")])


def test_bwa():
	"""MAQ-/BWA-compatible output"""
	run("-c -e 0.12 -a 330201030313112312 -x 552: --maq", "solidmaq.fastq", "solid.csfasta", 'solid.qual')


def test_bfast():
	"""BFAST-compatible output"""
	run("-c -e 0.12 -a 330201030313112312 -x abc: --strip-f3", "solidbfast.fastq", "solid.csfasta", 'solid.qual')


def test_trim_095():
	"""some reads properly trimmed since cutadapt 0.9.5"""
	run("-c -e 0.122 -a 330201030313112312", "solid.fasta", "solid.fasta")


def test_solid():
	run("-c -e 0.122 -a 330201030313112312", "solid.fastq", "solid.fastq")


def test_solid_basespace_adapter():
	"""colorspace adapter given in basespace"""
	run("-c -e 0.122 -a CGCCTTGGCCGTACAGCAG", "solid.fastq", "solid.fastq")


def test_solid5p():
	"""test 5' colorspace adapter"""
	# this is not a real adapter, just a random string
	# in colorspace: C0302201212322332333
	run("-c -e 0.1 --trim-primer -g CCGGAGGTCAGCTCGCTATA", "solid5p.fasta", "solid5p.fasta")


def test_solid5p_prefix_notrim():
	"""test anchored 5' colorspace adapter, no primer trimming"""
	run("-c -e 0.1 -g ^CCGGAGGTCAGCTCGCTATA", "solid5p-anchored.notrim.fasta", "solid5p.fasta")


def test_solid5p_prefix():
	"""test anchored 5' colorspace adapter"""
	run("-c -e 0.1 --trim-primer -g ^CCGGAGGTCAGCTCGCTATA", "solid5p-anchored.fasta", "solid5p.fasta")


def test_solid5p_fastq():
	"""test 5' colorspace adapter"""
	# this is not a real adapter, just a random string
	# in colorspace: C0302201212322332333
	run("-c -e 0.1 --trim-primer -g CCGGAGGTCAGCTCGCTATA", "solid5p.fastq", "solid5p.fastq")


def test_solid5p_prefix_notrim_fastq():
	"""test anchored 5' colorspace adapter, no primer trimming"""
	run("-c -e 0.1 -g ^CCGGAGGTCAGCTCGCTATA", "solid5p-anchored.notrim.fastq", "solid5p.fastq")


def test_solid5p_prefix_fastq():
	"""test anchored 5' colorspace adapter"""
	run("-c -e 0.1 --trim-primer -g ^CCGGAGGTCAGCTCGCTATA", "solid5p-anchored.fastq", "solid5p.fastq")


def test_sra_fastq():
	"""test SRA-formatted colorspace FASTQ"""
	run("-c -e 0.1 --format sra-fastq -a CGCCTTGGCCGTACAGCAG", "sra.fastq", "sra.fastq")


def test_no_zero_cap():
	run("--no-zero-cap -c -e 0.122 -a CGCCTTGGCCGTACAGCAG", "solid-no-zerocap.fastq", "solid.fastq")
