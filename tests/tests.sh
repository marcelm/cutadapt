#!/bin/bash
# Test cases for cutadapt. Run this script within the tests/ directory.
set -e

# path to the script to test
CA="../cutadapt"

function test_cutadapt() {
	# parameters:
	# 1. command-line parameters to cutadapt
	# 2. name of expected output
	# 3. input file
	# 4. optional: second input file (.qual file)

	params="$1"
	if [ x$4 != x ]; then
		second="data/$4"
	else
		second=""
	fi
	echo "Running $CA $1 data/$3 ${second}"
	if ! $CA $1 "data/$3" -o tmp.fastaq ${second} > tmp.log; then
		cat tmp.log
		exit 1
	fi
	sed -i '/Total time/d;/Time per read/d;/cutadapt version/d;/^Command line /d' tmp.log
	diff -u cut/$2 tmp.fastaq
	diff -u tmp.log log/$2.log
	rm tmp.fastaq tmp.log
	#tmp.log0
}

test_cutadapt "-b ADAPTER" example.fa example.fa

test_cutadapt "-b TTAGACATATCTCCGTCG" small.fastq small.fastq

# empty input
test_cutadapt "-a TTAGACATATCTCCGTCG" empty.fastq empty.fastq

# DOS/Windows newlines
test_cutadapt "-e 0.12 -b TTAGACATATCTCCGTCG" dos.fastq dos.fastq

# lower case adapter
test_cutadapt "-b ttagacatatctccgtcg" lowercase.fastq small.fastq

# -r/--rest-file
test_cutadapt "-b ADAPTER -r rest.tmp" rest.fa rest.fa
diff -u rest.tmp data/rest.txt
rm rest.tmp

# --discard
test_cutadapt "-b TTAGACATATCTCCGTCG --discard" discard.fastq small.fastq

# test if sequence name after the "+" is retained
test_cutadapt "-e 0.12 -b TTAGACATATCTCCGTCG" plus.fastq plus.fastq

# automatic recognition of "_sequence.txt.gz" extension
test_cutadapt "-b TTAGACATATCTCCGTCG" s_1_sequence.txt s_1_sequence.txt.gz

# the -f/--format parameter
test_cutadapt "-f fastq -b TTAGACATATCTCCGTCG" small.fastq small.myownextension

# -m/--minimum-length
test_cutadapt "-c -m 5 -a 330201030313112312" minlen.fa minlen.fa

# --too-short-output
test_cutadapt "-c -m 5 -a 330201030313112312 --too-short-output tooshort.tmp.fa" minlen.fa minlen.fa
diff -u data/tooshort.fa tooshort.tmp.fa
rm tooshort.tmp.fa

# -M/--maximum-length
test_cutadapt "-c -M 5 -a 330201030313112312" maxlen.fa maxlen.fa

# 454 data; -n and --length-tag
test_cutadapt "-n 3 -e 0.1 --length-tag length=
	-b TGAGACACGCAACAGGGGAAAGGCAAGGCACACAGGGGATAGG
	-b TCCATCTCATCCCTGCGTGTCCCATCTGTTCCCTCCCTGTCTCA" 454.fa 454.fa

# -O/--overlap with -a (-c omitted on purpose)
test_cutadapt "-O 10 -a 330201030313112312" overlapa.fa overlapa.fa

# -O/--overlap with -b
test_cutadapt "-O 10 -b TTAGACATATCTCCGTCG" overlapb.fa overlapb.fa

# -q with low qualities
test_cutadapt "-q 10 -a XXXXXX" lowqual.fastq lowqual.fastq

# -q with low qualities, using ascii(quality+64) encoding
test_cutadapt "-q 10 --quality-base 64 -a XXXXXX" illumina64.fastq illumina64.fastq

# two adapters
test_cutadapt "-b CTCGAGAATTCTGGATCCTC -b GAGGATCCAGAATTCTCGAGTT" twoadapters.fasta twoadapters.fasta

# MAQ-compatible output
test_cutadapt "-c -e 0.12 -a 330201030313112312 -x 552: --maq" solidmaq.fastq solid.csfasta solid.qual

# BFAST-compatible output
test_cutadapt "-c -e 0.12 -a 330201030313112312 -x abc: --strip-f3" solidbfast.fastq solid.csfasta solid.qual

# poly-A tails
test_cutadapt "-m 24 -O 10 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" polya.fasta polya.fasta

# some reads properly trimmed since cutadapt 0.9.5
test_cutadapt "-c -e 0.122 -a 330201030313112312" solid.fasta solid.fasta

# compressed gz file with multiple blocks (created by concatenating two .gz files)
test_cutadapt "-b TTAGACATATCTCCGTCG" small.fastq multiblock.fastq.gz

# -y/--suffix parameter, combined with _F3
test_cutadapt "-c -e 0.12 -a 330201030313112312 -y _my_suffix --strip-f3" suffix.fastq solid.csfasta solid.qual

# compressed gz file with multiple blocks (created by concatenating two .gz files)
test_cutadapt "-b TTAGACATATCTCCGTCG" small.fastq multiblock.fastq.gz

# -y/--suffix parameter, combined with _F3
test_cutadapt "-c -e 0.12 -a 330201030313112312 -y _my_suffix --strip-f3" suffix.fastq solid.csfasta solid.qual

# test read wildcards
test_cutadapt "--match-read-wildcards -b ACGTACGT" wildcard.fa wildcard.fa

# test adapter wildcards
test_cutadapt "--wildcard-file - -b ACGTNNNACGT" wildcard_adapter.fa wildcard_adapter.fa

# test 'N' wildcard matching with no allowed errors=
test_cutadapt "-e 0 -a GGGGGGG --match-read-wildcards" wildcardN.fa wildcardN.fa

# test adapter in front
test_cutadapt "--front ADAPTER" examplefront.fa example.fa

echo "Tests passed"
