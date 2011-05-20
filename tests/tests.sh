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
	if ! $CA $1 "data/$3" ${second} > tmp.fastaq 2> tmp.log; then
		cat tmp.log
		exit 1
	fi
	sed -i '/Total time/d;/Time per read/d;/cutadapt version/d;/^Command line /d' tmp.log
	diff -u cut/$2 tmp.fastaq
	diff -u tmp.log log/$2.log
	rm tmp.fastaq tmp.log
}

test_cutadapt "-b TTAGACATATCTCCGTCG" small.fastq small.fastq

# empty input
test_cutadapt "-a TTAGACATATCTCCGTCG" empty.fastq empty.fastq

# lower case adapter
test_cutadapt "-b ttagacatatctccgtcg" lowercase.fastq small.fastq

# test -r/--rest-file
test_cutadapt "-b ADAPTER -r rest.tmp" rest.fa rest.fa
diff -u rest.tmp data/rest.txt

# test --discard
test_cutadapt "-b TTAGACATATCTCCGTCG --discard" discard.fastq small.fastq

# test if sequence name after the "+" is retained
test_cutadapt "-e 0.12 -b TTAGACATATCTCCGTCG" plus.fastq plus.fastq

# test the -f/--format parameter
test_cutadapt "-f fastq -b TTAGACATATCTCCGTCG" small.fastq small.myownextension

# test -m/--minimum-length
test_cutadapt "-c -m 5 -a 330201030313112312" minlen.fa minlen.fa

# test --too-short-output
test_cutadapt "-c -m 5 -a 330201030313112312 --too-short-output tooshort.tmp.fa" minlen.fa minlen.fa
diff -u data/tooshort.fa tooshort.tmp.fa
rm tooshort.tmp.fa

# test -M/--maximum-length
test_cutadapt "-c -M 5 -a 330201030313112312" maxlen.fa maxlen.fa

# test 454 data; test -n and --length-tag
test_cutadapt "-n 3 -e 0.1 --length-tag length=
	-b TGAGACACGCAACAGGGGAAAGGCAAGGCACACAGGGGATAGG
	-b TCCATCTCATCCCTGCGTGTCCCATCTGTTCCCTCCCTGTCTCA" 454.fa 454.fa

echo "Tests passed"
