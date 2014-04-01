cutadapt
========

cutadapt removes adapter sequences from DNA high-throughput sequencing data.
This is necessary when the reads are longer than the molecule that is sequenced,
such as in microRNA data.

cutadapt is implemented in Python. It comes with an extension module
written in Cython that implements the alignment algorithm.

License
-------

cutadapt is licensed under the MIT license (see the file `LICENSE`).


Project homepage
----------------

See http://code.google.com/p/cutadapt/ . Please use the Google code issue
tracker for bug reports and feature requests.


Dependencies
------------

cutadapt needs Python 2.6 or later (this includes Python 3). Python 2.6 supported
is not tested thoroughly and is also slower than later versions. For installation
from sources, a C compiler needs to be installed. The program has been developed
and tested on Ubuntu and OpenSuSE.


Installation
------------

Replace "python" with "python3" in the following lines to install the Python 3
version.

    python setup.py build
    python setup.py install

If you get an error about a missing "Python.h" file, then make sure that
the python-dev package is installed (or python3-dev for Python 3).


Use without installation
------------------------

Build the C extension module (you can try to skip this step -- a compiled
version of the module for Linux x86 is already included):

    python setup.py build_ext -i

Then simply run the script from where it is, similar to this:

    bin/cutadapt --help

If you get any errors, first try to explicitly request a specific Python
version by running cutadapt like this:

    python2.7 bin/cutadapt --help


Galaxy
------

If you want to use cutadapt within the web-based Galaxy platform
(http://galaxy.psu.edu/), please see the README file in the galaxy/ subfolder.
Galaxy support was contributed by Lance Parsons.


How to use, examples
====================

Please also see the command-line help:

    cutadapt --help

The basic command-line for cutadapt is:

    cutadapt -a AACCGGTT input.fastq > output.fastq

The adapter sequence is given with the `-a` option. Replace
`AACCGGTT` with your actual adapter sequence.

input.fastq is a file with reads. The result will be written
to standard output. Use redirection with `>` (or the `-o` option) to
write the output to a file.

cutadapt also writes a report after it has finished processing the
reads. If you use `-o`, the report is sent to standard output and to
stderr otherwise.

By default, the output file contains all reads, even those
that did not contain an adapter. (See also the `--discard` option.)

The following examples refer to basespace reads. See the "Colorspace"
section on how to use cutadapt with SOLiD reads.


Illumina data
-------------

Assuming your sequencing data is available as a FASTQ file, use this
command line:

    cutadapt -a ADAPTER-SEQUENCE input.fastq > output.fastq

gz-compressed input is supported:

    cutadapt -a ADAPTER-SEQUENCE input.fastq.gz > output.fastq

gz-compressed output is also supported, but the -o parameter (output file) needs
to be used (gzip compression is auto-detected by looking at the file name):

    cutadapt -a ADAPTER-SEQUENCE -o output.fastq.gz input.fastq.gz

If your Python installation includes support for bzip2 compression, then
bzip2-compressed files are also supported and recognized by their extension
`.bz2`.


Named adapters
--------------

You can give names to the adapters. The names are shown in addition to
the sequences themselves in the statistics overview when the program
has finished trimming the reads. You can use it like this:

    cutadapt -a My_adapter=ACGTAA input.fastq > output.fastq

Here, the actual adapter sequence is `ACGTAA` and the name assigned
to it is `My_adapter`.


FASTA file
----------

Cut an adapter from reads given in a FASTA file. Try to remove an adapter three times
(this is usually not needed), use the default error rate of 10%, write result to `output.fa`:

    cutadapt -n 3 -a TGAGACACGCAACAGGGGAAAGGCAAGGCACACAGGGGATAGG input.fa > output.fa


Multiple adapters
-----------------

As many adapters as desired can be given to the program by using the `-a`, `-b` or `-g`
in any combination, for example, five `-a` adapters and two `-g` adapters. All
adapters will be searched for, but only the best matching one will be trimmed
from each read (but see the `--times` option).

    cutadapt -b TGAGACACGCA -g AGGCACACAGGG input.fastq > output.fastq


Quality trimming
----------------

The `-q` (or `--trim-qualities`) parameter can be used to trim low-quality ends
from reads before adapter removal. For this to work correctly, the quality
values must be encoded as ascii(phred quality + 33). If they are encoded as
ascii(phred quality + 64), you need to add `--quality-base=64` to the command line.

The trimming algorithm is the same as the one used by BWA. That is: Subtract
the given cutoff from all qualities; compute partial sums from all indices to
the end of the sequence; cut sequence at the index at which the sum is minimal.


Paired-end adapter trimming
---------------------------

Cutadapt supports paired-end trimming, but currently two passes over the data are
required.

Assume the input is in `reads.1.fastq` and `reads.2.fastq` and that
`ADAPTER_FWD` should be trimmed from the forward reads (first file) and
`ADAPTER_REV` from the reverse reads (second file).

If you do not use any of the filtering options that discard reads,
such as `--discard`, `--minimum-length` or `--maximum-length`, then
run cutadapt on each file separately:

    cutadapt -a ADAPTER_FWD -o trimmed.1.fastq reads1.fastq
    cutadapt -a ADAPTER_REV -o trimmed.2.fastq reads2.fastq

You can use the options that are listed under 'Additional modifications' in cutadapt's help
output without problems. For example, if you want to quality-trim the first read in each
pair with a threshold of 10, and the second read in each pair with a threshold of 15, then
the commands could be:

    cutadapt -q 10 -a ADAPTER_FWD -o trimmed.1.fastq reads1.fastq
    cutadapt -q 15 -a ADAPTER_REV -o trimmed.2.fastq reads2.fastq

However, if you use one of the filtering options that discard reads, then you need to give
both input read files to cutadapt and the `--paired-output` option is needed to keep the two
files synchronized. First trim the forward read, writing output to temporary files (we
also add some quality trimming):

    cutadapt -q 10 -a ADAPTER_FWD --minimum-length 20 --paired-output tmp.2.fastq -o tmp.1.fastq reads.1.fastq reads.2.fastq

Then trim the reverse read, using the temporary files as input:

    cutadapt -q 15 -a ADAPTER_REV --minimum-length 20 --paired-output trimmed.1.fastq -o trimmed.2.fastq tmp.2.fastq tmp.1.fastq

Finally, remove the temporary files:

    rm tmp.1.fastq tmp.2.fastq

In each call to cutadapt, the read-modifying options such as `-q` only apply to the first
file (first reads.1.fastq, then tmp.2.fastq in this example). Reads in the second file
are not affected by those options, but by the filtering options: If a read in
the first file is discarded, then the matching read in the second file is also filtered
and not written to the output given by `--paired-output` in order to keep both output
files synchronized.


Illumina TruSeq
---------------

If you have reads containing Illumina TruSeq adapters, follow these steps.

Trim read 1 with `A` + the “TruSeq Indexed Adapter”. Use only the prefix of the adapter
sequence that is common to all Indexed Adapter sequences:

    cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o trimmed.1.fastq.gz reads.1.fastq.gz

Trim read 2 with the reverse complement of the ”TruSeq Universal Adapter”:

    cutadapt -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o trimmed.2.fastq.gz reads.2.fastq.gz

See also the section about paired-end adapter trimming above.

If you want to simplify this a bit, you can also use `AGATCGGAAGAGC` as the adapter
sequence in both cases:

    cutadapt -a AGATCGGAAGAGC -o trimmed.1.fastq.gz reads.1.fastq.gz
    cutadapt -a AGATCGGAAGAGC -o trimmed.2.fastq.gz reads.2.fastq.gz

The adapter sequences can be found in the document [Illumina TruSeq Adapters De-Mystified][1].

[1]: http://tucf-genomics.tufts.edu/documents/protocols/TUCF_Understanding_Illumina_TruSeq_Adapters.pdf


Adapters
========

These are some 454 adapters:

    A1:   5'- TCCATCTCATCCCTGCGTGTCCCATCTGTTCCCTCCCTGTCTCA
    A2:   5'- TGAGACAGGGAGGGAACAGATGGGACACGCAGGGATGAGATGGA
    B1:   5'- CCTATCCCCTGTGTGCCTTGCCTATCCCCTGTTGCGTGTCTCA
    B2:   5'- TGAGACACGCAACAGGGGAAAGGCAAGGCACACAGGGGATAGG

This is an AB SOLiD adapter (in color space) used in the SREK protocol:

    330201030313112312


Algorithm
=========

cutadapt uses a modified semi-global alignment algorithm. For speed, the
algorithm is implemented as a Cython extension module in `_align.pyx`.

Cutadapt’s processing speed is currently not dominated by the alignment
algorithm, but by parsing the input and writing the output.


Partial adapter matches
-----------------------

Cutadapt correctly deals with partial adapter matches, and also with any trailing
sequences after the adapter. As an example, suppose your adapter sequence is
"ADAPTER" (specified via the `-a` or `--adapter` command-line parameter).
If you have these input sequences:

    MYSEQUENCEADAPTER
    MYSEQUENCEADAP
    MYSEQUENCEADAPTERSOMETHINGELSE

All of them will be trimmed to "MYSEQUENCE". If the sequence starts with an
adapter, like this:

    ADAPTERSOMETHING

It will be empty after trimming.

When the allowed error rate is sufficiently high (set with parameter `-e`), errors in
the adapter sequence are allowed. For example, `ADABTER` (1 mismatch), `ADAPTR` (1 deletion),
and `ADAPPTER` (1 insertion) will all be recognized if the error rate is set to 0.15.


Anchoring 5' adapters
---------------------

If you specify an adapter with the `-g` (`--front`) parameter, the adapter may
overlap the beginning of the read or occur anywhere within it. If it appears
within the read, the sequence that precedes it will also be trimmed in addition
to the adapter. For example, with `-g ADAPTER`, these sequences:

    HELLOADAPTERTHERE
    APTERTHERE

will both be trimmed to `THERE`. To avoid this, you can prefix the adapter with the
character `^`. This will restrict the search, forcing the adapter to be a prefix
of the read. With `-g ^ADAPTER`, only reads like this will be trimmed:

    ADAPTERHELLO


Allowing adapters anywhere
--------------------------

Cutadapt assumes that any adapter specified via the `-a` (or `--adapter`) parameter
was ligated to the 3' end of the sequence. This is the correct assumption for
at least the SOLiD and Illumina small RNA protocols and probably others.
The assumption is enforced by the alignment algorithm, which only finds the adapter
when its starting position is within the read. In other words, the 5' base of
the adapter must appear within the read. The adapter and all bases following
it are removed.

If, on the other hand, your adapter can also be ligated to the 5' end (on
purpose or by accident), you should tell cutadapt so by using the `-b` (or
`--anywhere`) parameter. It will then use a slightly different alignment algorithm
(so-called semiglobal alignment), which allows any type of overlap between the
adapter and the sequence. In particular, the adapter may appear only partially
in the beginning of the read, like this:

    PTERMYSEQUENCE

The decision which part of the read to remove is made as follows: If there is at
least one base before the found adapter, then the adapter is considered to be
a 3' adapter and the adapter itself and everything following it is removed.
Otherwise, the adapter is considered to be a 5' adapter and it is removed from
the read.

Here are some examples, which may make this clearer (left: read, right: trimmed
read):

    MYSEQUENCEADAPTER -> MYSEQUENCE (3' adapter)
    MADAPTER -> M (3' adapter)
    ADAPTERMYSEQUENCE -> MYSEQUENCE (5' adapter)
    PTERMYSEQUENCE -> MYSEQUENCE (5' adapter)

The regular algorithm (`-a`) would trim the first two examples in the same way,
but trim the third to an empty sequence and trim the fourth not at all.

The `-b` parameter currently does not work with color space data.


Interpreting the statistics output
==================================

After every run, cutadapt prints out per-adapter statistics. The output starts
with something like this:

	Adapter 'ACGTACGTACGTTAGCTAGC', length 20, was trimmed 2402 times.

The meaning of this should be obvious.

The next piece of information is this:

	No. of allowed errors:
	0-9 bp: 0; 10-19 bp: 1; 20 bp: 2

The adapter has, as was conveniently shown above, a length of 20 characters.
We are using the default error rate of 0.1. What this implies is shown above:
Matches up to a length of 9 bp are allowed to have no errors. Matches of
lengths 10-19 bp are allowd to have 1 error and matches of length 20 can have
2 errors.

Finally, a table is output that gives more detailed information about the
lengths of the removed sequences. The following is only an excerpt; some
rows are left out:

	Overview of removed sequences
	length  count   expect  max.err error counts
	3       140     156.2   0       140
	4       57      39.1    0       57
	5       50      9.8     0       50
	6       35      2.4     0       35
	...
	100     397     0.0     3       358 36 3

The first row tells us the following: Three bases were removed in 140 reads;
randomly, one would expect this to occur 156.2 times;
the maximum number of errors at that match length is 0 (this is actually redundant since
we know already that no errors are allowed at lengths 0-9bp).

The last column shows the number of reads that had 0, 1, 2 ... errors.
In the last row, for example, 358 reads matched the adapter with zero errors, 36 with
1 error, and 3 matched with 2 errors.

The "expect" column gives only a rough estimate of the number of sequences
that is expected to match randomly (it assumes a GC content of 50%,
for example), but it can help to estimate whether the matches that were found are true
adapter matches or if they are due to chance. At lengths 6, for example, only 2.4 reads
are expected, but 35 do match, which hints that most of these matches are due to actual
adapters.

Note that the "length" column refers to the length of the removed sequence. That is,
the actual length of the match in the above row at length 100 is 20 since that is
the adapter length. Assuming the read length is 100, the adapter was found in the
beginning of 397 reads and therefore those reads were trimmed to a length of zero.

The table may also be useful in case the given adapter sequence contains an
error. In that case, it may look like this:

	...
	length  count   expect  max.err error counts
	10      53      0.0     1       51 2
	11      45      0.0     1       42 3
	12      51      0.0     1       48 3
	13      39      0.0     1       0 39
	14      40      0.0     1       0 40
	15      36      0.0     1       0 36
	...

We can see that no matches longer than 12 have zero errors. In this case, it
indicates that the 13th base of the given adapter sequence is incorrect.


Format of the info file
=======================

When the `--info-file` command-line parameter is given, detailed information
about the found adapters is written to the given file. The output is a
tab-separated text file. Each line corresponds to one read of the input file.
The columns are:
1. Read name
2. Number of errors
3. 0-based start coordinate of the adapter match
4. 0-based end coordinate of the adapter match
5. Sequence of the read before the adapter match
6. Sequence of the read that was matched to the adapter
7. Sequence of the read after the adapter match
8. Name of the found adapter.

The concatenation of the fields 5-6 yields the full read sequence. The adapter
name that should be used in column 8 can be given by writing `-a name=sequence`
instead of just `-a sequence`. Adapters without a name are numbered starting from 1.

If no adapter was found, the format is as follows:
* Read name
* The value -1
* The read sequence

When parsing that file, be aware that additional columns may be added in the
future. Also, in the current version, when the `--times` option is set to a
value other than 1 (the default value), multiple lines are written to the info
file for each read.


Using a "configuration file"
============================

Cutadapt currently does not support using a configuration file in which, for
example, a list of adapters can be specified. If you have many adapters that
you want to seach for and want to avoid typing all of them on the command line,
then you can use so-called "command substitution" of your Unix shell. With
Bash, this works as follows.

First, create a configuration file cutadapt.conf that contains lines like
this:

    -a AACCGGTT
    -a GTAATAACCGGTT
    -e 0.05
The file may contain line breaks (they will be replaced by spaces).

Then run cutadapt like this:

    cutadapt $(<cutadapt.conf) input.fastq > output.fastq

The Bash shell will replace the `$(<...)` with the content of the given file.


Colorspace
==========

Cutadapt was designed to work with colorspace reads from the ABi SOLiD sequencer.
Colorspace trimming is activated by the `--colorspace` option (or use `-c` for short).
The input reads can be given either:
* in a FASTA file
* in a FASTQ file
* in a `.csfasta` and a `.qual` file (this is the native SOLiD format).

In all cases, the colors must be represented by the characters 0, 1, 2, 3.
Example input files are in the cutadapt distribution at `tests/data/solid.*`.
The `.csfasta`/`.qual` file format is automatically assumed if two input
files are given to cutadapt.

In colorspace mode, the adapter sequences given to the `-a`, `-b` and `-g` options
can be given both as colors or as nucleotides. If given as nucleotides, they will
automatically be converted to colorspace. For example, to trim an adapter from
`solid.csfasta` and `solid.qual`, use this command-line:

	cutadapt -c -a CGCCTTGGCCGTACAGCAG solid.csfasta solid.qual > output.fastq

In case you know the colorspace adapter sequence, you can also write `330201030313112312`
instead of `CGCCTTGGCCGTACAGCAG` and the result is the same.


Ambiguity in colorspace
-----------------------

The ambiguity of colorspace encoding leads to some effects to be aware of when trimming
3' adapters from colorspace reads. For example, when trimming the adapter `AACTC`, cutadapt
searches for its colorspace-encoded version `0122`. But also `TTGAG`, `CCAGA` and `GGTCT`
have an encoding of `0122`. This means that effectively four different adapter sequences
are searched and trimmed at the same time. There is no way around this, unless
the decoded sequence were available, but that is usually only the case after read mapping.

The effect should usually be quite small. The number of false positives is multiplied by
four, but with a sufficiently large overlap (3 or 4 is already enough), this is still
only around 0.2 bases lost per read on average. If inspecting k-mer frequencies or using
small overlaps, you need to be aware of the effect, however.


Double-encoding, BWA and MAQ
----------------------------

The read mappers MAQ and BWA (and possibly others) need their colorspace
input reads to be in a so-called "double encoding". This simply means that
they cannot deal with the characters 0, 1, 2, 3 in the reads, but require
that the letters A, C, G, T be used for colors. For example, the colorspace
sequence `0011321` would be `AACCTGC` in double-encoded form. This is not
the same as conversion to basespace! The read is still in colorspace, only
letters are used instead of digits. If that sounds confusing, that is
because it is.

Note that MAQ is unmaintained and should not be used in new projects.

BWA’s colorspace support was dropped in versions more recent than 0.5.9, but
that version works well.

When you want to trim reads that will be mapped with BWA or MAQ, you can
use the `--bwa` option, which enables colorspace mode (`-c`), double-encoding
(`-d`), primer trimming (`-t`), all of which are required for BWA, in
addition to some other useful options.

The `--maq` option is an alias for `--bwa`.


Colorspace examples
-------------------

To cut an adapter from SOLiD data given in `solid.csfasta` and `solid.qual`,
to produce MAQ- and BWA-compatible output, allow the default of 10% errors
and write the resulting FASTQ file to output.fastq:

    cutadapt --bwa -a CGCCTTGGCCGTACAGCAG solid.csfasta solid.qual > output.fastq

Instead of redirecting standard output with `>`, the `-o` option can be used. This
also shows that you can give the adapter in colorspace and how to use a different
error rate:

    cutadapt --bwa -e 0.15 -a 330201030313112312 -o output.fastq solid.csfasta solid.qual

This does the same as above, but produces BFAST-compatible output, strips the _F3 suffix
from read names and adds the prefix "abc:" to them:

    cutadapt -c -e 0.15 -a 330201030313112312 -x abc: --strip-f3 solid.csfasta solid.qual > output.fastq


Bowtie
------

Quality values of colorspace reads are sometimes negative. Bowtie gets
confused and prints this message:
> Encountered a space parsing the quality string for read xyz
To avoid this problem, use the `--zero-cap` option (or the short
version `-z`), which converts negative quality values to zero. Since
BWA has a similar problem (it crashes) the option is automatically
enabled when `--bwa` is used.


Changes
=======

v1.5
----
* Cutadapt is again compatible with Python 3.
* U characters in the adapter sequence are automatically converted to T.

v1.4
----

* This release of cutadapt reduces the overhead of reading and writing files.
  On my test data set, a typical run of cutadapt (with a single adapter) takes
  40% less time due to the following two changes.
* Reading and writing of FASTQ files is faster (thanks to Cython).
* Reading and writing of gzipped files is faster (up to 2x) on systems
  where the `gzip` program is available.
* The quality trimming function is four times faster (also due to Cython).
* Fix the statistics output for 3' colorspace adapters: The reported lengths were one
  too short. Thanks to Frank Wessely for reporting this.
* Support the `--no-indels` option. This disallows insertions and deletions while
  aligning the adapter. Currently, the option is only available for anchored 5' adapters.
  This fixes issue 69.
* As a sideeffect of implementing the --no-indels option: For colorspace, the
  length of a read (for `--minimum- and --maximum-length`) is now computed after
  primer base removal (when `--trim-primer` is specified).
* Added one column to the info file that contains the name of the found adapter.
* Add an explanation about colorspace ambiguity to the README

v1.3
----

* Preliminary paired-end support with the --paired-output option (contributed by
  James Casbon). See the README section on how to use it.
* Improved statistics.
* Fix incorrectly reported amount of quality-trimmed Mbp (issue 57, fix by Chris Penkett)
* Add the `--too-long-output` option.
* Add the `--no-trim` option, contributed by Dave Lawrence.
* Port handwritten C alignment module to Cython.
* Fix the `--rest-file` option (issue 56)
* Slightly speed up alignment of 5' adapters.
* Support bzip2-compressed files.

v1.2
----
* At least 25% faster processing of .csfasta/.qual files due to faster parser.
* Between 10% and 30% faster writing of gzip-compressed output files.
* Support 5' adapters in color space, even when no primer trimming is requested.
* Add the `--info-file' option, which has a line for each found adapter.
* Named adapters are possible. Usage: `-a My_Adapter=ACCGTA` assigns the name "My_adapter".
* Improve alignment algorithm for better poly-A trimming when there are sequencing errors.
  Previously, not the longest possible poly-A tail would be trimmed.
* James Casbon contributed the `--discard-untrimmed` option.

v1.1
----
* Allow to "anchor" 5' adapters (`-g`), forcing them to be a prefix of the read.
  To use this, add the special character `^` to the beginning of the adapter sequence.
* Add the "-N" option, which allows 'N' characters within adapters to match literally.
* Speedup of approx. 25% when reading from .gz files and using Python 2.7.
* Allow to only trim qualities when no adapter is given on the command-line.
* Add a patch by James Casbon:
    * include read names (ids) in rest file
* Use nosetest for testing. To run, install nose and run "nosetests".
* When using cutadapt without installing it, you now need to run `bin/cutadapt` due to
  a new directory layout.
* Allow to give a colorspace adapter in basespace (gets automatically converted).
* Allow to search for 5' adapters (those specified with `-g`) in colorspace.
* Speed up the alignment by a factor of at least 3 by using Ukkonen's algorithm.
  The total runtime decreases by about 30% in the tested cases.
* allow to deal with colorspace FASTQ files from the SRA that contain a fake
  additional quality in the beginning (use --format sra-fastq)

v1.0
----
* ASCII-encoded quality values were assumed to be encoded as ascii(quality+33).
  With the new parameter `--quality-base`, this can be changed to ascii(quality+64),
  as used in some versions of the Illumina pipeline. (Fixes issue 7.)
* Allow to specify that adapters were ligated to the 5' end of reads. This change
  is based on a patch contributed by James Casbon.
* Due to cutadapt being published in EMBnet.journal, I found it appropriate
  to call this release version 1.0. Please see
  http://journal.embnet.org/index.php/embnetjournal/article/view/200 for the
  article and I would be glad if you cite it.
* Add Galaxy support, contributed by Lance Parsons.
* Patch by James Casbon: Allow N wildcards in read or adapter or both.
  Wildcard matching of 'N's in the adapter is always done. If 'N's within reads
  should also match without counting as error, this needs to be explicitly
  requested via --match-read-wildcards.

v0.9.5
------
* Fix issue 20: Make the report go to standard output when `-o`/`--output` is
  specified.
* Recognize .fq as an extension for FASTQ files
* many more unit tests
* The alignment algorithm has changed. It will now find some adapters that
  previously were missed. Note that this will produce different output than
  older cutadapt versions!

  Before this change, finding an adapter would work as follows:
  - Find an alignment between adapter and read -- longer alignments are
    better.
  - If the number of errors in the alignment (divided by length) is above the
    maximum error rate, report the adapter as not being found.
  Sometimes, the long alignment that is found had too many errors, but a
  shorter alignment would not. The adapter was then incorrectly seen as "not
  found". The new alignment algorithm checks the error rate while aligning and only
  reports alignments that do not have too many errors.

v0.9.4
------
* now compatible with Python 3
* Add the `--zero-cap` option, which changes negative quality values to zero.
  This is a workaround to avoid segmentation faults in BWA. The option is now
  enabled by default when `--bwa`/`--maq` is used.
* Lots of unit tests added. Run them with `cd tests && ./tests.sh` .
* Fix issue 16: `--discard-trimmed` did not work.
* Allow to override auto-detection of input file format with the new `-f`/`--format`
  parameter. This mostly fixes issue 12.
* Don't break when input file is empty.

v0.9.2
------
* Install a single 'cutadapt' Python package instead of multiple Python
  modules. This avoids cluttering the global namespace and should lead to less
  problems with other Python modules. Thanks to Steve Lianoglou for
  pointing this out to me!
* ignore case (ACGT vs acgt) when comparing the adapter with the read sequence
* .FASTA/.QUAL files (not necessarily color space) can now be read (some
  454 software uses this format)
* Move some functions into their own modules
* lots of refactoring: replace the fasta module with a much nicer seqio module.
* allow to input FASTA/FASTQ on standard input (also FASTA/FASTQ is
  autodetected)

v0.9
----
 * add `--too-short-output` and `--untrimmed-output`, based on patch by Paul Ryvkin (thanks!)
 * add `--maximum-length` parameter: discard reads longer than a specified length
 * group options by category in `--help` output
 * add `--length-tag` option. allows to fix read length in FASTA/Q comment lines
  (e.g., `length=123` becomes `length=58` after trimming) (requested by Paul Ryvkin)
 * add `-q`/`--quality-cutoff` option for trimming low-quality ends (uses the same algorithm
  as BWA)
 * some refactoring
 * the filename `-` is now interpreted as standard in or standard output

v0.8
----
* Change default behavior of searching for an adapter: The adapter is now assumed to
  be an adapter that has been ligated to the 3' end. This should be the correct behavior
  for at least the SOLiD small RNA protocol (SREK) and also for the Illumina protocol.
  To get the old behavior, which uses a heuristic to determine whether the adapter was
  ligated to the 5' or 3' end and then trimmed the read accordingly, use the new
  `-b` (`--anywhere`) option.
* Clear up how the statistics after processing all reads are printed.
* Fix incorrect statistics. Adapters starting at pos. 0 were correctly trimmed,
  but not counted.
* Modify scoring scheme: Improves trimming (some reads that should have been
  trimmed were not). Increases no. of trimmed reads in one of our SOLiD data sets
  from 36.5 to 37.6%.
* Speed improvements (20% less runtime on my test data set).

v0.7
----
* Useful exit codes
* Better error reporting when malformed files are encountered
* Add `--minimum-length` parameter for discarding reads that are shorter than
  a specified length after trimming.
* Generalize the alignment function a bit. This is preparation for
  supporting adapters that are specific to either the 5' or 3' end.
* pure Python fallback for alignment function for when the C module cannot
  be used.

v0.6
----
* Support gzipped input and output.
* Print timing information in statistics.

v0.5
----
* add --discard option which makes cutadapt discard reads in which an adapter occurs

v0.4
----
 * (more) correctly deal with multiple adapters: If a long adapter matches with lots of
   errors, then this could lead to a a shorter adapter matching with few errors getting ignored.

v0.3
----
 * fix huge memory usage (entire input file was unintentionally read into memory)

v0.2
----
* allow FASTQ input

v0.1
----
* initial release


To Do / Ideas
=============

  * show average error rate
  * In color space and probably also for Illumina data, gapped alignment
    is not necessary
  * use `str.format` instead of `%`
  * allow to change scores at runtime (using command-line parameters)
  * multi-threading
  * `--progress`
  * run pylint, pychecker
  * length histogram
  * refactor read_sequences (use classes)
  * put write_read into a Fast(a|q)Writer class?
  * allow .txt input/output
  * test on Windows
  * check whether input is FASTQ although -f fasta is given
  * close on StopIteration
  * search for adapters in the order in which they are given on the command line
  * more tests for the alignment algorithm
