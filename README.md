cutadapt
========

cutadapt removes adapter sequences from DNA high-throughput
sequencing data. This is usually necessary when the read length of the
machine is longer than the molecule that is sequenced, such as in
microRNA data.

cutadapt is implemented in Python, with an extension module,
written in C, that implements the alignment algorithm.


Project homepage
================

See http://code.google.com/p/cutadapt/ . Please use the Google code issue
tracker for bug reports and feature requests.


License
=======

(This is the MIT license.)

Copyright (c) 2010-2012 Marcel Martin <marcel.martin@tu-dortmund.de>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.


Dependencies
============

cutadapt needs Python 2.6 or later, including Python 3 (tested with
Python 3.1.2 and 3.2.2). Using Python 2.7 is recommended since it works fastest.
For installation from sources, a C compiler needs to be installed. The program
has been developed and tested on Ubuntu and OpenSuSE.


Installation
============

Replace "python" with "python3" in the following lines to install the Python 3
version.

    python setup.py build
    python setup.py install

If you get an error about a missing "Python.h" file, then make sure that
the python-dev package is installed (or python3-dev for Python 3).

Use without installation
========================

Build the C extension module (you can try to skip this step -- a compiled
version of the module for Linux x86 is already included):

    python setup.py build_ext -i

Then simply run the script from where it is, similar to this:

    bin/cutadapt --help

If you get any errors, first try to explicitly request a specific Python
version by running cutadapt like this:

    python2.7 bin/cutadapt --help


Galaxy
======

If you want to use cutadapt within the web-based Galaxy platform
(http://galaxy.psu.edu/), please see the README file in the galaxy/ subfolder.
Galaxy support was contributed by Lance Parsons.


How to use, examples
====================

Please also see the command-line help:
    cutadapt --help

The basic command-line for cutadapt looks like this:

    cutadapt -a AACCGGTT input.fastq > output.fastq

The adapter sequence is given with the `-a` option. Replace
AACCGGTT with your actual adapter sequence.

input.fastq is a file with reads. The result will be written
to standard output. Use redirection with `>` (or the `-o` option) to
write the output to a file.

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
(this is usually not needed), use the default error rate of 10%, write result to output.fa:

    cutadapt -n 3 -a TGAGACACGCAACAGGGGAAAGGCAAGGCACACAGGGGATAGG input.fa > output.fa


Multiple adapters
-----------------

As many adapters as desired can be given to the program by using the `-a`, `-b` or `-g`
in any combination, for example, five `-a` adapters and two `-g` adapters. All
adapters will be searched for, but only the best matching one will be trimmed
from each read (but see the `--times` option).

    cutadapt -b TGAGACACGCA -g AGGCACACAGGG input.fastq > output.fastq


Quality Trimming
----------------

The `-q` (or `--trim-qualities`) parameter can be used to trim low-quality ends
from reads before adapter removal. For this to work correctly, the quality
values must be encoded as ascii(phred quality + 33). If they are encoded as
ascii(phred quality + 64), you need to add `--quality-base=64` to the command line.

The trimming algorithm is the same as the one used by BWA. That is: Subtract
the given cutoff from all qualities; compute partial sums from all indices to
the end of the sequence; cut sequence at the index at which the sum is minimal.


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
algorithm is implemented as a Python extension module in `calignmodule.c`.

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
it are remved.

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
can be given both as colors or as nucleotides. In the latter case, they will
simply be converted automatically. For example, to trim an adapter from
`solid.csfasta` and `solid.qual`, use this command-line:

	cutadapt -c -a CGCCTTGGCCGTACAGCAG solid.csfasta solid.qual > output.fastq

In case you know the colorspace adapter sequence, you can also write `330201030313112312`
instead of `CGCCTTGGCCGTACAGCAG` and the result is the same.


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
(`-d`) and primer trimming (`-t`), all of which are required for BWA, in
addition to some other useful options.

There is also the `--maq` option, which is simply another name for the
`--bwa` option.


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


Changes
=======

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
  * show table of length vs. max errors
  * In color space and probably also for Illumina data, gapped alignment
    is not necessary
  * bzip2 support
  * use `str.format` instead of `%`
  * allow to change scores at runtime (using command-line parameters)
  * multi-threading
  * `--progress`
  * run pylint, pychecker
  * print adapter fragments in statistics
  * `--discard-uncut` (discard sequences in which no adapter was found)
  * no. of trimmed nucleotides
  * length histogram
  * refactor read_sequences (use classes)
  * put write_read into a Fast(a|q)Writer class?
  * allow .txt input/output
  * test on Windows
  * check whether input is FASTQ although -f fasta is given
  * close on StopIteration
  * search for adapters in the order in which they are given on the command line
  * re-write alignment in Cython

