=======
Changes
=======

v1.10
-----

* Added a new “linked adapter” type, which can be used to search for a 5' and a
  3' adapter at the same time. Use ``-a ADAPTER1...ADAPTER2` to search
  for a linked adapter. ADAPTER1 is interpreted as an anchored 5' adapter, which
  is searched for first. Only if ADAPTER1 is found will ADAPTER2 be searched
  for, which is a regular 3' adapter.
* Added experimental ``--nextseq-trim`` option for quality trimming of NextSeq
  data. This is necessary because that machine cannot distinguish between G and
  reaching the end of the fragment (it encodes G as 'black').
* Even when trimming FASTQ files, output can now be FASTA (quality values are
  simply dropped). Use the ``-o``/``-p`` options with a file name that ends in
  ``.fasta`` or ``.fa`` to enable this.
* Cutadapt does not bundle pre-compiled C extension modules (``.so`` files)
  anymore. This affects only users that run cutadapt directly from an unpacked
  tarball. Install through ``pip`` or ``conda`` instead.
* Fix issue #167: Option ``--quiet`` was not entirely quiet.
* Fix issue #199: Be less strict when checking for properly-paired reads.
* This is the last version of cutadapt to support Python 2.6. Future versions
  will require at least Python 2.7.

v1.9.1
------

* Added ``--pair-filter`` option, which :ref:`modifies how filtering criteria
  apply to paired-end reads <filtering-paired>`
* Add ``--too-short-paired-output`` and ``--too-long-paired-output`` options.
* Fix incorrect number of trimmed bases reported if ``--times`` option was used.

v1.9
----

* Indels in the alignment can now be disabled for all adapter types (use
  ``--no-indels``).
* Quality values are now printed in the info file (``--info-file``)
  when trimming FASTQ files. Fixes issue #144.
* Options ``--prefix`` and ``--suffix``, which modify read names, now accept the
  placeholder ``{name}`` and will replace it with the name of the found adapter.
  Fixes issue #104.
* Interleaved FASTQ files: With the ``--interleaved`` switch, paired-end reads
  will be read from and written to interleaved FASTQ files. Fixes issue #113.
* Anchored 5' adapters can now be specified by writing ``-a SEQUENCE...`` (note
  the three dots).
* Fix ``--discard-untrimmed`` and ``--discard-trimmed`` not working as expected
  in paired-end mode (issue #146).
* The minimum overlap is now automatically reduced to the adapter length if it
  is too large. Fixes part of issue #153.
* Thanks to Wolfgang Gerlach, there is now a Dockerfile.
* The new ``--debug`` switch makes cutadapt print out the alignment matrix.

v1.8.3
------

* Fix issue #95: Untrimmed reads were not listed in the info file.
* Fix issue #138: pip install cutadapt did not work with new setuptools versions.
* Fix issue #137: Avoid a hang when writing to two or more gzip-compressed
  output files in Python 2.6.

v1.8.1
------

* Fix #110: Counts for 'too short' and 'too long' reads were swapped in statistics.
* Fix #115: Make ``--trim-n`` work also on second read for paired-end data.

v1.8
----

* Support single-pass paired-end trimming with the new ``-A``/``-G``/``-B``/``-U``
  parameters. These work just like their -a/-g/-b/-u counterparts, but they
  specify sequences that are removed from the *second read* in a pair.

  Also, if you start using one of those options, the read modification options
  such as ``-q`` (quality trimming) are applied to *both* reads. For backwards
  compatibility, read modifications are applied to the first read only if
  neither of ``-A``/``-G``/``-B``/``-U`` is used. See `the
  documentation <http://cutadapt.readthedocs.org/en/latest/guide.html#paired-end>`_
  for details.

  This feature has not been extensively tested, so please give feedback if
  something does not work.
* The report output has been re-worked in order to accomodate the new paired-end
  trimming mode. This also changes the way the report looks like in single-end
  mode. It is hopefully now more accessible.
* Chris Mitchell contributed a patch adding two new options: ``--trim-n``
  removes any ``N`` bases from the read ends, and the ``--max-n`` option can be
  used to filter out reads with too many ``N``.
* Support notation for repeated bases in the adapter sequence: Write ``A{10}``
  instead of ``AAAAAAAAAA``. Useful for poly-A trimming: Use ``-a A{100}`` to
  get the longest possible tail.
* Quality trimming at the 5' end of reads is now supported. Use ``-q 15,10`` to
  trim the 5' end with a cutoff of 15 and the 3' end with a cutoff of 10.
* Fix incorrectly reported statistics (> 100% trimmed bases) when ``--times``
  set to a value greater than one.
* Support .xz-compressed files (if running in Python 3.3 or later).
* Started to use the GitHub issue tracker instead of Google Code. All old issues
  have been moved.

v1.7
----
* IUPAC characters are now supported. For example, use ``-a YACGT`` for an
  adapter that matches both ``CACGT`` and ``TACGT`` with zero errors. Disable
  with ``-N``. By default, IUPAC characters in the read are not interpreted in
  order to avoid matches in reads that consist of many (low-quality) ``N``
  bases. Use ``--match-read-wildcards`` to enable them also in the read.
* Support for demultiplexing was added. This means that reads can be written to
  different files depending on which adapter was found. See `the section in the
  documentation <http://cutadapt.readthedocs.org/en/latest/guide.html#demultiplexing>`_
  for how to use it. This is currently only supported for single-end reads.
* Add support for anchored 3' adapters. Append ``$`` to the adapter sequence to
  force the adapter to appear in the end of the read (as a suffix). Closes
  issue #81.
* Option ``--cut`` (``-u``) can now be specified twice, once for each end of the
  read. Thanks to Rasmus Borup Hansen for the patch!
* Options ``--minimum-length``/``--maximum-length`` (``-m``/``-M``) can be used
  standalone. That is, cutadapt can be used to filter reads by length without
  trimming adapters.
* Fix bug: Adapters read from a FASTA file can now be anchored.

v1.6
----
* Fix bug: Ensure ``--format=...`` can be used even with paired-end input.
* Fix bug: Sometimes output files would be incomplete because they were not
  closed correctly.
* Alignment algorithm is a tiny bit faster.
* Extensive work on the documentation. It's now available at
  https://cutadapt.readthedocs.org/ .
* For 3' adapters, statistics about the bases preceding the trimmed adapter
  are collected and printed. If one of the bases is overrepresented, a warning
  is shown since this points to an incomplete adapter sequence. This happens,
  for example, when a TruSeq adapter is used but the A overhang is not taken
  into account when running cutadapt.
* Due to code cleanup, there is a change in behavior: If you use
  ``--discard-trimmed`` or ``--discard-untrimmed`` in combination with
  ``--too-short-output`` or ``--too-long-output``, then cutadapt now writes also
  the discarded reads to the output files given by the ``--too-short`` or
  ``--too-long`` options. If anyone complains, I will consider reverting this.
* Galaxy support files are now in `a separate
  repository <https://bitbucket.org/lance_parsons/cutadapt_galaxy_wrapper>`_.

v1.5
----
* Adapter sequences can now be read from a FASTA file. For example, write
  ``-a file:adapters.fasta`` to read 3' adapters from ``adapters.fasta``. This works
  also for ``-b`` and ``-g``.
* Add the option ``--mask-adapter``, which can be used to not remove adapters,
  but to instead mask them with ``N`` characters. Thanks to Vittorio Zamboni
  for contributing this feature!
* U characters in the adapter sequence are automatically converted to T.
* Do not run Cython at installation time unless the --cython option is provided.
* Add the option -u/--cut, which can be used to unconditionally remove a number
  of bases from the beginning or end of each read.
* Make ``--zero-cap`` the default for colorspace reads.
* When the new option ``--quiet`` is used, no report is printed after all reads
  have been processed.
* When processing paired-end reads, cutadapt now checks whether the reads are
  properly paired.
* To properly handle paired-end reads, an option --untrimmed-paired-output was
  added.

v1.4
----
* This release of cutadapt reduces the overhead of reading and writing files.
  On my test data set, a typical run of cutadapt (with a single adapter) takes
  40% less time due to the following two changes.
* Reading and writing of FASTQ files is faster (thanks to Cython).
* Reading and writing of gzipped files is faster (up to 2x) on systems
  where the ``gzip`` program is available.
* The quality trimming function is four times faster (also due to Cython).
* Fix the statistics output for 3' colorspace adapters: The reported lengths were one
  too short. Thanks to Frank Wessely for reporting this.
* Support the ``--no-indels`` option. This disallows insertions and deletions while
  aligning the adapter. Currently, the option is only available for anchored 5' adapters.
  This fixes issue 69.
* As a sideeffect of implementing the --no-indels option: For colorspace, the
  length of a read (for ``--minimum-length`` and ``--maximum-length``) is now computed after
  primer base removal (when ``--trim-primer`` is specified).
* Added one column to the info file that contains the name of the found adapter.
* Add an explanation about colorspace ambiguity to the README

v1.3
----
* Preliminary paired-end support with the ``--paired-output`` option (contributed by
  James Casbon). See the README section on how to use it.
* Improved statistics.
* Fix incorrectly reported amount of quality-trimmed Mbp (issue 57, fix by Chris Penkett)
* Add the ``--too-long-output`` option.
* Add the ``--no-trim`` option, contributed by Dave Lawrence.
* Port handwritten C alignment module to Cython.
* Fix the ``--rest-file`` option (issue 56)
* Slightly speed up alignment of 5' adapters.
* Support bzip2-compressed files.

v1.2
----
* At least 25% faster processing of .csfasta/.qual files due to faster parser.
* Between 10% and 30% faster writing of gzip-compressed output files.
* Support 5' adapters in colorspace, even when no primer trimming is requested.
* Add the ``--info-file`` option, which has a line for each found adapter.
* Named adapters are possible. Usage: ``-a My_Adapter=ACCGTA`` assigns the name "My_adapter".
* Improve alignment algorithm for better poly-A trimming when there are sequencing errors.
  Previously, not the longest possible poly-A tail would be trimmed.
* James Casbon contributed the ``--discard-untrimmed`` option.

v1.1
----
* Allow to "anchor" 5' adapters (``-g``), forcing them to be a prefix of the read.
  To use this, add the special character ``^`` to the beginning of the adapter sequence.
* Add the "-N" option, which allows 'N' characters within adapters to match literally.
* Speedup of approx. 25% when reading from .gz files and using Python 2.7.
* Allow to only trim qualities when no adapter is given on the command-line.
* Add a patch by James Casbon: include read names (ids) in rest file
* Use nosetest for testing. To run, install nose and run "nosetests".
* When using cutadapt without installing it, you now need to run ``bin/cutadapt`` due to
  a new directory layout.
* Allow to give a colorspace adapter in basespace (gets automatically converted).
* Allow to search for 5' adapters (those specified with ``-g``) in colorspace.
* Speed up the alignment by a factor of at least 3 by using Ukkonen's algorithm.
  The total runtime decreases by about 30% in the tested cases.
* allow to deal with colorspace FASTQ files from the SRA that contain a fake
  additional quality in the beginning (use ``--format sra-fastq``)

v1.0
----
* ASCII-encoded quality values were assumed to be encoded as ascii(quality+33).
  With the new parameter ``--quality-base``, this can be changed to ascii(quality+64),
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
  requested via ``--match-read-wildcards``.

v0.9.5
------
* Fix issue 20: Make the report go to standard output when ``-o``/``--output`` is
  specified.
* Recognize `.fq` as an extension for FASTQ files
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
* Add the ``--zero-cap`` option, which changes negative quality values to zero.
  This is a workaround to avoid segmentation faults in BWA. The option is now
  enabled by default when ``--bwa``/``--maq`` is used.
* Lots of unit tests added. Run them with ``cd tests && ./tests.sh``.
* Fix issue 16: ``--discard-trimmed`` did not work.
* Allow to override auto-detection of input file format with the new ``-f``/``--format``
  parameter. This mostly fixes issue 12.
* Don't break when input file is empty.

v0.9.2
------
* Install a single ``cutadapt`` Python package instead of multiple Python
  modules. This avoids cluttering the global namespace and should lead to less
  problems with other Python modules. Thanks to Steve Lianoglou for
  pointing this out to me!
* ignore case (ACGT vs acgt) when comparing the adapter with the read sequence
* .FASTA/.QUAL files (not necessarily colorspace) can now be read (some
  454 software uses this format)
* Move some functions into their own modules
* lots of refactoring: replace the fasta module with a much nicer seqio module.
* allow to input FASTA/FASTQ on standard input (also FASTA/FASTQ is
  autodetected)

v0.9
----
* add ``--too-short-output`` and ``--untrimmed-output``, based on patch by Paul Ryvkin (thanks!)
* add ``--maximum-length`` parameter: discard reads longer than a specified length
* group options by category in ``--help`` output
* add ``--length-tag`` option. allows to fix read length in FASTA/Q comment lines
  (e.g., ``length=123`` becomes ``length=58`` after trimming) (requested by Paul Ryvkin)
* add ``-q``/``--quality-cutoff`` option for trimming low-quality ends (uses the same algorithm
  as BWA)
* some refactoring
* the filename ``-`` is now interpreted as standard in or standard output

v0.8
----
* Change default behavior of searching for an adapter: The adapter is now assumed to
  be an adapter that has been ligated to the 3' end. This should be the correct behavior
  for at least the SOLiD small RNA protocol (SREK) and also for the Illumina protocol.
  To get the old behavior, which uses a heuristic to determine whether the adapter was
  ligated to the 5' or 3' end and then trimmed the read accordingly, use the new
  ``-b`` (``--anywhere``) option.
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
* Add ``--minimum-length`` parameter for discarding reads that are shorter than
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
* add ``--discard`` option which makes cutadapt discard reads in which an adapter occurs

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
