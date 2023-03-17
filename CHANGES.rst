
=========
Changelog
=========

v4.3 (2023-03-17)
-----------------

* :pr:`663`: Cutadapt became significantly faster due to an added runtime
  heuristic that avoids running the full alignment algorithm if it can be
  proven that it cannot succeed. Thanks to @rhpvorderman for this great
  improvement!
* :issue:`665`: 5' adapters did not allow partial matches in the beginning
  when the :ref:`rightmost <rightmost>` adapter-search parameter was used.
* :issue:`662`: Fixed assertion error when ``--discard-untrimmed`` was used
  together with ``--json`` and demultiplexing.
* :issue:`674`: When reading 3' adapters from an external file, they can now
  all be anchored by using the syntax ``-a file$:adapters.fasta`` (note the
  ``$`` in ``file$:``).
* :issue:`669`: The ``--rename`` option now understands the ``\t`` escape
  sequence and will insert a tab character in its place. This is useful when
  transferring FASTQ header comments to SAM tags.

v4.2 (2022-12-09)
-----------------

* :issue:`654`: When determining the error rate for a partial match of an
  adapter with ``N`` wildcards, the number of non-N bases was not computed
  correctly, which could lead to matches not being found.
* :issue:`546`: Automatically replace ``I`` in adapter sequences with ``N``.
  ``I`` is used to encode inosine, which matches any base. Contributed by @peterjc.
* :issue:`528`: Cutadapt should now no long hang in multicore mode when an error
  was raised in a worker process (for example, when an incorrectly formatted
  FASTQ file was encountered).

v4.1 (2022-06-07)
-----------------

* :issue:`624`: You can now combine reading adapter sequences from an external file
  with the search parameter notation. For example,
  ``-a "file:adapters.fasta;min_overlap=5"`` sets the minimum overlap to 5 for all
  adapters in ``adapters.fasta``.
* :issue:`361`: When reading 5' adapters from an external file, they can now
  all be anchored by using the syntax ``-g ^file:adapters.fasta``
  (note the ``^`` before ``file:``).
* :issue:`254`: Finding the *rightmost* 5' adapter occurrence is now supported by using the
  ``rightmost`` search parameter (the default is to find the leftmost occurrence).
* :issue:`615`: Fix linked adapter statistics for 5' and 3' end not
  being reported separated correctly.
* :issue:`616`: Report correct number of quality-trimmed bases when
  both ``-q`` and ``--nextseq-trim`` are used.

v4.0 (2022-04-13)
-----------------

* :issue:`604`, :pr:`608`: The :ref:`alignment algorithm was tweaked <algorithm-indel-scores>`
  to penalize indels more and to more accurately pick the leftmost adapter
  occurrence if there are multiple. This will normally affect very few
  reads, but should generally lead to fewer surprising results in cases
  where it matters. Because this changes trimming results, it was appropriate
  to bump the major version to 4.
* :issue:`607`: Print an error when an output file was specified
  multiple times (for example, for ``--untrimmed-output`` and
  ``--too-short-output``). Sending output from different filters to
  the same file is not supported at the moment.
* :issue:`603`: When ``-e`` was used with an absolute number of errors
  and there were ``N`` wildcards in the sequence, the actual number of
  allowed errors was too low.
* Speed up quality trimming (both ``-q`` and ``--nextseq-trim``) somewhat.
* Python 3.6 is no longer supported as it is end-of-life.

v3.7 (2022-02-23)
-----------------

* :issue:`600`: Fixed ``{match_sequence}`` placeholder not working when
  renaming paired-end reads.

v3.6 (2022-02-18)
---------------------

* :issue:`437`: Add ``{match_sequence}`` to the placeholders that ``--rename``
  accepts. This allows to add the sequence matching an adapter (including
  errors) to the read header. An empty string is inserted if there is no match.
* :issue:`589`: Windows wheels are now available on PyPI. That is,
  ``pip install`` will no longer attempt to compile things, but just install
  a pre-compiled binary.
* :issue:`592`: Clarify in documentation and error messages that anchored
  adapters need to match in full and that therefore setting an explict
  minimum overlap (``min_overlap=``, ``o=``) for them is not possible.

v3.5 (2021-09-29)
-----------------

* :issue:`555`: Add support for dumping statistics in JSON format using ``--json``.
* :issue:`541`: Add a "Read fate breakdown" section heading to the report, and also
  add statistics for reads discarded because of ``--discard-untrimmed`` and
  ``--discard-trimmed``. With this, the numbers in that section should add up to 100%.
* Add option ``-Q``, which allows to specify a quality-trimming threshold for R2 that is
  different from the one for R1.
* :issue:`567`: Add ``noindels`` adapter-trimming parameter. You can now write
  ``-a "ADAPTER;noindels"`` to disallow indels for a single adapter only.
* :issue:`570`: Fix ``--pair-adapters`` not finding some pairs when reads contain
  more than one adapter.
* :issue:`524`: Fix a memory leak when using ``--info-file`` with multiple cores.
* :issue:`559`: Fix adjacent base statistics not being shown for linked adapters.

v3.4 (2021-03-30)
-----------------

* :issue:`481`: An experimental single-file Windows executable of Cutadapt
  is `available for download on the GitHub "releases"
  page <https://github.com/marcelm/cutadapt/releases>`_.
* :issue:`517`: Report correct sequence in info file if read was reverse complemented
* :issue:`517`: Added a column to the info file that shows whether the read was
  reverse-complemented (if ``--revcomp`` was used)
* :issue:`320`: Fix (again) "Too many open files" when demultiplexing

v3.3 (2021-03-04)
-----------------

* :issue:`504`: Fix a crash on Windows.
* :issue:`490`: When ``--rename`` is used with ``--revcomp``, disable adding the
  ``rc`` suffix to reads that were reverse-complemented.
* Also, there is now a ``{rc}`` template variable for the ``--rename`` option, which
  is replaced with "rc" if the read was reverse-complemented (and the empty string if not).
* :issue:`512`: Fix issue :issue:`128` once more (the “Reads written” figure in the report
  incorrectly included both trimmed and untrimmed reads if ``--untrimmed-output`` was used).
* :issue:`515`: The report is now sent to stderr if any output file is
  written to stdout

v3.2 (2021-01-07)
-----------------

* :issue:`437`: Implement a ``--rename`` option for :ref:`flexible read
  name modifications <read-renaming>` such as moving a barcode sequence
  into the read name.
* :issue:`503`: The index for demultiplexing is now created a lot faster
  (within seconds instead of minutes) when allowing indels.
* :issue:`499`: Fix combinatorial demultiplexing not working when using
  multiple cores.

v3.1 (2020-12-03)
-----------------

* :issue:`443`: With ``--action=retain``, it is now possible to trim reads while
  leaving the adapter sequence itself in the read. That is, only the sequence
  before (for 5’ adapters) or after (for 3’ adapters) is removed. With linked
  adapters, both adapters are retained.
* :issue:`495`: Running with multiple cores did not work using macOS and Python 3.8+.
  To prevent problems like these in the future, automated testing has been extended
  to also run on macOS.
* :issue:`482`: Print statistics for ``--discard-casava`` and ``--max-ee`` in the
  report.
* :issue:`497`: The changelog for 3.0 previously forgot to mention that the following
  options, which were deprecated in version 2.0, have now been removed, and
  using them will lead to an error: ``--format``, ``--colorspace``, ``-c``, ``-d``,
  ``--double-encode``, ``-t``, ``--trim-primer``, ``--strip-f3``, ``--maq``,
  ``--bwa``, ``--no-zero-cap``. This frees up some single-character options,
  allowing them to be re-purposed for future Cutadapt features.

v3.0 (2020-11-10)
-----------------

* Demultiplexing on multiple cores is now supported. This was the last feature that
  only ran single-threaded.
* :issue:`478`: Demultiplexing now always generates all possible output files.
* :issue:`358`: You can now use ``-e`` also :ref:`to specify the maximum number of
  errors <error-tolerance>` (instead of the maximum error rate). For example, write
  ``-e 2`` to allow two errors over a full-length adapter match.
* :pr:`486`: Trimming many anchored adapters (for example when demultiplexing)
  is now faster by using an index even when indels are allowed. Previously, Cutadapt
  would only be able to build an index with ``--no-indels``.
* :issue:`469`: Cutadapt did not run under Python 3.8 on recent macOS versions.
* :issue:`425`: Change the default compression level for ``.gz`` output files from 6
  to 5. This reduces the time used for compression by about 50% while increasing file
  size by less than 10%. To get the old behavior, use ``--compression-level=6``.
  If you use Cutadapt to create intermediate files that are deleted anyway,
  consider also using the even faster option ``-Z`` (same as ``--compression-level=1``).
* :pr:`485`: Fix that, under some circumstances, in particular when trimming a
  5' adapter and there was a mismatch in its last nucleotide(s), not the entire adapter
  sequence would be trimmed from the read. Since fixing this required changed the
  alignment algorithm slightly, this is a backwards incompatible change.
* Fix that the report did not include the number of reads that are too long, too short
  or had too many ``N``. (This unintentionally disappeared in a previous version.)
* :issue:`487`: When demultiplexing, the reported number of written pairs was
  always zero.
* :issue:`497`: The following options, which were deprecated in version 2.0, have
  been removed, and using them will lead to an error:
  ``--format``, ``--colorspace``, ``-c``, ``-d``, ``--double-encode``,
  ``-t``, ``--trim-primer``, ``--strip-f3``, ``--maq``, ``--bwa``, ``--no-zero-cap``.
  This frees up some single-character options,
  allowing them to be re-purposed for future Cutadapt features.
* Ensure Cutadapt runs under Python 3.9.
* Drop support for Python 3.5.

v2.10 (2020-04-22)
------------------

* Fixed a performance regression introduced in version 2.9.
* :pr:`449`: ``--action=`` could not be used with ``--pair-adapters``.
  Fix contributed by wlokhorst.
* :issue:`450`: ``--untrimmed-output``, ``--too-short-output`` and ``--too-long-output`` can
  now be written interleaved.
* :issue:`453`: Fix problem that ``N`` wildcards in adapters did not match ``N`` characters
  in the read. ``N`` characters now match any character in the read, independent of whether
  ``--match-read-wildcards`` is used or not.
* With ``--action=lowercase``/``mask``, print which sequences would have been
  removed in the “Overview of removed sequences” statistics. Previously, it
  would show that no sequences have been removed.

v2.9 (2020-03-18)
-----------------

* :issue:`441`: Add a ``--max-ee`` (or ``--max-expected-errors``) option
  for filtering reads whose number of expected errors exceeds the given
  threshold. The idea comes from
  `Edgar et al. (2015) <https://academic.oup.com/bioinformatics/article/31/21/3476/194979>`_.
* :issue:`438`: The info file now contains the `` rc`` suffix that is added to
  the names of reverse-complemented reads (with ``--revcomp``).
* :issue:`448`: ``.bz2`` and ``.xz`` output wasn’t possible in multi-core mode.

v2.8 (2020-01-13)
-----------------

* :issue:`220`: With option ``--revcomp``, Cutadapt now searches both the read
  and its reverse complement for adapters. The version that matches best is
  kept. This can be used to “normalize” strandedness.
* :issue:`430`: ``--action=lowercase`` now works with linked adapters
* :issue:`431`: Info files can now be written even for linked adapters.

v2.7 (2019-11-22)
-----------------

* :issue:`427`: Multicore is now supported even when using ``--info-file``,
  ``--rest-file`` or ``--wildcard-file``. The only remaining feature that
  still does not work with multicore is now demultiplexing.
* :issue:`290`: When running on a single core, Cutadapt no longer spawns
  external ``pigz`` processes for writing gzip-compressed files. This is a first
  step towards ensuring that using ``--cores=n`` uses only at most *n* CPU
  cores.
* This release adds support for Python 3.8.

v2.6 (2019-10-26)
-----------------

* :issue:`395`: Do not show animated progress when ``--quiet`` is used.
* :issue:`399`: When two adapters align to a read equally well (in terms
  of the number of matches), prefer the alignment that has fewer errors.
* :issue:`401` Give priority to adapters given earlier on the command
  line. Previously, the priority was: All 3' adapters, all 5' adapters,
  all anywhere adapters. In rare cases this could lead to different results.
* :issue:`404`: Fix an issue preventing Cutadapt from being used on Windows.
* This release no longer supports Python 3.4 (which has reached end of life).


v2.5 (2019-09-04)
-----------------

* :issue:`391`: Multicore is now supported even when using
  ``--untrimmed-output``, ``--too-short-output``, ``--too-long-output``
  or the corresponding ``...-paired-output`` options.
* :issue:`393`: Using ``--info-file`` no longer crashes when processing
  paired-end data. However, the info file itself will only contain results
  for R1.
* :issue:`394`: Options ``-e``/``--no-indels``/``-O`` were ignored for
  linked adapters
* :issue:`320`: When a “Too many open files” error occurs during
  demultiplexing, Cutadapt can now automatically raise the limit and
  re-try if the limit is a “soft” limit.


v2.4 (2019-07-09)
-----------------

* :issue:`292`: Implement support for demultiplexing paired-end reads that use
  :ref:`combinatorial indexing (“combinatorial demultiplexing”)
  <combinatorial-demultiplexing>`.
* :pr:`384`: Speed up reading compressed files by requiring an xopen version
  that uses an external pigz process even for *reading* compressed input files
  (not only for writing).
* :issue:`381`: Fix ``--report=minimal`` not working.
* :issue:`380`: Add a ``--fasta`` option for forcing that FASTA is written
  to standard output even when input is FASTQ. Previously, forcing
  FASTA was only possible by providing an output file name.


v2.3 (2019-04-25)
-----------------

* :issue:`378`: The ``--pair-adapters`` option, added in version 2.1, was
  not actually usable for demultiplexing.


v2.2 (2019-04-20)
---------------------

* :issue:`376`: Fix a crash when using anchored 5' adapters together with
  ``--no-indels`` and trying to trim an empty read.
* :issue:`369`: Fix a crash when attempting to trim an empty read using a ``-g``
  adapter with wildcards.

v2.1 (2019-03-15)
-----------------

* :issue:`366`: Fix problems when combining ``--cores`` with
  reading from standard input or writing to standard output.
* :issue:`347`: Support :ref:`“paired adapters” <paired-adapters>`. One use case is
  demultiplexing Illumina *Unique Dual Indices* (UDI).

v2.0 (2019-03-06)
-----------------

This is a major new release with lots of bug fixes and new features, but
also some backwards-incompatible changes. These should hopefully
not affect too many users, but please make sure to review them and
possibly update your scripts!

Backwards-incompatible changes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* :issue:`329`: Linked adapters specified with ``-a ADAPTER1...ADAPTER2``
  are no longer anchored by default. To get results consist with the old
  behavior, use ``-a ^ADAPTER1...ADAPTER2`` instead.
* Support for colorspace data was removed. Thus, the following command-line
  options can no longer be used: ``-c``, ``-d``, ``-t``, ``--strip-f3``,
  ``--maq``, ``--bwa``, ``--no-zero-cap``.
* “Legacy mode” has been removed. This mode was enabled under certain
  conditions and would change the behavior such that the read-modifying options
  such as ``-q`` would only apply to the forward/R1 reads. This was necessary
  for compatibility with old Cutadapt versions, but became increasingly
  confusing.
* :issue:`360`: Computation of the error rate of an adapter match no longer
  counts the ``N`` wildcard bases. Previously, an adapter like ``N{18}CC``
  (18 ``N`` wildcards followed by ``CC``) would effectively match
  anywhere because the default error rate of 0.1 (10%) would allow for
  two errors. The error rate of a match is now computed as
  the number of non-``N`` bases in the matching part of the adapter
  divided by the number of errors.
* This release of Cutadapt requires at least Python 3.4 to run. Python 2.7
  is no longer supported.

Features
~~~~~~~~

* A progress indicator is printed while Cutadapt is working. If you redirect
  standard error to a file, the indicator is disabled.
* Reading of FASTQ files has gotten faster due to a new parser. The FASTA
  and FASTQ reading/writing functions are now available as part of the
  `dnaio library <https://github.com/marcelm/dnaio/>`_. This is a separate
  Python package that can be installed independently from Cutadapt.
  There is one regression at the moment: FASTQ files that use a second
  header (after the "+") will have that header removed in the output.
* Some other performance optimizations were made. Speedups of up to 15%
  are possible.
* Demultiplexing has become a lot faster :ref:`under certain conditions <speed-up-demultiplexing>`.
* :issue:`335`: For linked adapters, it is now possible to
  :ref:`specify which of the two adapters should be required <linked-override>`,
  overriding the default.
* :issue:`166`: By specifying ``--action=lowercase``, it is now possible
  to not trim adapters, but to instead convert the section of the read
  that would have been trimmed to lowercase.

Bug fixes
~~~~~~~~~

* Removal of legacy mode fixes also :issue:`345`: ``--length`` would not enable
  legacy mode.
* The switch to ``dnaio`` also fixed :issue:`275`: Input files with
  non-standard names now no longer lead to a crash. Instead the format
  is now recognized from the file content.
* Fix :issue:`354`: Sequences given using ``file:`` can now be unnamed.
* Fix :issue:`257` and :issue:`242`: When only R1 or only R2 adapters are given, the
  ``--pair-filter`` setting is now forced to ``both`` for the
  ``--discard-untrimmed`` (and ``--untrimmed-(paired-)output``) filters.
  Otherwise, with the default ``--pair-filter=any``, all pairs would be
  considered untrimmed because one of the reads in the pair is always
  untrimmed.

Other
~~~~~

* :issue:`359`: The ``-f``/``--format`` option is now ignored and a warning
  will be printed if it is used. The input file format is always
  auto-detected.


v1.18 (2018-09-07)
------------------

Features
~~~~~~~~

* Close :issue:`327`: Maximum and minimum lengths can now be specified
  separately for R1 and R2 with ``-m LENGTH1:LENGTH2``. One of the
  lengths can be omitted, in which case only the length of the other
  read is checked (as in ``-m 17:`` or ``-m :17``).
* Close :issue:`322`: Use ``-j 0`` to auto-detect how many cores to run on.
  This should even work correctly on cluster systems when Cutadapt runs as
  a batch job to which fewer cores than exist on the machine have been
  assigned. Note that the number of threads used by ``pigz`` cannot be
  controlled at the moment, see :issue:`290`.
* Close :issue:`225`: Allow setting the maximum error rate and minimum overlap
  length per adapter. A new :ref:`syntax for adapter-specific
  parameters <trimming-parameters>` was added for this. Example:
  ``-a "ADAPTER;min_overlap=5"``.
* Close :issue:`152`: Using the new syntax for adapter-specific parameters,
  it is now possible to allow partial matches of a 3' adapter at the 5' end
  (and partial matches of a 5' adapter at the 3' end) by specifying the
  ``anywhere`` parameter (as in ``-a "ADAPTER;anywhere"``).
* Allow ``--pair-filter=first`` in addition to ``both`` and ``any``. If
  used, a read pair is discarded if the filtering criterion applies to R1;
  and R2 is ignored.
* Close :issue:`112`: Implement a ``--report=minimal`` option for printing
  a succinct two-line report in tab-separated value (tsv) format. Thanks
  to :user:`jvolkening` for coming up with an initial patch!

Bug fixes
~~~~~~~~~

* Fix :issue:`128`: The “Reads written” figure in the report incorrectly
  included both trimmed and untrimmed reads if ``--untrimmed-output`` was used.

Other
~~~~~

* The options ``--no-trim`` and ``--mask-adapter`` should now be written as
  ``--action=mask`` and ``--action=none``. The old options still work.
* This is the last release to support `colorspace data <https://cutadapt.readthedocs.io/en/v1.18/colorspace.html>`_
* This is the last release to support Python 2.


v1.17 (2018-08-20)
------------------

* Close :issue:`53`: Implement adapters :ref:`that disallow internal matches <non-internal>`.
  This is a bit like anchoring, but less strict: The adapter sequence
  can appear at different lengths, but must always be at one of the ends.
  Use ``-a ADAPTERX`` (with a literal ``X``) to disallow internal matches
  for a 3' adapter. Use ``-g XADAPTER`` to disallow for a 5' adapter.
* :user:`klugem` contributed PR :issue:`299`: The ``--length`` option (and its
  alias ``-l``) can now be used with negative lengths, which will remove bases
  from the beginning of the read instead of from the end.
* Close :issue:`107`: Add a ``--discard-casava`` option to remove reads
  that did not pass CASAVA filtering (this is possibly relevant only for
  older datasets).
* Fix :issue:`318`: Cutadapt should now be installable with Python 3.7.
* Running Cutadapt under Python 3.3 is no longer supported (Python 2.7 or
  3.4+ are needed)
* Planned change: One of the next Cutadapt versions will drop support for
  Python 2 entirely, requiring Python 3.

v1.16 (2018-02-21)
------------------

* Fix :issue:`291`: When processing paired-end reads with multiple cores, there
  could be errors about incomplete FASTQs although the files are intact.
* Fix :issue:`280`: Quality trimming statistics incorrectly show the same
  values for R1 and R2.

v1.15 (2017-11-23)
------------------

* Cutadapt can now run on multiple CPU cores in parallel! To enable
  it, use the option ``-j N`` (or the long form ``--cores=N``), where ``N`` is
  the number of cores to use. Multi-core support is only available on Python 3,
  and not yet with some command-line arguments. See
  :ref:`the new section about multi-core in the documentation <multicore>`
  for details. When writing ``.gz`` files, make sure you have ``pigz`` installed
  to get the best speedup.
* The plan is to make multi-core the default (automatically using as many cores as
  are available) in future releases, so please test it and `report an
  issue <https://github.com/marcelm/cutadapt/issues/>`_ if you find problems!
* Issue :issue:`256`: ``--discard-untrimmed`` did not
  have an effect on non-anchored linked adapters.
* Issue :issue:`118`: Added support for demultiplexing of paired-end data.


v1.14 (2017-06-16)
------------------

* Fix: Statistics for 3' part of a linked adapter were reported incorrectly
* Fix `issue #244 <https://github.com/marcelm/cutadapt/issues/244>`_:
  Quality trimming with ``--nextseq-trim`` would not apply to R2 when
  trimming paired-end reads.
* ``--nextseq-trim`` now disables legacy mode.
* Fix `issue #246 <https://github.com/marcelm/cutadapt/issues/246>`_: installation
  failed on non-UTF8 locale

v1.13 (2017-03-16)
------------------

* The 3' adapter of linked adapters can now be anchored. Write
  ``-a ADAPTER1...ADAPTER2$`` to enable this. Note that the
  5' adapter is always anchored in this notation.
* Issue #224: If you want the 5' part of a linked adapter *not* to be
  anchored, you can now write ``-g ADAPTER...ADAPTER2`` (note ``-g``
  instead of ``-a``). This feature is experimental and may change behavior
  in the next release.
* Issue #236: For more accurate statistics, it is now possible to specify the
  GC content of the input reads with ``--gc-content``. This does
  not change trimming results, only the number in the "expect"
  column of the report. Since this is probably not needed by many
  people, the option is not listed when running ``cutadapt --help``.
* Issue #235: Adapter sequences are now required to contain only
  valid IUPAC codes (lowercase is also allowed, ``U`` is an alias
  for ``T``). This should help to catch hard-to-find bugs, especially
  in scripts. Use option ``-N`` to match characters literally
  (possibly useful for amino acid sequences).
* Documentation updates and some refactoring of the code

v1.12 (2016-11-28)
------------------

* Add read modification option ``--length`` (short: ``--l``), which will
  shorten each read to the given length.
* Cutadapt will no longer complain that it has nothing to do when you do not
  give it any adapters. For example, you can use this to convert file formats:
  ``cutadapt -o output.fasta input.fastq.gz`` converts FASTQ to FASTA.
* The ``xopen`` module for opening compressed files was moved to a `separate
  package on PyPI <https://pypi.python.org/pypi/xopen>`_.

v1.11 (2016-08-16)
------------------

* The ``--interleaved`` option no longer requires that both input and output
  is interleaved. It is now possible to have two-file input and interleaved
  output, and to have interleaved input and two-file output.
* Fix issue #202: First and second FASTQ header could get out of sync when
  options modifying the read name were used.

v1.10 (2016-05-19)
------------------

* Added a new “linked adapter” type, which can be used to search for a 5' and a
  3' adapter at the same time. Use ``-a ADAPTER1...ADAPTER2`` to search
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

v1.9.1 (2015-12-02)
-------------------

* Added ``--pair-filter`` option, which :ref:`modifies how filtering criteria
  apply to paired-end reads <filtering-paired>`
* Add ``--too-short-paired-output`` and ``--too-long-paired-output`` options.
* Fix incorrect number of trimmed bases reported if ``--times`` option was used.

v1.9 (2015-10-29)
-----------------

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

v1.8.3 (2015-07-29)
-------------------

* Fix issue #95: Untrimmed reads were not listed in the info file.
* Fix issue #138: pip install cutadapt did not work with new setuptools versions.
* Fix issue #137: Avoid a hang when writing to two or more gzip-compressed
  output files in Python 2.6.

v1.8.2 (2015-07-24)
-------------------

v1.8.1 (2015-04-09)
-------------------

* Fix #110: Counts for 'too short' and 'too long' reads were swapped in statistics.
* Fix #115: Make ``--trim-n`` work also on second read for paired-end data.

v1.8 (2015-03-14)
-----------------

* Support single-pass paired-end trimming with the new ``-A``/``-G``/``-B``/``-U``
  parameters. These work just like their -a/-g/-b/-u counterparts, but they
  specify sequences that are removed from the *second read* in a pair.

  Also, if you start using one of those options, the read modification options
  such as ``-q`` (quality trimming) are applied to *both* reads. For backwards
  compatibility, read modifications are applied to the first read only if
  neither of ``-A``/``-G``/``-B``/``-U`` is used. See `the
  documentation <http://cutadapt.readthedocs.io/en/latest/guide.html#paired-end>`_
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

v1.7 (2014-11-25)
-----------------

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

v1.6 (2014-10-07)
-----------------

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

v1.5 (2014-08-05)
-----------------

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

v1.4 (2014-03-13)
-----------------

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

v1.3 (2013-11-08)
-----------------

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

v1.2 (2012-11-30)
-----------------

* At least 25% faster processing of .csfasta/.qual files due to faster parser.
* Between 10% and 30% faster writing of gzip-compressed output files.
* Support 5' adapters in colorspace, even when no primer trimming is requested.
* Add the ``--info-file`` option, which has a line for each found adapter.
* Named adapters are possible. Usage: ``-a My_Adapter=ACCGTA`` assigns the name "My_adapter".
* Improve alignment algorithm for better poly-A trimming when there are sequencing errors.
  Previously, not the longest possible poly-A tail would be trimmed.
* James Casbon contributed the ``--discard-untrimmed`` option.

v1.1 (2012-06-18)
-----------------

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

v1.0 (2011-11-04)
-----------------

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

v0.9.5 (2011-07-20)
-------------------

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

v0.9.4 (2011-05-20)
-------------------

* now compatible with Python 3
* Add the ``--zero-cap`` option, which changes negative quality values to zero.
  This is a workaround to avoid segmentation faults in BWA. The option is now
  enabled by default when ``--bwa``/``--maq`` is used.
* Lots of unit tests added. Run them with ``cd tests && ./tests.sh``.
* Fix issue 16: ``--discard-trimmed`` did not work.
* Allow to override auto-detection of input file format with the new ``-f``/``--format``
  parameter. This mostly fixes issue 12.
* Don't break when input file is empty.

v0.9.2 (2011-03-16)
-------------------

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

v0.9 (2011-01-10)
-----------------

* add ``--too-short-output`` and ``--untrimmed-output``, based on patch by Paul Ryvkin (thanks!)
* add ``--maximum-length`` parameter: discard reads longer than a specified length
* group options by category in ``--help`` output
* add ``--length-tag`` option. allows to fix read length in FASTA/Q comment lines
  (e.g., ``length=123`` becomes ``length=58`` after trimming) (requested by Paul Ryvkin)
* add ``-q``/``--quality-cutoff`` option for trimming low-quality ends (uses the same algorithm
  as BWA)
* some refactoring
* the filename ``-`` is now interpreted as standard in or standard output

v0.8 (2010-12-08)
-----------------

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

v0.7 (2010-12-03)
-----------------

* Useful exit codes
* Better error reporting when malformed files are encountered
* Add ``--minimum-length`` parameter for discarding reads that are shorter than
  a specified length after trimming.
* Generalize the alignment function a bit. This is preparation for
  supporting adapters that are specific to either the 5' or 3' end.
* pure Python fallback for alignment function for when the C module cannot
  be used.

v0.6 (2010-11-18)
-----------------

* Support gzipped input and output.
* Print timing information in statistics.

v0.5 (2010-11-17)
-----------------

* add ``--discard`` option which makes cutadapt discard reads in which an adapter occurs

v0.4 (2010-11-17)
-----------------

* (more) correctly deal with multiple adapters: If a long adapter matches with lots of
  errors, then this could lead to a a shorter adapter matching with few errors getting ignored.

v0.3 (2010-09-27)
-----------------

* fix huge memory usage (entire input file was unintentionally read into memory)

v0.2 (2010-09-14)
-----------------

* allow FASTQ input

v0.1 (2010-09-14)
-----------------

* initial release
