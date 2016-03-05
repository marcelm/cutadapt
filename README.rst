.. image:: https://travis-ci.org/marcelm/cutadapt.svg?branch=master
    :target: https://travis-ci.org/marcelm/cutadapt

.. image:: https://img.shields.io/pypi/v/cutadapt.svg?branch=master
    :target: https://pypi.python.org/pypi/cutadapt

=================
cutadapt-parallel
=================

This is a working implementation of a multi-threaded version of Cutadapt. It is based on a fork of Cutadapt 1.9.2_dev. Most of the changes here will be integrated into the main Cutadapt program, however some options/functionality may change in the process.

Architecture: There is a main thread that reads from the input file(s) and posts batches of reads to a queue. There are one or more worker threads that take batches from the queue, process them, and post the results on a result queue. Finally, there is a worker thread that reads batches of results and writes them to a file.

Performance: As a general rule, we see linear speed increases compared to single threaded mode with more than one worker thread (e.g. --threads=2 gives a 2x performance boost versus single-threaded mode, --threads=8 gives 8-10x boost, etc).

Optimization: The key to maximizing cutadapt-parallel performance is to keep all threads working as much as possible. This can be controlled through four parameters:

* Threads: the number of *worker* threads to use. There will additionally be a main thread and writer thread. So if your system only has one or two cores, you may be better off running in single-threaded mode, which has less overhead.
* Batch size: the number of reads in each batch. If it takes longer for the reader to read a batch than it does for a worker to process that batch, the workers will end up blocking waiting for additional batches to process. We set the default batch size to 1000, however we recommend tuning this parameter to your own system. You can efficiently do this by limiting the number of reads run by the program ('--max-reads 1M', for example), turning on debug-level logging, and altering the batch size larger and larger until you start to see "Worker waiting for batch" messages. For example, we found that on our test system (actually a node in a computing cluster), with 8 threads, the optimal batch size was 5000.
* Read and result queue sizes: this is a trade-off between memory usage and probability of threads blocking (i.e., the larger the queue size, the greater the memory usage but lesser the chance of reader or worker threads blocking). We set this by default to 10 x number of threads, which seems to work well.

Note: A progress bar can be shown using the --progress option. This is off by default, because the progress bar involves some overhead (roughly 1 us/read penalty). In single-threaded mode, the ETA will be roughly accurate; however, in threaded mode, the progress bar will reach completion with some time remaining. This is because the progress bar tracks the number of reads loaded from the input file, but after reading is complete there will still be some tasks on the queue. Thus, the ETA remaining on the progress bar is roughly reflective of the time the program will take to complete after progress reaches 100%.

========
cutadapt
========

Cutadapt finds and removes adapter sequences, primers, poly-A tails and other
types of unwanted sequence from your high-throughput sequencing reads.

Cleaning your data in this way is often required: Reads from small-RNA
sequencing contain the 3’ sequencing adapter because the read is longer than
the molecule that is sequenced. Amplicon reads start with a primer sequence.
Poly-A tails are useful for pulling out RNA from your sample, but often you
don’t want them to be in your reads.

Cutadapt helps with these trimming tasks by finding the adapter or primer
sequences in an error-tolerant way. It can also modify and filter reads in
various ways. Adapter sequences can contain IUPAC wildcard characters. Also,
paired-end reads and even colorspace data is supported. If you want, you can
also just demultiplex your input data, without removing adapter sequences at all.

Cutadapt comes with an extensive suite of automated tests and is available under
the terms of the MIT license.

If you use cutadapt, please cite
`DOI:10.14806/ej.17.1.200 <http://dx.doi.org/10.14806/ej.17.1.200>`_ .


Links
-----

* `Documentation <https://cutadapt.readthedocs.org/>`_
* `Source code <https://github.com/marcelm/cutadapt/>`_
* `Report an issue <https://github.com/marcelm/cutadapt/issues>`_
* `Project page on PyPI (Python package index) <https://pypi.python.org/pypi/cutadapt/>`_
* `Follow @marcelm_ on Twitter <https://twitter.com/marcelm_>`_
* `Wrapper for the Galaxy platform <https://bitbucket.org/lance_parsons/cutadapt_galaxy_wrapper>`_
