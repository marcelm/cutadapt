.. image:: https://travis-ci.org/marcelm/cutadapt.svg?branch=master
    :target: https://travis-ci.org/marcelm/cutadapt

========
cutadapt
========

Cutadapt finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequence from your high-throughput sequencing reads.

Cleaning your data in this way is often required: Reads from small-RNA sequencing contain the 3’ sequencing adapter because the read is longer than the molecule that is sequenced. Amplicon reads start with a primer sequence. Poly-A tails are useful for pulling out RNA from your sample, but often you don’t want them to be in your reads.

Cutadapt helps with these trimming tasks by finding the adapter or primer sequences in an error-tolerant way. It can also filter reads by length and do quality trimming. Adapter sequences can contain IUPAC wildcard characters. Also, paired-end reads and even colorspace data is supported. If you want, you can also just demultiplex your input data, without removing adapter sequences at all.

Cutadapt comes with an extensive suite of automated tests and is available under the terms of the MIT license.


Links
-----

* `Project homepage <http://code.google.com/p/cutadapt/>`_
* `Github page <https://github.com/marcelm/cutadapt/>`_
* `Read the documentation online <https://cutadapt.readthedocs.org/>`_.
  Offline documentation is available in the ``doc/`` subdirectory in the repository
  and in the downloaded tar distribution.
* `Report issues to the Google code bug tracker <https://code.google.com/p/cutadapt/issues/list>`_
* `A cutadapt wrapper for the Galaxy platform <https://bitbucket.org/lance_parsons/cutadapt_galaxy_wrapper>`_,
  written by Lance Parsons.
