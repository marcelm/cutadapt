.. image:: https://github.com/marcelm/cutadapt/workflows/CI/badge.svg
    :alt:

.. image:: https://img.shields.io/pypi/v/cutadapt.svg?branch=master
    :target: https://pypi.python.org/pypi/cutadapt
    :alt:

.. image:: https://codecov.io/gh/marcelm/cutadapt/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/marcelm/cutadapt
    :alt:

.. image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat
    :target: http://bioconda.github.io/recipes/cutadapt/README.html
    :alt: install with bioconda


========
Cutadapt
========

Cutadapt finds and removes adapter sequences, primers, poly-A tails and other
types of unwanted sequence from your high-throughput sequencing reads.

Cleaning your data in this way is often required: Reads from small-RNA
sequencing contain the 3’ sequencing adapter because the read is longer than
the molecule that is sequenced. Amplicon reads start with a primer sequence.
Poly-A tails are useful for pulling out RNA from your sample, but often you
don’t want them to be in your reads.

Cutadapt helps with these trimming tasks by finding the adapter or primer
sequences in an error-tolerant way. It can also modify and filter single-end
and paired-end reads in various ways. Adapter sequences can contain IUPAC
wildcard characters. Cutadapt can also demultiplex your reads.

Cutadapt is available under the terms of the MIT license.

Cutadapt development was started at `TU Dortmund University <https://www.tu-dortmund.de>`_
in the group of `Prof. Dr. Sven Rahmann <https://www.rahmannlab.de/>`_.
It is currently being developed within
`NBIS (National Bioinformatics Infrastructure Sweden) <https://nbis.se/>`_.

If you use Cutadapt, please cite
`DOI:10.14806/ej.17.1.200 <http://dx.doi.org/10.14806/ej.17.1.200>`_ .


Links
-----

* `Documentation <https://cutadapt.readthedocs.io/>`_
* `Source code <https://github.com/marcelm/cutadapt/>`_
* `Report an issue <https://github.com/marcelm/cutadapt/issues>`_
* `Project page on PyPI (Python package index) <https://pypi.python.org/pypi/cutadapt/>`_
* `Wrapper for the Galaxy platform <https://github.com/galaxyproject/tools-iuc/tree/master/tools/cutadapt>`_
