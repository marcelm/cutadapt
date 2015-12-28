Colorspace reads
================

Cutadapt was designed to work with colorspace reads from the ABi SOLiD
sequencer. Colorspace trimming is activated by the ``--colorspace``
option (or use ``-c`` for short). The input reads can be given either:

-  in a FASTA file
-  in a FASTQ file
-  in a ``.csfasta`` and a ``.qual`` file (this is the native SOLiD
   format).

In all cases, the colors must be represented by the characters 0, 1, 2,
3. Example input files are in the cutadapt distribution at
``tests/data/solid.*``. The ``.csfasta``/``.qual`` file format is
automatically assumed if two input files are given to cutadapt.

In colorspace mode, the adapter sequences given to the ``-a``, ``-b``
and ``-g`` options can be given both as colors or as nucleotides. If
given as nucleotides, they will automatically be converted to
colorspace. For example, to trim an adapter from ``solid.csfasta`` and
``solid.qual``, use this command-line::

    cutadapt -c -a CGCCTTGGCCGTACAGCAG solid.csfasta solid.qual > output.fastq

In case you know the colorspace adapter sequence, you can also write
``330201030313112312`` instead of ``CGCCTTGGCCGTACAGCAG`` and the result
is the same.

Ambiguity in colorspace
-----------------------

The ambiguity of colorspace encoding leads to some effects to be aware
of when trimming 3' adapters from colorspace reads. For example, when
trimming the adapter ``AACTC``, cutadapt searches for its
colorspace-encoded version ``0122``. But also ``TTGAG``, ``CCAGA`` and
``GGTCT`` have an encoding of ``0122``. This means that effectively four
different adapter sequences are searched and trimmed at the same time.
There is no way around this, unless the decoded sequence were available,
but that is usually only the case after read mapping.

The effect should usually be quite small. The number of false positives
is multiplied by four, but with a sufficiently large overlap (3 or 4 is
already enough), this is still only around 0.2 bases lost per read on
average. If inspecting k-mer frequencies or using small overlaps, you
need to be aware of the effect, however.


Double-encoding, BWA and MAQ
----------------------------

The read mappers MAQ and BWA (and possibly others) need their colorspace
input reads to be in a so-called "double encoding". This simply means
that they cannot deal with the characters 0, 1, 2, 3 in the reads, but
require that the letters A, C, G, T be used for colors. For example, the
colorspace sequence ``0011321`` would be ``AACCTGC`` in double-encoded
form. This is not the same as conversion to basespace! The read is still
in colorspace, only letters are used instead of digits. If that sounds
confusing, that is because it is.

Note that MAQ is unmaintained and should not be used in new projects.

BWA’s colorspace support was dropped in versions more recent than 0.5.9,
but that version works well.

When you want to trim reads that will be mapped with BWA or MAQ, you can
use the ``--bwa`` option, which enables colorspace mode (``-c``),
double-encoding (``-d``), primer trimming (``-t``), all of which are
required for BWA, in addition to some other useful options.

The ``--maq`` option is an alias for ``--bwa``.


Colorspace examples
-------------------

To cut an adapter from SOLiD data given in ``solid.csfasta`` and
``solid.qual``, to produce MAQ- and BWA-compatible output, allow the
default of 10% errors and write the resulting FASTQ file to
output.fastq::

    cutadapt --bwa -a CGCCTTGGCCGTACAGCAG solid.csfasta solid.qual > output.fastq

Instead of redirecting standard output with ``>``, the ``-o`` option can
be used. This also shows that you can give the adapter in colorspace and
how to use a different error rate::

    cutadapt --bwa -e 0.15 -a 330201030313112312 -o output.fastq solid.csfasta solid.qual

This does the same as above, but produces BFAST-compatible output,
strips the \_F3 suffix from read names and adds the prefix "abc:" to
them::

    cutadapt -c -e 0.15 -a 330201030313112312 -x abc: --strip-f3 solid.csfasta solid.qual > output.fastq


Bowtie
------

Quality values of colorspace reads are sometimes negative. Bowtie gets
confused and prints this message::

    Encountered a space parsing the quality string for read xyz

BWA also has a problem with such data. Cutadapt therefore converts
negative quality values to zero in colorspace data. Use the option
``--no-zero-cap`` to turn this off.

.. _sra-fastq:

Sequence Read Archive
---------------------

The Sequence Read Archive provides files in a special "SRA" file format. When
the ``fastq-dump`` program from the sra-toolkit package is used to convert
these ``.sra`` files to FASTQ format, colorspace reads will get an extra
quality value in the beginning of each read. You may get an error like this::

    cutadapt: error: In read named 'xyz': length of colorspace quality
    sequence (36) and length of read (35) do not match (primer is: 'T')

To make cutadapt ignore the extra quality base, add ``--format=sra-fastq`` to
your command-line, as in this example::

    cutadapt -c --format=sra-fastq -a CGCCTTGGCCG sra.fastq > trimmed.fastq

When you use ``--format=sra-fastq``, the spurious quality value will be removed
from all reads in the file.
