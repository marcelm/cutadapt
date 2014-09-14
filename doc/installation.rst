Installation
============

Dependencies
------------

cutadapt needs Python 2.6 or later (this includes Python 3). Python 2.6
supported is not tested thoroughly and is also slower than later
versions. For installation from sources, a C compiler needs to be
installed. The program has been developed and tested on Ubuntu and
OpenSuSE.

Installation
------------

The easiest way to install cutadapt is via the ``pip`` command if that
is available on your system::

    pip install --user cutadapt

This installs the cutadapt binary into ``$HOME/.local/bin``. Make sure
that this directory is on your ``$PATH``.

If you have already downloaded and unpacked the ``.tar.gz`` file, then
installation is done like this (replace "python" with "python3" to
install the Python 3 version)::


    python setup.py install --user

If you get an error about a missing ``Python.h`` file, then make sure
that the ``python-dev`` package is installed (or ``python3-dev`` for
Python 3).

Use without installation
------------------------

Build the C extension module (you can try to skip this step -- a
compiled version of the module for Linux x86\_64 is already included)::

    python setup.py build_ext -i

Then simply run the script from where it is, similar to this::

    bin/cutadapt --help

If you get any errors, first try to explicitly request a specific Python
version by running cutadapt like this::

    python2.7 bin/cutadapt --help
