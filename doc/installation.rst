============
Installation
============

Cutadapt is being developed and tested under Linux. Users have run it
successfully under macOS and Windows.


Quick installation
------------------

The easiest way to install cutadapt is to use ``pip`` on the command line::

    pip install --user --upgrade cutadapt

This will download the software from `PyPI (the Python packaging
index) <https://pypi.python.org/pypi/cutadapt/>`_, and
install the cutadapt binary into ``$HOME/.local/bin``. If an old version of
cutadapt exists on your system, the ``--upgrade`` parameter is required in order
to install a newer version. You can then run the program like this::

    ~/.local/bin/cutadapt --help

If you want to avoid typing the full path, add the directory
``$HOME/.local/bin`` to your ``$PATH`` environment variable.


Installation with conda
-----------------------

Alternatively, cutadapt is available as a conda package from the
`bioconda channel <https://bioconda.github.io/>`_. If you do not have conda,
`install miniconda <http://conda.pydata.org/miniconda.html>`_ first.
Then install cutadapt like this::

    conda install -c bioconda cutadapt

If neither `pip` nor `conda` installation works, keep reading.


.. _dependencies:

Dependencies
------------

Cutadapt installation requires this software to be installed:

* Python 2.7 or at least Python 3.3
* A C compiler.

Under Ubuntu, you may need to install the packages ``build-essential`` and
``python-dev`` (or ``python3-dev``) to get a C compiler.

On Windows, you need `Microsoft Visual C++ Compiler for
Python 2.7 <https://www.microsoft.com/en-us/download/details.aspx?id=44266>`_.


Installation
------------

If you have already downloaded and unpacked the ``.tar.gz`` file, then
installation is done like this (replace "python" with "python3" to
install the Python 3 version)::

    python setup.py install --user

If you get an error message::

    error: command 'gcc' failed with exit status 1

Then check the entire error message. If it says something about a missing ``Python.h``
file, then you need to install the Python development packages. The
appropriate package is called ``python-dev`` in Ubuntu (or ``python3-dev``
for Python 3).


System-wide installation (root required)
----------------------------------------

If you have root access, then you can install cutadapt system-wide by running::

    sudo pip install cutadapt

This installs cutadapt into `/usr/local/bin`.

If you want to upgrade from an older version, use this command instead::

    sudo pip install --upgrade cutadapt


Uninstalling
------------

Type  ::

    pip uninstall cutadapt

and confirm with ``y`` to remove the package. Under some circumstances, multiple
versions may be installed at the same time. Repeat the above command until you
get an error message in order to make sure that all versions are removed.


Shared installation (on a cluster)
----------------------------------

If you have a larger installation and want to provide cutadapt as a module
that can be loaded and unloaded (with the Lmod system, for example), we
recommend that you create a virtual environment and 'pip install' cutadapt into
it. These instructions work on our SLURM cluster that uses the Lmod system
(replace ``1.9.1`` with the actual version you want to use)::

    BASE=/software/cutadapt-1.9.1
    virtualenv $BASE/venv
    $BASE/venv/bin/pip install --install-option="--install-scripts=$BASE/bin" cutadapt==1.9.1

The ``install-option`` part is important. It ensures that a second, separate
``bin/`` directory is created (``/software/cutadapt-1.9.1/bin/``) that *only*
contains the ``cutadapt`` script and nothing else. To make cutadapt available to
the users, that directory (``$BASE/bin``) needs to be added to the ``$PATH``.

Make sure you *do not* add the ``bin/`` directory within the ``venv`` directory
to the ``$PATH``! Otherwise, a user trying to run ``python`` who also has the
cutadapt module loaded would get the python from the virtual environment,
which leads to confusing error messages.

A simple module file for the Lmod system matching the above example could look
like this::

    conflict("cutadapt")
    whatis("adapter trimming tool")
    prepend_path("PATH", "/software/cutadapt-1.9.1/bin")

Please note that there is no need to “activate” the virtual environment:
Activation merely adds the ``bin/`` directory to the ``$PATH``, so the
``prepend_path`` directive is equivalent to activating the virtual environment.


Installing the development version
----------------------------------

We recommend that you install cutadapt into a so-called virtual environment if
you decide to use the development version. The virtual environment is a single
directory that contains everything needed to run the software. Nothing else on
your system is changed, so you can simply uninstall this particular version of
cutadapt by removing the directory with the virtual environment.

The following instructions work on Linux using Python 3. Make sure you have
installed the :ref:`dependencies <dependencies>` (``python3-dev`` and
``build-essential`` on Ubuntu)!

First, choose where you want to place the directory with the virtual
environment and what you want to call it. Let us assume you chose the path
``~/cutadapt-venv``. Then use these commands for the installation::

    python3 -m venv ~/cutadapt-venv
    ~/cutadapt-venv/bin/pip install Cython
    ~/cutadapt-venv/bin/pip install git+https://github.com/marcelm/cutadapt.git/

To run cutadapt and see the version number, type ::

    ~/cutadapt-venv/bin/cutadapt --version

The reported version number will be something like ``1.14+65.g5610275``. This
means that you are now running a cutadapt version that contains 65 additional
changes (*commits*) since version 1.14.
