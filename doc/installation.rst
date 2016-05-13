============
Installation
============

Quickstart
----------

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

Alternatively, cutadapt is also available as a conda package from the
`bioconda channel <https://bioconda.github.io/>`_. If you do not have conda,
`install miniconda <http://conda.pydata.org/miniconda.html>`_ first.
Then install cutadapt like this::

    conda install -c bioconda cutadapt

If neither `pip` nor `conda` installation works, keep reading.


Dependencies
------------

Cutadapt requires this software to be installed:

* One of Python 2.6, 2.7, 3.3, 3.4 or 3.5. Python 2.7 is a bit faster than the
  other versions.
* A C compiler.

Under Ubuntu, you may need to install the packages ``build-essential`` and
``python-dev`` (or ``python3-dev``).


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
