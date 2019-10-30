============
Installation
============

Cutadapt is being developed and tested under Linux. Users have run it
successfully under macOS and Windows.


Quick installation
------------------

The easiest way to install Cutadapt is to use ``pip3`` on the command line::

    pip3 install --user --upgrade cutadapt

This will download the software from `PyPI (the Python packaging
index) <https://pypi.python.org/pypi/cutadapt/>`_, and
install the cutadapt binary into ``$HOME/.local/bin``. If an old version of
Cutadapt exists on your system, the ``--upgrade`` parameter is required in order
to install a newer version. You can then run the program like this::

    ~/.local/bin/cutadapt --help

If you want to avoid typing the full path, add the directory
``$HOME/.local/bin`` to your ``$PATH`` environment variable.


Installation with conda
-----------------------

Alternatively, Cutadapt is available as a conda package from the
`bioconda channel <https://bioconda.github.io/>`_. If you do not have conda,
`install miniconda <http://conda.pydata.org/miniconda.html>`_ first.
Then install Cutadapt like this::

    conda install -c bioconda cutadapt

If neither ``pip`` nor ``conda`` installation works, keep reading.


Installation on a Debian-based Linux distribution
-------------------------------------------------

Cutadapt is also included in Debian-based Linux distributions, such as Ubuntu.
Simply use your favorite package manager to install Cutadapt. On the
command-line, this should work ::

    sudo apt install cutadapt

or possibly ::

    sudo apt install python3-cutadapt

Please be aware that this will likely give you an old version of Cutadapt. If
you encounter unexpected behavior, please use one of the other installation
methods to get an up-to-date version before reporting bugs.


.. _dependencies:

Dependencies
------------

Cutadapt installation requires this software to be installed:

* Python 3.4 or newer
* Possibly a C compiler. For Linux, Cutadapt packages are provided as
  so-called “wheels” (``.whl`` files) which come pre-compiled.

Under Ubuntu, you may need to install the packages ``build-essential`` and
``python3-dev`` to get a C compiler.

If you get an error message::

    error: command 'gcc' failed with exit status 1

Then check the entire error message. If it says something about a missing
``Python.h`` file, then the problem is that you are missing Python development
packages (``python3-dev`` in Ubuntu).


System-wide installation (root required)
----------------------------------------

If you have root access, then you can install Cutadapt system-wide by running::

    sudo python3 -m pip install cutadapt

This installs cutadapt into ``/usr/local/bin``.

If you want to upgrade from an older version, use this command instead::

    sudo python3 -m pip install --upgrade cutadapt


If the above does not work for you, then you can try to install Cutadapt
into a virtual environment. This may lead to fewer conflicts with
system-installed packages::

    sudo python3 -m venv /usr/local/cutadapt
    sudo /usr/local/cutadapt/bin/pip install cutadapt
    cd /usr/local/bin/
    sudo ln -s ../cutadapt/bin/cutadapt


Uninstalling
------------

Type  ::

    pip3 uninstall cutadapt

and confirm with ``y`` to remove the package. Under some circumstances, multiple
versions may be installed at the same time. Repeat the above command until you
get an error message in order to make sure that all versions are removed.


Shared installation (on a cluster)
----------------------------------

If you have a larger installation and want to provide Cutadapt as a module
that can be loaded and unloaded (with the Lmod system, for example), we
recommend that you create a virtual environment and 'pip install' cutadapt into
it. These instructions work on our SLURM cluster that uses the Lmod system
(replace ``1.9.1`` with the actual version you want to use)::

    BASE=/software/cutadapt-1.9.1
    virtualenv $BASE/venv
    $BASE/venv/bin/pip install --install-option="--install-scripts=$BASE/bin" cutadapt==1.9.1

The ``install-option`` part is important. It ensures that a second, separate
``bin/`` directory is created (``/software/cutadapt-1.9.1/bin/``) that *only*
contains the ``cutadapt`` script and nothing else. To make Cutadapt available to
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

We recommend that you install Cutadapt into a so-called virtual environment if
you decide to use the development version. The virtual environment is a single
directory that contains everything needed to run the software. Nothing else on
your system is changed, so you can simply uninstall this particular version of
Cutadapt by removing the directory with the virtual environment.

The following instructions work on Linux using Python 3. Make sure you have
installed the :ref:`dependencies <dependencies>` (``python3-dev`` and
``build-essential`` on Ubuntu)!

First, choose where you want to place the directory with the virtual
environment and what you want to call it. Let us assume you chose the path
``~/cutadapt-venv``. Then use these commands for the installation::

    python3 -m venv ~/cutadapt-venv
    ~/cutadapt-venv/bin/python3 -m pip install --upgrade pip
    ~/cutadapt-venv/bin/pip install git+https://github.com/marcelm/cutadapt.git#egg=cutadapt

To run Cutadapt and see the version number, type ::

    ~/cutadapt-venv/bin/cutadapt --version

The reported version number will be something like ``2.2.dev5+gf564208``. This
means that you are now running the version of Cutadapt that will become 2.2, and that it contains
5 changes (*commits*) since the previous release (2.1 in this case).
