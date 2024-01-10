============
Installation
============

Because Cutadapt development happens on Linux, this is the best supported
platform, but it should also run on macOS and Windows.


Installation with Conda
-----------------------

Cutadapt is available as a Conda package from the
`Bioconda channel <https://bioconda.github.io/>`_.

1. Install Conda. For example, by
   `installing miniforge <https://github.com/conda-forge/miniforge#install>`_.

2. Configure the Bioconda channel by following the
   `Bioconda setup instructions <https://bioconda.github.io/#usage>`_.
   In short::

     conda config --add channels bioconda
     conda config --add channels conda-forge
     conda config --set channel_priority strict

   (The Bioconda instructions mention the ``defaults`` channel,
   but it is not needed.)

3. Install Cutadapt into a new Conda environment::

     conda create -n cutadapt cutadapt

   The first ``cutadapt`` in this command is the name of the Conda environment.
   You can choose a different name.

   If you are on macOS and your machine uses an M1/M2 processor (Apple Silicon),
   you may need to run this command instead::

     CONDA_SUBDIR=osx-64 conda create -n cutadapt cutadapt

   (If you have problems, see `this issue for troubleshooting
   <https://github.com/marcelm/cutadapt/issues/672>`_.)

4. Activate the Conda environment. This needs to be done every time you open a
   new shell in order to be able to use Cutadapt::

     conda activate cutadaptenv

5. Finally, check whether the installation was successful::

     cutadapt --version

   This should show the Cutadapt version number.


Installation with pipx
----------------------

This works on Ubuntu 20.04 and later::

    sudo apt install pipx python3-venv
    pipx install cutadapt
    cutadapt --version


Installation with pip
---------------------

Ensure you have virtualenv installed. On Ubuntu/Debian::

    sudo apt install python3-virtualenv

Create a new virtual environment and install Cutadapt into it::

    virtualenv cutadapt-venv
    cutadapt-venv/bin/pip --upgrade pip
    cutadapt-venv/bin/pip install cutadapt

Cutadapt is now available as `cutadapt-venv/bin/cutadapt`::

    cutadapt-venv/bin/cutadapt --version

Optionally, you can *activate* the virtual environment, which allows you to
just type `cutadapt` without the full path::

    source cutadapt-venv/bin/activate
    cutadapt --version

Activation must be re-done whenever you open a new shell/terminal window.


Installation on Debian/Ubuntu
-----------------------------

Cutadapt is also included in Debian-based Linux distributions, such as Ubuntu.
Simply use your favorite package manager to install Cutadapt. On the
command-line, this should work ::

    sudo apt install cutadapt

or possibly ::

    sudo apt install python3-cutadapt

Please be aware that distribution packages are very likely to be outdated.
If you encounter unexpected behavior or need newer features, please use one
of the other installation methods to get an up-to-date version before
reporting bugs.


Installation on Windows
-----------------------

For some releases of Cutadapt, a single-file executable (``cutadapt.exe``)
is made available on the
`GitHub releases page <https://github.com/marcelm/cutadapt/releases>`_. Try that
first, and if it does not work for you, please report the issue.

To install Cutadapt manually, keep reading.

There is no Bioconda package for Windows because Bioconda does not produce
Windows packages. To install Cutadapt, you can use ``pip``, but because
Cutadapt contains components that need to be compiled, you also need to install
a compiler.

1. Download a recent version (at least 3.7) of Python for Windows from
   <https://www.python.org/> and install it.
2. Download and install “Build Tools for Visual Studio 2019” from
   <https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2019>.
   (There are many similarly named downloads on that page, ensure you get the
   right one.)

   During installation, when the dialog about which components to install pops
   up, ensure that “C++ Build tools” is ticked. The download is quite big and
   can take a long time.
3. Open the command line (``cmd.exe``) and run ``py -m pip install cutadapt``.
4. Test whether it worked by running ``py -m cutadapt --version``. You should
   see the version number of Cutadapt.

When running Cutadapt this way, you will need to remember to write
``py -m cutadapt`` instead of just ``cutadapt``.


Shared installation (on a cluster)
----------------------------------

If you have a larger installation and want to provide Cutadapt as a module
that can be loaded and unloaded (with the Lmod system, for example), we
recommend that you create a virtual environment and 'pip install' Cutadapt into
it. These instructions work on a SLURM cluster that uses the Lmod system
(replace ``3.1`` with the actual version you want to use)::

    BASE=/software/cutadapt-3.1
    virtualenv $BASE/venv
    $BASE/venv/bin/pip install cutadapt==3.1
    mkdir $BASE/bin
    cd $BASE/bin
    ln -s ../venv/bin/cutadapt

Then add the directory ``$BASE/bin/`` to the ``$PATH`` when a user loads the
module, somewhat like this (this is for the Lmod system)::

    conflict("cutadapt")
    whatis("adapter trimming tool")
    prepend_path("PATH", "/software/cutadapt-3.1/bin")

Make sure that you **do not** add ``$BASE/venv/bin/`` to the ``$PATH``!
Otherwise, a user trying to run ``python`` who also has the
cutadapt module loaded would get the python from the virtual environment,
which leads to confusing error messages. The ``$BASE/bin/`` directory only
contains the ``cutadapt`` script and nothing else, avoiding this problem.

Please note that there is no need to “activate” virtual environments.


Installing the development version
----------------------------------

We recommend that you install Cutadapt into a so-called virtual environment if
you decide to use the development version. The virtual environment is a single
directory that contains everything needed to run the software. Nothing else on
your system is changed, so you can uninstall this particular version of
Cutadapt by just removing the directory with the virtual environment.

The following instructions work on Linux using Python 3. Make sure you have
installed the ``python3-dev`` and ``build-essential`` packages on Ubuntu.

First, choose where you want to place the directory with the virtual
environment and what you want to call it. Let us assume you chose the path
``~/cutadapt-venv``. Then use these commands for the installation::

    python3 -m venv ~/cutadapt-venv
    ~/cutadapt-venv/bin/python3 -m pip install --upgrade pip
    ~/cutadapt-venv/bin/pip install git+https://github.com/marcelm/cutadapt.git

To run Cutadapt and see the version number, type ::

    ~/cutadapt-venv/bin/cutadapt --version

The reported version number will be something like ``2.2.dev5+gf564208``. This
means that you are now running the version of Cutadapt that will become 2.2,
and that it contains 5 changes (*commits*) since the previous release (2.1 in this case).
