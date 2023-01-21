============
Installation
============

Because Cutadapt development happens on Linux, this is the best supported
platform, but it should also run on macOS and Windows.


Installation with Conda
-----------------------

Cutadapt is available as a Conda package from the
`Bioconda channel <https://bioconda.github.io/>`_.
`Install miniconda <http://conda.pydata.org/miniconda.html>`_ if
you don’t have Conda. Then follow the `Bioconda installation
instructions <https://bioconda.github.io/user/install.html>`_ (in particular,
make sure you have both `bioconda` and `conda-forge` in your channels list).

To then install Cutadapt into a new Conda environment, use this command::

    conda create -n cutadaptenv cutadapt

Here, ``cutadaptenv`` is the name of the Conda environment. You can
choose a different name.

If you are on macOS and your machine uses an M1/M2 processor (Apple Silicon),
you may need to run this command instead::

       CONDA_SUBDIR=osx-64 conda create -n cutadaptenv cutadapt

(If you have problems, see `this issue for troubleshooting
<https://github.com/marcelm/cutadapt/issues/672>`_.)

Then activate the environment. This needs to be done every time you open a
new shell before you can use Cutadapt::

    conda activate cutadaptenv

Finally, check whether it worked::

    cutadapt --version

This should show the Cutadapt version number.


Installation with pipx
----------------------

This works on Ubuntu 20.04 and later:

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


.. _dependencies:

Dependencies
------------

Cutadapt installation requires this software to be installed:

* Python 3.7 or newer
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

Generally, using ``sudo`` can be dangerous and the above methods that don’t
require it are preferred. That said, if you have root access, you can install
Cutadapt system-wide by running::

    sudo python3 -m pip install cutadapt

This installs cutadapt into ``/usr/local/bin``.

If you want to upgrade from an older version, use this command instead::

    sudo python3 -m pip install --upgrade cutadapt

If the above does not work for you, then you can try to install Cutadapt
into a virtual environment. This leads to fewer conflicts with
system-installed packages::

    sudo python3 -m venv /usr/local/cutadapt
    sudo /usr/local/cutadapt/bin/pip install cutadapt
    cd /usr/local/bin/
    sudo ln -s ../cutadapt/bin/cutadapt


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


Uninstalling
------------

Type ::

    pip3 uninstall cutadapt

and confirm with ``y`` to remove the package. Under some circumstances, multiple
versions may be installed at the same time. Repeat the above command until you
get an error message in order to make sure that all versions are removed.


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
    ~/cutadapt-venv/bin/pip install git+https://github.com/marcelm/cutadapt.git

To run Cutadapt and see the version number, type ::

    ~/cutadapt-venv/bin/cutadapt --version

The reported version number will be something like ``2.2.dev5+gf564208``. This
means that you are now running the version of Cutadapt that will become 2.2, and that it contains
5 changes (*commits*) since the previous release (2.1 in this case).
