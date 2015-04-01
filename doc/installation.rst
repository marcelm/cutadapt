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

If the above does not work, keep reading.


Dependencies
------------

Cutadapt requires this software to be installed:

* One of Python 2.6, 2.7, 3.3 or 3.4. Python 2.7 is a bit faster than the other
  versions.
* A C compiler.

Under Ubuntu, you may need to install the packages ``build-essential`` and
``python-dev``.


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


System-wide installation
------------------------

If you have root access, then you can install cutadapt system-wide by running::

    sudo pip install cutadapt

This installs cutadapt into `/usr/local/bin`.

If you want to upgrade from an older version, use this command instead::

    sudo pip install --upgrade cutadapt


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
