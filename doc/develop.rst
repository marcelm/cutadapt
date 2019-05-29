Developing
==========

The `Cutadapt source code is on GitHub <https://github.com/marcelm/cutadapt/>`_.
Cutadapt is written in Python 3 with some extension modules that are written
in Cython. Support for Python 2 has been dropped.


Development installation
------------------------

For development, make sure that you install Cython and tox. We also recommend
using a virtualenv. This sequence of commands should work::

    git clone https://github.com/marcelm/cutadapt.git  # or clone your own fork
    cd cutadapt
    python3 -m venv venv
    venv/bin/pip3 install Cython pytest nose tox
    venv/bin/pip3 install -e .

Then you can run Cutadapt like this (or activate the virtualenv and omit the
``venv/bin`` part)::

    venv/bin/cutadapt --help

The tests can then be run like this::

    venv/bin/pytest

Or with tox (but then you will need to have binaries for all tested Python
versions installed)::

    venv/bin/tox


Making a release
----------------

Since version 1.17, Travis CI is used to automatically deploy a new Cutadapt release
(both as an sdist and as wheels) whenever a new tag is pushed to the Git repository.

Cutadapt uses `setuptools_scm <https://github.com/pypa/setuptools_scm>`_ to automatically manage
version numbers. This means that the version is not stored in the source code but derived from
the most recent Git tag. The following procedure can be used to bump the version and make a new
release.

#. Update ``CHANGES.rst`` (version number and list of changes)

#. Ensure you have no uncommitted changes in the working copy.

#. Run a ``git pull``.

#. Run ``tox``, ensuring all tests pass.

#. Tag the current commit with the version number (there must be a ``v`` prefix)::

       git tag v0.1

   To release a development version, use a ``dev`` version number such as ``v1.17.dev1``.
   Users will not automatically get these unless they use ``pip install --pre``.

#. Push the tag::

       git push --tags

#. Wait for Travis to finish and to deploy to PyPI.

#. The `bioconda recipe <https://github.com/bioconda/bioconda-recipes/blob/master/recipes/cutadapt/meta.yaml>`_
   also needs to be updated, but the bioconda bot will likely do this automatically
   if you just wait a little while.

   Ensure that the list of dependencies (the ``requirements:``
   section in the recipe) is in sync with the ``setup.py`` file.

If something went wrong *after* a version has already been tagged and published to
PyPI, fix the problem and tag a new version. Do not change a version that has already
been uploaded.


.. include:: ../CONTRIBUTING.rst
