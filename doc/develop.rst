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


Contributing
------------

Contributions to Cutadapt in the form of source code or documentation
improvements or helping out with responding to issues are welcome!

To contribute to Cutadapt development, send in a pull request (PR) on GitHub.

* Limit a PR to a single topic. Submit multiple PRs if necessary. This way, it
  is easier to discuss the changes individually, and in case we find that one
  of them should not go in, the others can still be accepted.
* For larger changes, consider opening an issue first to plan what you want to
  do.
* Include appropriate unit or integration tests. Sometimes, tests are hard to
  write or don’t make sense. If you think this is the case, just leave the tests
  out initially and we can discuss whether to add any.
* Add documentation and a changelog entry if appropriate.


Code style
~~~~~~~~~~

* Cutadapt tries to follow PEP8, except that the allowed line length is 100
  characters, not 80. But try to wrap comments after 80 characters.
* There are inconsistencies in the current code base since it’s a few years old
  already. New code should follow the current rules, however.
* At the moment, no automatic code formatting is done, but one idea might be to
  switch to the `black <https://black.readthedocs.io/>` code formatter at some
  point. If you’re familiar with its style, you can use that already now for
  new code to make the diff smaller.
* Prefer double quotation marks in new code. This will also make the diff smaller
  when and if we eventually switch to black.
* Using an IDE is beneficial (PyCharm, for example). It helps to catch lots of
  style issues early (unused imports, spacing etc.).
* Avoid unnecessary abbreviations for variable names. Code is more often read
  than written.
* When writing a help text for a new command-line option, look at the output of
  ``cutadapt --help`` and try to make it look nice and short.
* In comments and documentation, capitalize FASTQ, BWA, CPU etc.


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
