Contributing
------------

Contributions to Cutadapt in the form of source code or documentation
improvements or helping out with responding to issues are welcome!

To contribute to Cutadapt development, it is easiest to send in a pull request
(PR) on GitHub.

Here are some guidelines for how to do this. They are not strict rules. When in
doubt, send in a PR and we will sort it out.

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

* The source code needs to be formatted with
  `black <https://black.readthedocs.io/>`_.
  If you install `pre-commit <https://pre-commit.com>`_,
  the formatting will be done for you.
* There are inconsistencies in the current code base since it’s a few years old
  already. New code should follow the current rules, however.
* Using an IDE is beneficial (PyCharm, for example). It helps to catch lots of
  style issues early (unused imports, spacing etc.).
* Use `Google-style docstrings <https://www.sphinx-doc.org/en/master/usage/extensions/example_google.html>`_
  (this is not PyCharm’s default setting).
* Avoid unnecessary abbreviations for variable names. Code is more often read
  than written.
* When writing a help text for a new command-line option, look at the output of
  ``cutadapt --help`` and try to make it look nice and short.
* In comments and documentation, capitalize FASTQ, BWA, CPU etc.
