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

* Cutadapt tries to follow PEP8, except that the allowed line length is 100
  characters, not 80. But try to wrap comments after 80 characters.
* There are inconsistencies in the current code base since it’s a few years old
  already. New code should follow the current rules, however.
* At the moment, no automatic code formatting is done, but one idea might be to
  switch to the `black <https://black.readthedocs.io/>`_ code formatter at some
  point. If you’re familiar with its style, you can use that already now for
  new code to make the diff smaller.
* Prefer double quotation marks in new code. This will also make the diff smaller
  if we eventually switch to black.
* Using an IDE is beneficial (PyCharm, for example). It helps to catch lots of
  style issues early (unused imports, spacing etc.).
* Avoid unnecessary abbreviations for variable names. Code is more often read
  than written.
* When writing a help text for a new command-line option, look at the output of
  ``cutadapt --help`` and try to make it look nice and short.
* In comments and documentation, capitalize FASTQ, BWA, CPU etc.

