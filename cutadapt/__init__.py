# coding: utf-8
from __future__ import print_function, division, absolute_import
import sys

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions


def check_importability():  # pragma: no cover
	try:
		import cutadapt._align
	except ImportError as e:
		if 'undefined symbol' in str(e):
			print("""
ERROR: A required extension module could not be imported because it is
incompatible with your system. A quick fix is to recompile the extension
modules with the following command:

    {0} setup.py build_ext -i

See the documentation for alternative ways of installing the program.

The original error message follows.
""".format(sys.executable))
		raise
