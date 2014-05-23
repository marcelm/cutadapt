"""
Build cutadapt.

If pre-generated .c extension files are found, then Cython is not run, even if it is installed, unless --cython is given
on the command line.

If .c files are not found, then Cython is always run. This is the case, for example, when using a fresh Git checkout.
"""
import sys
import os.path
from distutils.core import setup, Extension
from distutils.version import LooseVersion

from cutadapt import __version__

MIN_CYTHON_VERSION = '0.15'

if sys.version_info < (2, 6):
	sys.stdout.write("At least Python 2.6 is required.\n")
	sys.exit(1)


use_cython = not os.path.exists('cutadapt/_align.c')

if '--cython' in sys.argv:
	use_cython = True
	sys.argv.remove('--cython')

# Try to find out whether a recent enough Cython is installed.
# If it is not, fall back to using the pre-compiled C sources.
# Pre-compiled sources are available only in the official releases,
# not within the Git repository.
if use_cython:
	try:
		from Cython import __version__ as cyversion
	except ImportError:
		sys.stdout.write(
			"ERROR: Cython is not installed. Install at least Cython >= " + str(MIN_CYTHON_VERSION) +
		    " to continue.\n")
		sys.exit(1)
	if LooseVersion(cyversion) < LooseVersion(MIN_CYTHON_VERSION):
		sys.stdout.write(
			"Error: Your Cython is at version '" + str(cyversion) +
			"', but at least version " + str(MIN_CYTHON_VERSION) + " is required.\n")
		sys.exit(1)

	from Cython.Distutils import build_ext
	cmdclass = {'build_ext': build_ext}
else:
	cmdclass = {}

ext = '.pyx' if use_cython else '.c'

extensions = [
	Extension('cutadapt._align', sources=['cutadapt/_align' + ext]),
	Extension('cutadapt._qualtrim', sources=['cutadapt/_qualtrim' + ext]),
	Extension('cutadapt._seqio', sources=['cutadapt/_seqio' + ext]),
]

setup(
	name = 'cutadapt',
	version = __version__,
	author = 'Marcel Martin',
	author_email = 'marcel.martin@scilifelab.se',
	url = 'http://code.google.com/p/cutadapt/',
	description = 'trim adapters from high-throughput sequencing reads',
	license = 'MIT',
	ext_modules = extensions,
	cmdclass = cmdclass,
	packages = ['cutadapt', 'cutadapt.scripts'],
	scripts = ['bin/cutadapt'],
	classifiers = [
		"Development Status :: 5 - Production/Stable",
		"Environment :: Console",
		"Intended Audience :: Science/Research",
		"License :: OSI Approved :: MIT License",
		"Natural Language :: English",
		"Programming Language :: Cython",
		"Programming Language :: Python :: 2.6",
		"Programming Language :: Python :: 2.7",
		"Programming Language :: Python :: 3",
		"Topic :: Scientific/Engineering :: Bio-Informatics"
	]
)
