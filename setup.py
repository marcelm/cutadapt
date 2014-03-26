import sys
import os.path
from distutils.core import setup, Extension
from distutils.version import StrictVersion

from cutadapt import __version__

MIN_CYTHON_VERSION = '0.15'

if sys.version_info < (2, 6):
	sys.stdout.write("At least Python 2.6 is required.\n")
	sys.exit(1)

# Try to find out whether a recent enough Cython is installed.
# If it is not, fall back to using the pre-compiled C sources.
# Pre-compiled sources are available only in the official releases,
# not within the Git repository.
try:
	from Cython import __version__ as cyversion
	USE_CYTHON = StrictVersion(cyversion) >= StrictVersion(MIN_CYTHON_VERSION)
except ImportError:
	cyversion = None
	USE_CYTHON = False

if USE_CYTHON:
	from Cython.Distutils import build_ext
	cmdclass = { 'build_ext' : build_ext }
else:
	cmdclass = { }

ext = '.pyx' if USE_CYTHON else '.c'

extensions = [
	Extension('cutadapt._align', sources=['cutadapt/_align' + ext]),
	Extension('cutadapt._qualtrim', sources=['cutadapt/_qualtrim' + ext]),
	Extension('cutadapt._seqio', sources=['cutadapt/_seqio' + ext]),
]

if not USE_CYTHON:
	for extension in extensions:
		if not os.path.exists(extension.sources[0]):
			sys.stdout.write("Error:\n")
			sys.stdout.write("Cython version " + str(MIN_CYTHON_VERSION) + " or later is not installed and pre-compiled C sources\n")
			sys.stdout.write("are also not available. You need to install Cython to continue.\n")
			if cyversion:
				sys.stdout.write("Found Cython version: " + cyversion + "\n")
			sys.exit(1)

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
