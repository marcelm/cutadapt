from distutils.core import setup, Extension
import sys
import os.path

from cutadapt import __version__

if sys.version_info < (2, 6):
	sys.stdout.write("At least Python 2.6 is required.\n")
	sys.exit(1)

try:
	from Cython.Distutils import build_ext
	USE_CYTHON = True
except ImportError:
	USE_CYTHON = False
	cmdclass = { }
else:
	cmdclass = { 'build_ext' : build_ext }

ext = '.pyx' if USE_CYTHON else '.c'

extensions = [
	Extension('cutadapt._align', sources=['cutadapt/_align' + ext]),
	Extension('cutadapt._qualtrim', sources=['cutadapt/_qualtrim' + ext]),
]

if not USE_CYTHON:
	for extension in extensions:
		if not os.path.exists(extension.sources[0]):
			sys.stdout.write("Cython is not installed and pre-compiled C sources\n")
			sys.stdout.write("are also not available. You need to install Cython to continue.\n")
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
		"Programming Language :: Python :: 3.2",
		"Programming Language :: Python :: 3.3",
		"Topic :: Scientific/Engineering :: Bio-Informatics"
	]
)
