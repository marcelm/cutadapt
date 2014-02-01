from distutils.core import setup, Extension
import sys
import os.path

from cutadapt import __version__

if sys.version_info < (2, 6):
	sys.stdout.write("At least Python 2.6 is required.\n")
	sys.exit(1)

try:
	from Cython.Distutils import build_ext
except ImportError:
	# no Cython available
	cmdclass = { }
	align_sources = [ 'cutadapt/calign.c' ]
	qualtrim_sources = [ 'cutadapt/cqualtrim.c' ]
	if not os.path.exists(align_sources[0]) or not os.path.exists(qualtrim_sources[0]):
		sys.stdout.write("Cython is not installed and a pre-compiled alignment module\n")
		sys.stdout.write("is also not available. You need to install Cython to continue.\n")
		sys.exit(1)
else:
	cmdclass = { 'build_ext' : build_ext }
	align_sources = [ 'cutadapt/calign.pyx' ]
	qualtrim_sources = [ 'cutadapt/cqualtrim.pyx' ]

align_module = Extension('cutadapt.calign', sources=align_sources)
qualtrim_module = Extension('cutadapt.cqualtrim', sources=qualtrim_sources)

setup(
	name = 'cutadapt',
	version = __version__,
	author = 'Marcel Martin',
	author_email = 'marcel.martin@tu-dortmund.de',
	url = 'http://code.google.com/p/cutadapt/',
	description = 'trim adapters from high-throughput sequencing reads',
	license = 'MIT',
	ext_modules = [align_module, qualtrim_module],
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
