from distutils.core import setup, Extension
import sys

from cutadapt import __version__

if sys.version_info < (2, 6):
	sys.stdout.write("At least Python 2.6 is required.\n")
	sys.exit(1)

try:
	from Cython.Distutils import build_ext
except ImportError:
	# no Cython available
	cmdclass = { }
	align_module = Extension('cutadapt.calign', sources = [ 'cutadapt/calign.c' ])
else:
	cmdclass = { 'build_ext' : build_ext }
	align_module = Extension('cutadapt.calign', sources = [ 'cutadapt/calign.pyx' ])


setup(
	name = 'cutadapt',
	version = __version__,
	author = 'Marcel Martin',
	author_email = 'marcel.martin@tu-dortmund.de',
	url = 'http://code.google.com/p/cutadapt/',
	description = 'trim adapters from high-throughput sequencing reads',
	license = 'MIT',
	ext_modules = [align_module],
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

