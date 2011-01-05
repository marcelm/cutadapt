from distutils.core import setup, Extension
import sys

if sys.version_info < (2, 6):
	print "At least Python 2.6 is required."
	sys.exit(1)

module = Extension('calign', sources = [ 'calignmodule.c' ])

setup(name = 'cutadapt',
	version = '0.9',
	description = 'trim adapters from high-throughput sequencing reads',
	author = 'Marcel Martin',
	author_email = 'marcel.martin@tu-dortmund.de',
	url = 'http://cutadapt.googlecode.com/',
	license = 'MIT',
	ext_modules = [ module ],
	py_modules = ['fasta', 'align'],
	scripts = ['cutadapt']
)
