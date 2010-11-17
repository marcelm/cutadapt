from distutils.core import setup, Extension

module = Extension('calign', sources = [ 'calignmodule.c' ])

setup(name = 'cutadapt',
	version = '0.5'
	description = 'trim adapters from high-throughput sequencing reads',
	author = 'Marcel Martin',
	author_email = 'marcel.martin@tu-dortmund.de',
	url = 'http://cutadapt.googlecode.com/',
	license = 'MIT',
	ext_modules = [ module ],
	py_modules = ['fasta'],
	scripts = ['cutadapt']
)
