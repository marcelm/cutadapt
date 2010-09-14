from distutils.core import setup, Extension

module = Extension('calign', sources = [ 'calignmodule.c' ])

setup (name = 'cutadapt',
	version = '0.1',
	description = 'trim adapters from high-throughput sequencing reads',
	ext_modules = [ module ],
	py_modules = ['fasta' ],
	scripts = ['cutadapt' ]
)
