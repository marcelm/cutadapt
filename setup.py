from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

#calign = Extension('calign', sources = [ 'calignmodule.c' ])
pyxalign = Extension("align", ["align.pyx"])

setup(name = 'cutadapt',
	version = '0.5.1',
	description = 'trim adapters from high-throughput sequencing reads',
	author = 'Marcel Martin',
	author_email = 'marcel.martin@tu-dortmund.de',
	url = 'http://cutadapt.googlecode.com/',
	license = 'MIT',
	cmdclass = {'build_ext': build_ext},
	ext_modules = [ pyxalign ],
	py_modules = ['fasta'],
	scripts = ['cutadapt']
)
