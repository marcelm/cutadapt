"""
Build cutadapt.
"""
import sys
import os.path

from setuptools import setup, Extension
from distutils.version import LooseVersion
from distutils.command.sdist import sdist as _sdist
from distutils.command.build_ext import build_ext as _build_ext
import versioneer

MIN_CYTHON_VERSION = '0.24'

if sys.version_info < (2, 6):
	sys.stdout.write("At least Python 2.6 is required.\n")
	sys.exit(1)


def no_cythonize(extensions, **_ignore):
	"""
	Change file extensions from .pyx to .c or .cpp.

	Copied from Cython documentation
	"""
	for extension in extensions:
		sources = []
		for sfile in extension.sources:
			path, ext = os.path.splitext(sfile)
			if ext in ('.pyx', '.py'):
				if extension.language == 'c++':
					ext = '.cpp'
				else:
					ext = '.c'
				sfile = path + ext
			sources.append(sfile)
		extension.sources[:] = sources


def check_cython_version():
	"""Exit if Cython was not found or is too old"""
	try:
		from Cython import __version__ as cyversion
	except ImportError:
		sys.stdout.write(
			"ERROR: Cython is not installed. Install at least Cython version " +
			str(MIN_CYTHON_VERSION) + " to continue.\n")
		sys.exit(1)
	if LooseVersion(cyversion) < LooseVersion(MIN_CYTHON_VERSION):
		sys.stdout.write(
			"ERROR: Your Cython is at version '" + str(cyversion) +
			"', but at least version " + str(MIN_CYTHON_VERSION) + " is required.\n")
		sys.exit(1)


extensions = [
	Extension('cutadapt._align', sources=['cutadapt/_align.pyx']),
	Extension('cutadapt._qualtrim', sources=['cutadapt/_qualtrim.pyx']),
	Extension('cutadapt._seqio', sources=['cutadapt/_seqio.pyx']),
]

cmdclass = versioneer.get_cmdclass()
versioneer_build_ext = cmdclass.get('build_ext', _build_ext)
versioneer_sdist = cmdclass.get('sdist', _sdist)


class build_ext(versioneer_build_ext):
	def run(self):
		# If we encounter a PKG-INFO file, then this is likely a .tar.gz/.zip
		# file retrieved from PyPI that already includes the pre-cythonized
		# extension modules, and then we do not need to run cythonize().
		if os.path.exists('PKG-INFO'):
			no_cythonize(extensions)
		else:
			# Otherwise, this is a 'developer copy' of the code, and then the
			# only sensible thing is to require Cython to be installed.
			check_cython_version()
			from Cython.Build import cythonize
			self.extensions = cythonize(self.extensions)
		versioneer_build_ext.run(self)


class sdist(versioneer_sdist):
	def run(self):
		# Make sure the compiled Cython files in the distribution are up-to-date
		from Cython.Build import cythonize
		check_cython_version()
		cythonize(extensions)
		versioneer_sdist.run(self)


cmdclass['build_ext'] = build_ext
cmdclass['sdist'] = sdist


encoding_arg = {'encoding': 'utf-8'} if sys.version > '3' else dict()
with open('README.rst', **encoding_arg) as f:
	long_description = f.read()

setup(
	name = 'cutadapt',
	version = versioneer.get_version(),
	author = 'Marcel Martin',
	author_email = 'marcel.martin@scilifelab.se',
	url = 'https://cutadapt.readthedocs.io/',
	description = 'trim adapters from high-throughput sequencing reads',
	long_description = long_description,
	license = 'MIT',
	cmdclass = cmdclass,
	ext_modules = extensions,
	packages = ['cutadapt', 'cutadapt.scripts'],
	install_requires = ['xopen>=0.1.1'],
	entry_points = {'console_scripts': ['cutadapt = cutadapt.scripts.cutadapt:main']},
	classifiers = [
		"Development Status :: 5 - Production/Stable",
		"Environment :: Console",
		"Intended Audience :: Science/Research",
		"License :: OSI Approved :: MIT License",
		"Natural Language :: English",
		"Programming Language :: Cython",
		"Programming Language :: Python :: 2.7",
		"Programming Language :: Python :: 3",
		"Topic :: Scientific/Engineering :: Bio-Informatics"
	]
)
