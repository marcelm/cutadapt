"""
Build cutadapt.

Cython is run when
* no pre-generated C sources are found,
* or the pre-generated C sources are out of date,
* or when --cython is given on the command line.
"""
import sys
import os.path

from distutils.core import setup, Extension
from distutils.version import LooseVersion
from distutils.command.sdist import sdist as _sdist
from distutils.command.build_ext import build_ext as _build_ext

MIN_CYTHON_VERSION = '0.24'

if sys.version_info < (2, 6):
	sys.stdout.write("At least Python 2.6 is required.\n")
	sys.exit(1)


# set __version__
with open(os.path.join(os.path.dirname(__file__), 'cutadapt', '__init__.py')) as f:
	for line in f:
		if line.startswith('__version__'):
			exec(line)
			break


def out_of_date(extensions):
	"""
	Check whether any pyx source is newer than the corresponding generated
	C source or whether any C source is missing.
	"""
	for extension in extensions:
		for pyx in extension.sources:
			path, ext = os.path.splitext(pyx)
			if ext not in ('.pyx', '.py'):
				continue
			if extension.language == 'c++':
				csource = path + '.cpp'
			else:
				csource = path + '.c'
			# When comparing modification times, allow five seconds slack:
			# If the installation is being run from pip, modification
			# times are not preserved and therefore depends on the order in
			# which files were unpacked.
			if not os.path.exists(csource) or (
				os.path.getmtime(pyx) > os.path.getmtime(csource) + 5):
				return True
	return False


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


class build_ext(_build_ext):
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
		_build_ext.run(self)


class sdist(_sdist):
	def run(self):
		# Make sure the compiled Cython files in the distribution are up-to-date
		from Cython.Build import cythonize
		check_cython_version()
		cythonize(extensions)
		_sdist.run(self)


setup(
	name = 'cutadapt',
	version = __version__,
	author = 'Marcel Martin',
	author_email = 'marcel.martin@scilifelab.se',
	url = 'https://cutadapt.readthedocs.org/',
	description = 'trim adapters from high-throughput sequencing reads',
	license = 'MIT',
	cmdclass = {'sdist': sdist, 'build_ext': build_ext},
	ext_modules = extensions,
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
