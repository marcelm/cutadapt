"""
Build Cutadapt.
"""
import sys
import os.path

from setuptools import setup, Extension, find_packages
from distutils.version import LooseVersion
from distutils.command.sdist import sdist as _sdist
from distutils.command.build_ext import build_ext as _build_ext

MIN_CYTHON_VERSION = '0.28'

if sys.version_info[:2] < (3, 4):
    sys.stdout.write('You need at least Python 3.4\n')
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
    Extension('cutadapt._align', sources=['src/cutadapt/_align.pyx']),
    Extension('cutadapt.qualtrim', sources=['src/cutadapt/qualtrim.pyx']),
]


class BuildExt(_build_ext):
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
        super().run()


class SDist(_sdist):
    def run(self):
        # Make sure the compiled Cython files in the distribution are up-to-date
        from Cython.Build import cythonize
        check_cython_version()
        cythonize(extensions)
        super().run()


encoding_arg = {'encoding': 'utf-8'} if sys.version > '3' else dict()
with open('README.rst', **encoding_arg) as f:
    long_description = f.read()

setup(
    name='cutadapt',
    setup_requires=['setuptools_scm'],  # Support pip versions that don't know about pyproject.toml
    use_scm_version={'write_to': 'src/cutadapt/_version.py'},
    author='Marcel Martin',
    author_email='marcel.martin@scilifelab.se',
    url='https://cutadapt.readthedocs.io/',
    description='trim adapters from high-throughput sequencing reads',
    long_description=long_description,
    license='MIT',
    cmdclass={'build_ext': BuildExt, 'sdist': SDist},
    ext_modules=extensions,
    package_dir={'': 'src'},
    packages=find_packages('src'),
    entry_points={'console_scripts': ['cutadapt = cutadapt.__main__:main']},
    install_requires=['dnaio>=0.3', 'xopen>=0.5.0'],
    extras_require={
        'dev': ['Cython', 'pytest', 'pytest-timeout', 'sphinx', 'sphinx_issues'],
    },
    python_requires='>=3.4',
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Cython",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
)
