"""
Build Cutadapt.
"""
import sys

from setuptools import setup, Extension, find_packages
from distutils.version import LooseVersion
from Cython.Build import cythonize

MIN_CYTHON_VERSION = '0.28'


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


encoding_arg = {'encoding': 'utf-8'} if sys.version_info[0] >= 3 else dict()
with open('README.rst', **encoding_arg) as f:
    long_description = f.read()

extensions = [
    Extension('cutadapt._align', sources=['src/cutadapt/_align.pyx']),
    Extension('cutadapt.qualtrim', sources=['src/cutadapt/qualtrim.pyx']),
]

setup(
    name='cutadapt',
    use_scm_version={'write_to': 'src/cutadapt/_version.py'},
    author='Marcel Martin',
    author_email='marcel.martin@scilifelab.se',
    url='https://cutadapt.readthedocs.io/',
    description='trim adapters from high-throughput sequencing reads',
    long_description=long_description,
    long_description_content_type='text/x-rst',
    license='MIT',
    ext_modules=cythonize(extensions),
    package_dir={'': 'src'},
    packages=find_packages('src'),
    entry_points={'console_scripts': ['cutadapt = cutadapt.__main__:main_cli']},
    install_requires=[
        'dnaio~=0.5',
        'xopen~=1.1',
        "dataclasses>=0.8; python_version<'3.7'",
    ],
    extras_require={
        'dev': ['Cython', 'pytest', 'pytest-timeout', 'pytest-mock', 'sphinx', 'sphinx_issues'],
    },
    python_requires='>=3.6',
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
    project_urls={
        "Changelog": "https://cutadapt.readthedocs.io/en/stable/changes.html",
    },
)
