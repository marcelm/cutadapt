"""
Build Cutadapt.
"""
import sys

from setuptools import setup, Extension, find_packages
from Cython.Build import cythonize

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
