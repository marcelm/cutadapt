import os.path
import subprocess
import sys
from contextlib import contextmanager
from shutil import rmtree
from tempfile import mkdtemp

from cutadapt.__main__ import main


@contextmanager
def redirect_stderr():
    """Send stderr to stdout. Nose doesn't capture stderr, yet."""
    old_stderr = sys.stderr
    sys.stderr = sys.stdout
    yield
    sys.stderr = old_stderr


@contextmanager
def temporary_path(name):
    tempdir = mkdtemp(prefix='cutadapt-tests.')
    path = os.path.join(tempdir, name)
    try:
        yield path
    finally:
        rmtree(tempdir)


def datapath(path):
    return os.path.join(os.path.dirname(__file__), 'data', path)


def cutpath(path):
    return os.path.join(os.path.dirname(__file__), 'cut', path)


class FilesDifferent(Exception):
    pass


def assert_files_equal(path1, path2):
    try:
        subprocess.check_output(['diff', '-u', path1, path2], stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise FilesDifferent('\n' + e.output.decode()) from None


def run(params, expected, inpath, inpath2=None):
    if type(params) is str:
        params = params.split()
    with temporary_path(expected) as tmp_fastaq:
        params += ['-o', tmp_fastaq]  # TODO not parallelizable
        params += [datapath(inpath)]
        if inpath2:
            params += [datapath(inpath2)]
        assert main(params) is None
        # TODO redirect standard output
        assert_files_equal(cutpath(expected), tmp_fastaq)
    # TODO diff log files
