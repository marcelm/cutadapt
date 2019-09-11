import os.path
import subprocess
import sys
from contextlib import contextmanager


@contextmanager
def redirect_stderr():
    """Send stderr to stdout. Nose doesn't capture stderr, yet."""
    old_stderr = sys.stderr
    sys.stderr = sys.stdout
    yield
    sys.stderr = old_stderr


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


def binomial(n, k):
    """
    Return binomial coefficient ('n choose k').
    This implementation does not use factorials.
    """
    k = min(k, n - k)
    if k < 0:
        return 0
    r = 1
    for j in range(k):
        r *= n - j
        r //= j + 1
    return r
