import os
from cutadapt.files import ProxyTextFile


def test_proxy_test_file():
    newline = os.linesep.encode()
    pf = ProxyTextFile()
    print("hello", file=pf)
    assert pf.drain() == b"hello" + newline
    assert pf.drain() == b""

    print("world", file=pf, end="\n")
    print("foo", file=pf, end="\n")
    assert pf.drain() == b"world" + newline + b"foo" + newline


def test_proxy_test_file_pickleable():
    import pickle

    pf = ProxyTextFile()
    pickled = pickle.dumps(pf)

    unpickled = pickle.loads(pickled)
    assert isinstance(unpickled, ProxyTextFile)
