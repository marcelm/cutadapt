import os
import pickle

from cutadapt.files import ProxyTextFile, ProxyRecordWriter
from dnaio import SequenceRecord


def test_proxy_text_file():
    newline = os.linesep.encode()
    pf = ProxyTextFile()
    print("hello", file=pf)
    assert pf.drain() == [b"hello" + newline]
    assert pf.drain() == [b""]

    print("world", file=pf, end="\n")
    print("foo", file=pf, end="\n")
    assert pf.drain() == [b"world" + newline + b"foo" + newline]


def test_proxy_test_file_pickleable():
    pf = ProxyTextFile()
    pickled = pickle.dumps(pf)

    unpickled = pickle.loads(pickled)
    assert isinstance(unpickled, ProxyTextFile)


def test_proxy_record_writer():
    pw = ProxyRecordWriter(n_files=1, qualities=True)
    pw.write(SequenceRecord("name", "ACGT", qualities="####"))
    assert pw.drain() == [
        b"@name\nACGT\n+\n####\n",
    ]

    pw.write(SequenceRecord("foo", "AA", "HH"))
    pw.write(SequenceRecord("bar", "CC", ",,"))
    assert pw.drain() == [
        b"@foo\nAA\n+\nHH\n@bar\nCC\n+\n,,\n",
    ]


def test_proxy_record_writer_paired():
    pw = ProxyRecordWriter(n_files=2, qualities=True)
    pw.write(
        SequenceRecord("name", "ACGT", qualities="####"),
        SequenceRecord("name", "GGGG", qualities="!!!!"),
    )
    assert pw.drain() == [b"@name\nACGT\n+\n####\n", b"@name\nGGGG\n+\n!!!!\n"]

    pw.write(
        SequenceRecord("foo", "AA", "HH"),
        SequenceRecord("foo", "TT", "33"),
    )
    pw.write(
        SequenceRecord("bar", "CC", ",,"),
        SequenceRecord("bar", "GGG", "444"),
    )
    assert pw.drain() == [
        b"@foo\nAA\n+\nHH\n@bar\nCC\n+\n,,\n",
        b"@foo\nTT\n+\n33\n@bar\nGGG\n+\n444\n",
    ]


def test_proxy_record_writer_picklable():
    pw = ProxyRecordWriter(n_files=2, qualities=True)
    pickled = pickle.dumps(pw)

    unpickled = pickle.loads(pickled)
    assert isinstance(unpickled, ProxyRecordWriter)
    assert unpickled._n_files == 2
