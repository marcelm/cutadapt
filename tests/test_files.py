import os
import pickle

from cutadapt.files import ProxyTextFile, ProxyRecordWriter, OutputFiles,open_rb
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


class TestOutputFiles:
    def test_open_text(self, tmp_path):
        o = OutputFiles(
            proxied=False,
            qualities=False,
            interleaved=False,
        )
        path = tmp_path / "out.txt"
        f = o.open_text(path)
        print("Hello", file=f)
        o.close()
        assert path.read_text() == "Hello\n"

    def test_open_record_writer(self, tmp_path):
        o = OutputFiles(
            proxied=False,
            qualities=True,
            interleaved=False,
        )
        path = tmp_path / "out.fastq"
        f = o.open_record_writer(path)
        f.write(SequenceRecord("r", "ACGT", "####"))
        o.close()
        assert path.read_text() == "@r\nACGT\n+\n####\n"

    def test_paired_record_writer(self, tmp_path):
        o = OutputFiles(
            proxied=False,
            qualities=True,
            interleaved=False,
        )
        path1 = tmp_path / "out.1.fastq"
        path2 = tmp_path / "out.2.fastq"
        f = o.open_record_writer(path1, path2)
        f.write(
            SequenceRecord("r", "AACC", "####"), SequenceRecord("r", "GGTT", "####")
        )
        o.close()
        assert path1.read_text() == "@r\nAACC\n+\n####\n"
        assert path2.read_text() == "@r\nGGTT\n+\n####\n"

    def test_interleaved_record_writer(self, tmp_path):
        o = OutputFiles(
            proxied=False,
            qualities=True,
            interleaved=True,
        )
        path = tmp_path / "out.1.fastq"
        f = o.open_record_writer(path, interleaved=True)
        f.write(
            SequenceRecord("r", "AACC", "####"), SequenceRecord("r", "GGTT", "####")
        )
        o.close()
        assert path.read_text() == "@r\nAACC\n+\n####\n@r\nGGTT\n+\n####\n"


def test_open_rb_local_file(tmp_path):
        
    # Create a local file
    file_path = tmp_path / "test.txt"
    file_path.write_text("Hello, World!")
    # Test opening a local file
    file = open_rb(str(file_path))
    assert file.read() == b"Hello, World!"
    file.close()

def test_open_rb_remote_file():
    # Test opening a remote file over https
    file = open_rb("https://raw.githubusercontent.com/marcelm/cutadapt/main/tests/data/454.fa")
    assert file.readline() == b">000163_1255_2627 length=52 uaccno=E0R4ISW01DCIQD\n"
    file.close()

def test_open_rb_s3_file():
    # Test opening a remote file on s3 
    file = open_rb("s3://platinum-genomes/2017-1.0/md5sum.txt")
    assert file.readline() == b"2e6aa26b42283bbbc4ca03686f427dc2  ./hg38/small_variants/ConfidentRegions.bed.gz\n"
    file.close()
    # - test force fasta
    # - test qualities
    # - test proxied
    # - test complaint about duplicate file names
