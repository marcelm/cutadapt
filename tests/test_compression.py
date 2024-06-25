import pytest

from cutadapt.cli import main
from utils import assert_files_equal, datapath, cutpath


def test_gz_multiblock(run):
    """compressed gz file with multiple blocks (created by concatenating two .gz files)"""
    run("-b TTAGACATATCTCCGTCG", "small.fastq", "multiblock.fastq.gz")


def test_extensiontxtgz(run):
    """automatic recognition of "_sequence.txt.gz" extension"""
    run("-b TTAGACATATCTCCGTCG", "s_1_sequence.txt", "s_1_sequence.txt.gz")


def test_compressed_fasta(run):
    run("", "simple.fasta", "simple.fasta.gz")


def test_write_compressed_fastq(cores, tmp_path):
    main(
        [
            "--cores",
            str(cores),
            "-o",
            str(tmp_path / "out.fastq.gz"),
            datapath("small.fastq"),
        ]
    )


@pytest.mark.parametrize("extension", ["bz2", "xz", "gz"])
def test_compressed_output(tmp_path, cores, extension):
    out_path = str(tmp_path / ("small.fastq." + extension))
    params = [
        "--cores",
        str(cores),
        "-a",
        "TTAGACATATCTCCGTCG",
        "-o",
        out_path,
        datapath("small.fastq"),
    ]
    main(params)


def test_bzip2_multiblock(run):
    run("-b TTAGACATATCTCCGTCG", "small.fastq", "multiblock.fastq.bz2")


def test_xz(run):
    run("-b TTAGACATATCTCCGTCG", "small.fastq", "small.fastq.xz")


def test_bzip2_input(run, cores):
    run(
        ["--cores", str(cores), "-a", "TTAGACATATCTCCGTCG"],
        "small.fastq",
        "small.fastq.bz2",
    )


def test_underscore(run):
    """File name ending in _fastq.gz (issue #275)"""
    run("-b TTAGACATATCTCCGTCG", "small.fastq", "underscore_fastq.gz")
