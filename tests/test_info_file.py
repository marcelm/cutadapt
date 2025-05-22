import pytest

from cutadapt.cli import main
from utils import assert_files_equal, datapath, cutpath

# pytest.mark.timeout will not fail even if pytest-timeout is not installed
try:
    import pytest_timeout as _unused
except ImportError:  # pragma: no cover
    raise ImportError("pytest_timeout needs to be installed")
del _unused


def test_info_file(run, tmp_path, cores):
    # The true adapter sequence in the illumina.fastq.gz data set is
    # GCCTAACTTCTTAGACTGCCTTAAGGACGT (fourth base is different from the sequence shown here)
    info_path = tmp_path / "info.txt"
    run(
        [
            "--cores",
            str(cores),
            "--info-file",
            info_path,
            "-a",
            "adapt=GCCGAACTTCTTAGACTGCCTTAAGGACGT",
        ],
        "illumina.fastq",
        "illumina.fastq.gz",
    )
    assert_files_equal(
        cutpath("illumina.info.txt"), info_path, ignore_trailing_space=True
    )


def test_info_file_times(run, tmp_path, cores):
    info_path = tmp_path / "info.txt"
    run(
        [
            "--cores",
            str(cores),
            "--info-file",
            info_path,
            "--times",
            "2",
            "-a",
            "adapt=GCCGAACTTCTTA",
            "-a",
            "adapt2=GACTGCCTTAAGGACGT",
        ],
        "illumina5.fastq",
        "illumina5.fastq",
    )
    assert_files_equal(
        cutpath("illumina5.info.txt"), info_path, ignore_trailing_space=True
    )


def test_info_file_fasta(run, tmp_path, cores):
    info_path = tmp_path / "info.txt"
    # Just make sure that it runs
    run(
        [
            "--cores",
            str(cores),
            "--info-file",
            info_path,
            "-a",
            "TTAGACATAT",
            "-g",
            "GAGATTGCCA",
            "--no-indels",
        ],
        "no_indels.fasta",
        "no_indels.fasta",
    )


def test_info_file_revcomp(run, tmp_path):
    info_path = tmp_path / "info-rc.txt"
    main(
        [
            "--info-file",
            str(info_path),
            "-a",
            "adapt=GAGTCG",
            "--revcomp",
            "--rename={header}",
            "-o",
            str(tmp_path / "out.fasta"),
            datapath("info-rc.fasta"),
        ]
    )
    assert_files_equal(cutpath("info-rc.txt"), info_path)


def test_reverse_complement_and_info_file(run, tmp_path, cores):
    info_path = str(tmp_path / "info.txt")
    run(
        [
            "--revcomp",
            "--no-index",
            "-g",
            "^TTATTTGTCT",
            "-g",
            "^TCCGCACTGG",
            "--info-file",
            info_path,
        ],
        "revcomp-single-normalize.fastq",
        "revcomp.1.fastq",
    )
    with open(info_path) as f:
        lines = f.readlines()
    assert len(lines) == 6
    assert lines[0].split("\t")[0] == "read1/1"
    assert lines[1].split("\t")[0] == "read2/1 rc"


def test_linked_info_file(tmp_path):
    info_path = tmp_path / "info.txt"
    main(
        [
            "-a linkedadapter=^AAAAAAAAAA...TTTTTTTTTT",
            "--info-file",
            str(info_path),
            "-o",
            str(tmp_path / "out.fasta"),
            datapath("linked.fasta"),
        ]
    )
    assert_files_equal(
        cutpath("linked-info.txt"), info_path, ignore_trailing_space=True
    )


def test_info_file_paired_only_r1(tmp_path):
    info_path = tmp_path / "info.txt"
    params = [
        "--info-file",
        str(info_path),
        "-o",
        str(tmp_path / "out.1.fastq"),
        "-p",
        str(tmp_path / "out.2.fastq"),
        datapath("paired.1.fastq"),
        datapath("paired.2.fastq"),
    ]
    main(params)


@pytest.mark.timeout(0.5)
def test_issue_296(tmp_path):
    # Hang when using both --no-trim and --info-file together
    info_path = tmp_path / "info.txt"
    reads_path = tmp_path / "reads.fasta"
    out_path = tmp_path / "out.fasta"
    reads_path.write_text(">read\nCACAAA\n")
    main(
        [
            "--info-file",
            str(info_path),
            "--no-trim",
            "-g",
            "TTTCAC",
            "-o",
            str(out_path),
            str(reads_path),
        ]
    )
    # Output should be unchanged because of --no-trim
    assert_files_equal(reads_path, out_path)
