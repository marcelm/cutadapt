"""Tests that run the program in a subprocess"""

import subprocess
import sys
import os

import pytest

from utils import datapath, assert_files_equal, cutpath


def test_run_cutadapt_process():
    subprocess.check_call(["cutadapt", "--version"])


def test_run_as_module():
    """Check that "python3 -m cutadapt ..." works"""
    from cutadapt import __version__

    with subprocess.Popen(
        [sys.executable, "-m", "cutadapt", "--version"], stdout=subprocess.PIPE
    ) as py:
        assert py.communicate()[0].decode().strip() == __version__


@pytest.mark.skipif(sys.platform == "win32", reason="Perhaps this can be fixed")
def test_standard_input_pipe(tmp_path, cores):
    """Read FASTQ from standard input"""
    out_path = os.fspath(tmp_path / "out.fastq")
    in_path = datapath("small.fastq")
    # Simulate that no file name is available for stdin
    with subprocess.Popen(["cat", in_path], stdout=subprocess.PIPE) as cat:
        with subprocess.Popen(
            [
                sys.executable,
                "-m",
                "cutadapt",
                "--cores",
                str(cores),
                "-a",
                "TTAGACATATCTCCGTCG",
                "-o",
                out_path,
                "-",
            ],
            stdin=cat.stdout,
        ) as py:
            _ = py.communicate()
            cat.stdout.close()
            _ = py.communicate()[0]
    assert_files_equal(cutpath("small.fastq"), out_path)


def test_standard_output(tmp_path, cores):
    """Write FASTQ to standard output (not using --output/-o option)"""
    out_path = os.fspath(tmp_path / "out.fastq")
    with open(out_path, "w") as out_file:
        py = subprocess.Popen(
            [
                sys.executable,
                "-m",
                "cutadapt",
                "--cores",
                str(cores),
                "-a",
                "TTAGACATATCTCCGTCG",
                datapath("small.fastq"),
            ],
            stdout=out_file,
        )
        _ = py.communicate()
    assert_files_equal(cutpath("small.fastq"), out_path)


def test_explicit_standard_output(tmp_path, cores):
    """Write FASTQ to standard output (using "-o -")"""

    out_path = os.fspath(tmp_path / "out.fastq")
    with open(out_path, "w") as out_file:
        py = subprocess.Popen(
            [
                sys.executable,
                "-m",
                "cutadapt",
                "-o",
                "-",
                "--cores",
                str(cores),
                "-a",
                "TTAGACATATCTCCGTCG",
                datapath("small.fastq"),
            ],
            stdout=out_file,
        )
        _ = py.communicate()
    assert_files_equal(cutpath("small.fastq"), out_path)


def test_force_fasta_output(tmp_path, cores):
    """Write FASTA to standard output even on FASTQ input"""

    out_path = os.fspath(tmp_path / "out.fasta")
    with open(out_path, "w") as out_file:
        py = subprocess.Popen(
            [
                sys.executable,
                "-m",
                "cutadapt",
                "--fasta",
                "-o",
                "-",
                "--cores",
                str(cores),
                "-a",
                "TTAGACATATCTCCGTCG",
                datapath("small.fastq"),
            ],
            stdout=out_file,
        )
        _ = py.communicate()
    assert_files_equal(cutpath("small.fasta"), out_path)


@pytest.mark.skipif(sys.platform == "win32", reason="Maybe this can be made to work")
def test_non_utf8_locale():
    subprocess.check_call(
        [sys.executable, "-m", "cutadapt", "-o", os.devnull, datapath("small.fastq")],
        env={"LC_CTYPE": "C"},
    )
