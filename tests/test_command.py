"""Tests that run the program in a subprocess"""

import subprocess
import sys

from utils import datapath, assert_files_equal, cutpath


def test_run_cutadapt_process():
    subprocess.check_call(["cutadapt", "--version"])


def test_run_as_module():
    """Check that "python3 -m cutadapt ..." works"""
    from cutadapt import __version__
    with subprocess.Popen([sys.executable, "-m", "cutadapt", "--version"], stdout=subprocess.PIPE) as py:
        assert py.communicate()[0].decode().strip() == __version__


def test_standard_input_pipe(tmpdir, cores):
    """Read FASTQ from standard input"""
    out_path = str(tmpdir.join("out.fastq"))
    in_path = datapath("small.fastq")
    # Use 'cat' to simulate that no file name is available for stdin
    with subprocess.Popen(["cat", in_path], stdout=subprocess.PIPE) as cat:
        with subprocess.Popen([
            sys.executable, "-m", "cutadapt", "--cores", str(cores),
            "-a", "TTAGACATATCTCCGTCG", "-o", out_path, "-"],
            stdin=cat.stdout
        ) as py:
            _ = py.communicate()
            cat.stdout.close()
            _ = py.communicate()[0]
    assert_files_equal(cutpath("small.fastq"), out_path)


def test_standard_output(tmpdir, cores):
    """Write FASTQ to standard output (not using --output/-o option)"""
    out_path = str(tmpdir.join("out.fastq"))
    with open(out_path, "w") as out_file:
        py = subprocess.Popen([
            sys.executable, "-m", "cutadapt", "--cores", str(cores),
            "-a", "TTAGACATATCTCCGTCG", datapath("small.fastq")],
            stdout=out_file)
        _ = py.communicate()
    assert_files_equal(cutpath("small.fastq"), out_path)


def test_explicit_standard_output(tmpdir, cores):
    """Write FASTQ to standard output (using "-o -")"""

    out_path = str(tmpdir.join("out.fastq"))
    with open(out_path, "w") as out_file:
        py = subprocess.Popen([
            sys.executable, "-m", "cutadapt", "-o", "-", "--cores", str(cores),
            "-a", "TTAGACATATCTCCGTCG", datapath("small.fastq")],
            stdout=out_file)
        _ = py.communicate()
    assert_files_equal(cutpath("small.fastq"), out_path)


def test_force_fasta_output(tmpdir, cores):
    """Write FASTA to standard output even on FASTQ input"""

    out_path = str(tmpdir.join("out.fasta"))
    with open(out_path, "w") as out_file:
        py = subprocess.Popen([
            sys.executable, "-m", "cutadapt", "--fasta", "-o", "-", "--cores", str(cores),
            "-a", "TTAGACATATCTCCGTCG", datapath("small.fastq")],
            stdout=out_file)
        _ = py.communicate()
    assert_files_equal(cutpath("small.fasta"), out_path)


def test_non_utf8_locale():
    subprocess.check_call(
        [sys.executable, "-m", "cutadapt", "-o", "/dev/null", datapath("small.fastq")],
        env={"LC_CTYPE": "C"},
    )
