import os.path
import shutil
import tempfile
from itertools import product

import pytest

from cutadapt.__main__ import main
from utils import assert_files_equal, datapath, cutpath, redirect_stderr


@pytest.fixture
def run_paired(tmpdir):
    def _run(params, in1, in2, expected1, expected2, cores):
        if type(params) is str:
            params = params.split()
        params += ["--cores", str(cores), "--buffer-size=512"]
        path1 = str(tmpdir.join(expected1))
        path2 = str(tmpdir.join(expected2))
        params += ["-o", path1, "-p", path2]
        params += [datapath(in1), datapath(in2)]
        assert main(params) is None
        assert_files_equal(cutpath(expected1), path1)
        assert_files_equal(cutpath(expected2), path2)

    return _run


@pytest.fixture
def run_interleaved(tmpdir):
    """
    Interleaved input or output (or both)
    """
    def _run(params, inpath1, inpath2=None, expected1=None, expected2=None, cores=1):
        assert not (inpath1 and inpath2 and expected1 and expected2)
        assert not (expected2 and not expected1)
        assert not (inpath2 and not inpath1)
        if type(params) is str:
            params = params.split()
        params += ["--interleaved", "--cores", str(cores), "--buffer-size=512"]
        tmp1 = str(tmpdir.join("out1-" + expected1))
        params += ["-o", tmp1]
        paths = [datapath(inpath1)]
        if inpath2:
            paths += [datapath(inpath2)]
        if expected2:
            tmp2 = str(tmpdir.join("out2-" + expected2))
            params += ["-p", tmp2]
            assert main(params + paths) is None
            assert_files_equal(cutpath(expected2), tmp2)
        else:
            assert main(params + paths) is None
        assert_files_equal(cutpath(expected1), tmp1)

    return _run


def test_paired_end_no_legacy(run_paired, cores):
    """--paired-output, not using -A/-B/-G"""
    # the -m 14 filters out one read, which should then also be removed from the second file
    # Since legacy mode was removed, -q 10 should filter out an additional read which gets
    # quality-trimmed in file 2
    run_paired(
        "-a TTAGACATAT -m 14 -q 10",
        in1="paired.1.fastq", in2="paired.2.fastq",
        expected1="paired.m14.1.fastq", expected2="paired.m14.2.fastq",
        cores=cores
    )


def test_untrimmed_paired_output(tmpdir, run_paired):
    untrimmed1 = str(tmpdir.join("untrimmed.1.fastq"))
    untrimmed2 = str(tmpdir.join("untrimmed.2.fastq"))
    run_paired(
        ["-a", "TTAGACATAT", "--pair-filter=first",
            "--untrimmed-output", untrimmed1,
            "--untrimmed-paired-output", untrimmed2],
        in1="paired.1.fastq", in2="paired.2.fastq",
        expected1="paired-trimmed.1.fastq", expected2="paired-trimmed.2.fastq",
        cores=1
    )
    assert_files_equal(cutpath("paired-untrimmed.1.fastq"), untrimmed1)
    assert_files_equal(cutpath("paired-untrimmed.2.fastq"), untrimmed2)


def test_untrimmed_paired_output_automatic_pair_filter(tmpdir, run_paired):
    # When no R2 adapters are given, --pair-filter should be ignored for
    # --discard-untrimmed, --untrimmed-output, --untrimmed-paired-output
    # and always be "both" (with --pair-filter=any, all pairs would be
    # considered untrimmed because the R1 read is always untrimmed)
    untrimmed1 = str(tmpdir.join("untrimmed.1.fastq"))
    untrimmed2 = str(tmpdir.join("untrimmed.2.fastq"))
    run_paired(
        ["-a", "TTAGACATAT",
            "--untrimmed-output", untrimmed1,
            "--untrimmed-paired-output", untrimmed2],
        in1="paired.1.fastq", in2="paired.2.fastq",
        expected1="paired-trimmed.1.fastq", expected2="paired-trimmed.2.fastq",
        cores=1
    )
    assert_files_equal(cutpath("paired-untrimmed.1.fastq"), untrimmed1)
    assert_files_equal(cutpath("paired-untrimmed.2.fastq"), untrimmed2)


def test_explicit_format_with_paired(tmpdir, run_paired):
    # Use --format=fastq with input files whose extension is .txt
    txt1 = str(tmpdir.join("paired.1.txt"))
    txt2 = str(tmpdir.join("paired.2.txt"))
    shutil.copyfile(datapath("paired.1.fastq"), txt1)
    shutil.copyfile(datapath("paired.2.fastq"), txt2)
    run_paired(
        "--format=fastq -a TTAGACATAT -m 14 -q 10",
        in1=txt1, in2=txt2,
        expected1="paired.m14.1.fastq",
        expected2="paired.m14.2.fastq",
        cores=1
    )


def test_no_trimming_legacy():
    # make sure that this doesn"t divide by zero
    main([
        "-a", "XXXXX", "-o", "/dev/null", "-p", "/dev/null",
        datapath("paired.1.fastq"), datapath("paired.2.fastq")])


def test_no_trimming():
    # make sure that this doesn"t divide by zero
    main([
        "-a", "XXXXX", "-A", "XXXXX", "-o", "/dev/null", "-p", "/dev/null",
        datapath("paired.1.fastq"), datapath("paired.2.fastq")])


def test_missing_file(tmpdir):
    with redirect_stderr():
        with pytest.raises(SystemExit):
            main(["--paired-output", str(tmpdir.join("out.fastq")), datapath("paired.1.fastq")])


def test_first_too_short(tmpdir, cores):
    # Create a truncated file in which the last read is missing
    trunc1 = tmpdir.join("truncated.1.fastq")
    with open(datapath("paired.1.fastq")) as f:
        lines = f.readlines()
        lines = lines[:-4]
    trunc1.write("".join(lines))

    with redirect_stderr():
        with pytest.raises(SystemExit):
            main([
                "-o", "/dev/null",
                "--paired-output", str(tmpdir.join("out.fastq")),
                "--cores", str(cores),
                str(trunc1), datapath("paired.2.fastq")
            ])


def test_second_too_short(tmpdir, cores):
    # Create a truncated file in which the last read is missing
    trunc2 = tmpdir.join("truncated.2.fastq")
    with open(datapath("paired.2.fastq")) as f:
        lines = f.readlines()
        lines = lines[:-4]
    trunc2.write("".join(lines))

    with redirect_stderr():
        with pytest.raises(SystemExit):
            main([
                "-o", "/dev/null",
                "--paired-output", str(tmpdir.join("out.fastq")),
                "--cores", str(cores),
                datapath("paired.1.fastq"), str(trunc2)
            ])


def test_unmatched_read_names(tmpdir, cores):
    # Create a file in which reads 2 and 1 are swapped
    with open(datapath("paired.1.fastq")) as f:
        lines = f.readlines()
        lines = lines[0:4] + lines[8:12] + lines[4:8] + lines[12:]
    swapped = tmpdir.join("swapped.1.fastq")

    swapped.write("".join(lines))

    with pytest.raises(SystemExit):
        main([
            "-o", str(tmpdir.join("out1.fastq")),
            "--paired-output", str(tmpdir.join("out2.fastq")),
            "--cores", str(cores),
            str(swapped), datapath("paired.2.fastq")
        ])


def test_p_without_o(cores):
    """Option -p given but -o missing"""
    with pytest.raises(SystemExit):
        main("-a XX -p /dev/null".split()
            + ["--cores", str(cores)]
            + [datapath("paired.1.fastq"), datapath("paired.2.fastq")])


def test_paired_but_only_one_input_file(cores):
    """Option -p given but only one input file"""
    with pytest.raises(SystemExit):
        main("-a XX -o /dev/null -p /dev/null".split()
            + ["--cores", str(cores)]
            + [datapath("paired.1.fastq")])


def test_no_legacy_minlength(run_paired, cores):
    """Legacy mode was removed: Ensure -m is applied to second read in a pair"""
    run_paired(
        "-a XXX -m 27",
        in1="paired.1.fastq", in2="paired.2.fastq",
        expected1="paired-m27.1.fastq", expected2="paired-m27.2.fastq",
        cores=cores
    )


def test_paired_end(run_paired, cores):
    """single-pass paired-end with -m"""
    run_paired(
        "-a TTAGACATAT -A CAGTGGAGTA -m 14",
        in1="paired.1.fastq", in2="paired.2.fastq",
        expected1="paired.1.fastq", expected2="paired.2.fastq",
        cores=cores
    )


def test_paired_anchored_back_no_indels(run_paired):
    run_paired(
        "-a BACKADAPTER$ -A BACKADAPTER$ -N --no-indels",
        in1="anchored-back.fasta", in2="anchored-back.fasta",
        expected1="anchored-back.fasta", expected2="anchored-back.fasta",
        cores=1
    )


def test_paired_end_qualtrim(run_paired, cores):
    """single-pass paired-end with -q and -m"""
    run_paired(
        "-q 20 -a TTAGACATAT -A CAGTGGAGTA -m 14 -M 90",
        in1="paired.1.fastq", in2="paired.2.fastq",
        expected1="pairedq.1.fastq", expected2="pairedq.2.fastq",
        cores=cores
    )


def test_paired_end_qualtrim_swapped(run_paired, cores):
    """single-pass paired-end with -q and -m, but files swapped"""
    run_paired(
        "-q 20 -a CAGTGGAGTA -A TTAGACATAT -m 14",
        in1="paired.2.fastq", in2="paired.1.fastq",
        expected1="pairedq.2.fastq", expected2="pairedq.1.fastq",
        cores=cores
    )


def test_paired_end_cut(run_paired, cores):
    run_paired(
        "-u 3 -u -1 -U 4 -U -2",
        in1="paired.1.fastq", in2="paired.2.fastq",
        expected1="pairedu.1.fastq", expected2="pairedu.2.fastq",
        cores=cores
    )


def test_paired_end_upper_a_only(run_paired, cores):
    run_paired(
        "-A CAGTGGAGTA",
        in1="paired.1.fastq", in2="paired.2.fastq",
        expected1="paired-onlyA.1.fastq", expected2="paired-onlyA.2.fastq",
        cores=cores
    )


def test_discard_untrimmed(run_paired, cores):
    # issue #146
    # the first adapter is a sequence cut out from the first read
    run_paired(
        "-a CTCCAGCTTAGACATATC -A XXXXXXXX --discard-untrimmed",
        in1="paired.1.fastq", in2="paired.2.fastq",
        expected1="empty.fastq", expected2="empty.fastq",
        cores=cores
    )


def test_discard_trimmed(run_paired, cores):
    run_paired(
        "-A C -O 1 --discard-trimmed",  # applies everywhere
        in1="paired.1.fastq", in2="paired.2.fastq",
        expected1="empty.fastq", expected2="empty.fastq",
        cores=cores
    )


def test_interleaved_in_and_out(run_interleaved, cores):
    """Single-pass interleaved paired-end with -q and -m"""
    run_interleaved(
        "-q 20 -a TTAGACATAT -A CAGTGGAGTA -m 14 -M 90",
        inpath1="interleaved.fastq", expected1="interleaved.fastq",
        cores=cores
    )


def test_interleaved_in(run_interleaved, cores):
    """Interleaved input, two files output"""
    run_interleaved(
        "-q 20 -a TTAGACATAT -A CAGTGGAGTA -m 14 -M 90",
        inpath1="interleaved.fastq",
        expected1="pairedq.1.fastq", expected2="pairedq.2.fastq",
        cores=cores
    )


def test_interleaved_out(run_interleaved, cores):
    """Two files input, interleaved output"""
    run_interleaved(
        "-q 20 -a TTAGACATAT -A CAGTGGAGTA -m 14 -M 90",
        inpath1="paired.1.fastq", inpath2="paired.2.fastq",
        expected1="interleaved.fastq",
        cores=cores
    )


def test_interleaved_neither_nor(tmpdir):
    """Option --interleaved used, but pairs of files given for input and output"""
    p1 = str(tmpdir.join("temp-paired.1.fastq"))
    p2 = str(tmpdir.join("temp-paired.2.fastq"))
    params = "-a XX --interleaved".split()
    with redirect_stderr():
        params += ["-o", p1, "-p1", p2, "paired.1.fastq", "paired.2.fastq"]
        with pytest.raises(SystemExit):
            main(params)


def test_pair_filter_both(run_paired, cores):
    run_paired(
        "--pair-filter=both -a TTAGACATAT -A GGAGTA -m 14",
        in1="paired.1.fastq", in2="paired.2.fastq",
        expected1="paired-filterboth.1.fastq", expected2="paired-filterboth.2.fastq",
        cores=cores
    )


def test_pair_filter_first(run_paired, cores):
    run_paired(
        "--pair-filter=first -a TTAGACATAT -A GGAGTA -m 14",
        in1="paired.1.fastq", in2="paired.2.fastq",
        expected1="paired-filterfirst.1.fastq", expected2="paired-filterfirst.2.fastq",
        cores=cores
    )


def test_too_short_paired_output(run_paired, tmpdir):
    p1 = str(tmpdir.join("too-short.1.fastq"))
    p2 = str(tmpdir.join("too-short.2.fastq"))
    run_paired(
        "-a TTAGACATAT -A CAGTGGAGTA -m 14 --too-short-output "
        "{0} --too-short-paired-output {1}".format(p1, p2),
        in1="paired.1.fastq", in2="paired.2.fastq",
        expected1="paired.1.fastq", expected2="paired.2.fastq",
        cores=1
    )
    assert_files_equal(cutpath("paired-too-short.1.fastq"), p1)
    assert_files_equal(cutpath("paired-too-short.2.fastq"), p2)


def test_too_long_output(run_paired, tmpdir):
    p1 = str(tmpdir.join("too-long.1.fastq"))
    p2 = str(tmpdir.join("too-long.2.fastq"))
    run_paired(
        "-a TTAGACATAT -A CAGTGGAGTA -M 14 --too-long-output "
        "{0} --too-long-paired-output {1}".format(p1, p2),
        in1="paired.1.fastq", in2="paired.2.fastq",
        expected1="paired-too-short.1.fastq", expected2="paired-too-short.2.fastq",
        cores=1
    )
    assert_files_equal(cutpath("paired.1.fastq"), p1)
    assert_files_equal(cutpath("paired.2.fastq"), p2)


def test_too_short_output_paired_option_missing(run_paired, tmpdir):
    p1 = str(tmpdir.join("too-short.1.fastq"))
    with pytest.raises(SystemExit):
        run_paired(
            "-a TTAGACATAT -A CAGTGGAGTA -m 14 --too-short-output "
            "{0}".format(p1),
            in1="paired.1.fastq", in2="paired.2.fastq",
            expected1="paired.1.fastq", expected2="paired.2.fastq",
            cores=1
        )


def test_nextseq_paired(run_paired, cores):
    run_paired("--nextseq-trim 22", in1="nextseq.fastq", in2="nextseq.fastq",
        expected1="nextseq.fastq", expected2="nextseq.fastq",
        cores=cores)


def test_paired_demultiplex(tmpdir):
    multiout1 = str(tmpdir.join("demultiplexed.{name}.1.fastq"))
    multiout2 = str(tmpdir.join("demultiplexed.{name}.2.fastq"))
    params = [
        "-a", "first=AACATTAGACA", "-a", "second=CATTAGACATATCGG",
        "-A", "ignored=CAGTGGAGTA", "-A", "alsoignored=AATAACAGTGGAGTA",
        "-o", multiout1, "-p", multiout2,
        datapath("paired.1.fastq"), datapath("paired.2.fastq")]
    assert main(params) is None
    assert_files_equal(cutpath("demultiplexed.first.1.fastq"), multiout1.format(name="first"))
    assert_files_equal(cutpath("demultiplexed.second.1.fastq"), multiout1.format(name="second"))
    assert_files_equal(cutpath("demultiplexed.unknown.1.fastq"), multiout1.format(name="unknown"))
    assert_files_equal(cutpath("demultiplexed.first.2.fastq"), multiout2.format(name="first"))
    assert_files_equal(cutpath("demultiplexed.second.2.fastq"), multiout2.format(name="second"))
    assert_files_equal(cutpath("demultiplexed.unknown.2.fastq"), multiout2.format(name="unknown"))


@pytest.mark.parametrize("name_op,l1,l2,m", list(product(
    (("m", lambda x, y: x >= y), ("M", lambda x, y: x <= y)),
    range(1, 5),
    range(1, 5),
    [(2, 3), (2, None), (None, 3)]
)))
def test_separate_minmaxlength(tmpdir, name_op, l1, l2, m):
    """Separate minimum lengths for R1 and R2"""
    m1, m2 = m
    name, func = name_op
    inpath = str(tmpdir.join("separate_minlength.fasta"))
    expected = str(tmpdir.join("separate_minlength_expected.fasta"))
    outpath = str(tmpdir.join("out.fasta"))
    record = ">r{}:{}\n{}\n".format(l1, l2, "A" * l1)
    record += ">r{}:{}\n{}".format(l1, l2, "A" * l2)
    with open(inpath, "w") as f:
        print(record, file=f)
    with open(expected, "w") as f:
        if (m1 is None or func(l1, m1)) and (m2 is None or func(l2, m2)):
            print(record, file=f)

    assert os.path.exists(inpath)
    assert os.path.exists(expected)
    if m1 is None:
        m1 = ""
    if m2 is None:
        m2 = ""

    main(["--interleaved", "-o", outpath, "-" + name, "{}:{}".format(m1, m2), inpath])
    assert_files_equal(expected, outpath)


def test_separate_minlength_single():
    """Using separate minlengths for single-end data"""
    with pytest.raises(SystemExit):
        main(["-m", "5:7", datapath("small.fastq")])


def test_paired_end_minimal_report(run_paired, cores):
    run_paired(
        "-a TTAGACATAT -A CAGTGGAGTA -m 14 --report=minimal",
        in1="paired.1.fastq", in2="paired.2.fastq",
        expected1="paired.1.fastq", expected2="paired.2.fastq",
        cores=cores,
    )


def test_pair_adapters(run_paired, cores):
    run_paired(
        "--pair-adapters -a GTCTCCAGCT -A GACAAATAAC",
        in1="paired.1.fastq", in2="paired.2.fastq",
        expected1="pair-adapters.1.fastq", expected2="pair-adapters.2.fastq",
        cores=cores
    )


def test_pair_adapters_unequal_length(tmpdir):
    with pytest.raises(SystemExit):
        main([
            "--paired-adapters",
            "-a", "GTCTCCAGCT", "-a", "ACGTACGT",  # Two R1 adapters
            "-A", "TGCA",  # But only one R2 adapter
            "-o", str(tmpdir.join("out.1.fastq")),
            "-p", str(tmpdir.join("out.2.fastq")),
            datapath("paired.1.fastq"),
            datapath("paired.2.fastq"),
        ])


def test_pair_adapters_demultiplexing(tmpdir):
    params = "-g i1=AAAA -G i1=GGGG -g i2=CCCC -G i2=TTTT".split()
    params += ["--pair-adapters"]
    params += ["-o", str(tmpdir.join("dual-{name}.1.fastq"))]
    params += ["-p", str(tmpdir.join("dual-{name}.2.fastq"))]
    params += [datapath("dual-index.1.fastq"), datapath("dual-index.2.fastq")]
    assert main(params) is None
    for name in [
        "dual-i1.1.fastq",
        "dual-i1.2.fastq",
        "dual-i2.1.fastq",
        "dual-i2.2.fastq",
        "dual-unknown.1.fastq",
        "dual-unknown.2.fastq",
    ]:
        assert tmpdir.join(name).check()
        assert_files_equal(cutpath(name), str(tmpdir.join(name)))
