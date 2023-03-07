import subprocess
import sys
import os
from io import StringIO, BytesIO

import dnaio
import pytest

from cutadapt.__main__ import main
from utils import assert_files_equal, datapath, cutpath

# pytest.mark.timeout will not fail even if pytest-timeout is not installed
try:
    import pytest_timeout as _unused
except ImportError:  # pragma: no cover
    raise ImportError("pytest_timeout needs to be installed")
del _unused


def test_does_not_close_stdout():
    main([datapath("small.fastq")])
    assert not sys.stdout.closed


def test_help():
    with pytest.raises(SystemExit) as e:
        main(["--help"])
    assert e.value.args[0] == 0


def test_unknown_file_format(tmp_path):
    path = tmp_path / "unknown_format.txt"
    path.write_text("raw text")
    with pytest.raises(SystemExit):
        main([str(path)])


def test_cores_negative():
    with pytest.raises(SystemExit) as e:
        main(["--cores=-1", datapath("simple.fasta")])
    assert e.value.args[0] == 2
    # "cannot be negative"


def test_quiet_and_report():
    with pytest.raises(SystemExit) as e:
        main(["--quiet", "--report=minimal", datapath("simple.fasta")])
    assert e.value.args[0] == 2
    # "Options --quiet and --report cannot be used at the same time"


@pytest.mark.parametrize(
    "args",
    [
        ("--discard-trimmed", "--discard-untrimmed"),
        ("--discard-trimmed", "--untrimmed-output", os.devnull),
        ("--discard-untrimmed", "--untrimmed-output", os.devnull),
    ],
)
def test_only_one_of_discard_trimmed_discard_untrimmed_untrimmed_output(args):
    with pytest.raises(SystemExit) as e:
        main(["-o", os.devnull, *args, datapath("small.fastq")])
    assert e.value.args[0] == 2


def test_debug():
    main(["--debug", "--", datapath("small.fastq")])


def test_debug_trace():
    main(["--debug", "--debug", "-a", "ACGT", datapath("small.fastq")])


def test_example(run):
    run("-N -b ADAPTER", "example.fa", "example.fa")


def test_compressed_fasta(run):
    run("", "simple.fasta", "simple.fasta.gz")


def test_small(run):
    run("-a TTAGACATATCTCCGTCG", "small.fastq", "small.fastq")


def test_empty_fastq(run, cores):
    run("--cores {} -a TTAGACATATCTCCGTCG".format(cores), "empty.fastq", "empty.fastq")


def test_empty_fasta_input(run, cores):
    run(["--cores", str(cores)], "empty.fasta", "empty.fasta")


def test_no_read_only_comment_fasta_input(run, cores):
    run(["--cores", str(cores)], "empty.fasta", "onlycomment.fasta")


def test_newlines(run):
    """DOS/Windows newlines"""
    run("-e 0.12 -a TTAGACATATCTCCGTCG", "dos.fastq", "dos.fastq")


def test_lowercase(run):
    """lowercase adapter"""
    run("-a ttagacatatctccgtcg", "lowercase.fastq", "small.fastq")


def test_rest(run, tmp_path, cores):
    """-r/--rest-file"""
    rest = tmp_path / "rest.tmp"
    run(
        ["--cores", str(cores), "-b", "ADAPTER", "-N", "-r", rest], "rest.fa", "rest.fa"
    )
    assert_files_equal(datapath("rest.txt"), rest)


def test_restfront(run, tmp_path):
    path = tmp_path / "rest.txt"
    run(["-g", "ADAPTER", "-N", "-r", path], "restfront.fa", "rest.fa")
    assert_files_equal(datapath("restfront.txt"), path)


def test_discard(run):
    """--discard"""
    run("-b TTAGACATATCTCCGTCG --discard", "discard.fastq", "small.fastq")


def test_discard_untrimmed(run):
    """--discard-untrimmed"""
    run("-b CAAGAT --discard-untrimmed", "discard-untrimmed.fastq", "small.fastq")


def test_extensiontxtgz(run):
    """automatic recognition of "_sequence.txt.gz" extension"""
    run("-b TTAGACATATCTCCGTCG", "s_1_sequence.txt", "s_1_sequence.txt.gz")


def test_minimum_length(run):
    """-m/--minimum-length"""
    stats = run("-m 5 -a TTAGACATATCTCCGTCG", "minlen.fa", "lengths.fa")
    assert stats.written_bp[0] == 45
    assert stats.written == 6


def test_too_short(run, tmp_path, cores):
    too_short_path = tmp_path / "tooshort.fa"
    stats = run(
        [
            "--cores",
            str(cores),
            "-m",
            "5",
            "-a",
            "TTAGACATATCTCCGTCG",
            "--too-short-output",
            too_short_path,
        ],
        "minlen.fa",
        "lengths.fa",
    )
    assert_files_equal(datapath("tooshort.fa"), too_short_path)
    assert stats.filtered["too_short"] == 5


@pytest.mark.parametrize("redirect", (False, True))
def test_too_short_statistics(redirect):
    args = [
        "-a",
        "TTAGACATATCTCCGTCG",
        "-m",
        "24",
        "-o",
        os.devnull,
        datapath("small.fastq"),
    ]
    if redirect:
        args[:0] = ["--too-short-output", os.devnull]
    stats = main(args)
    assert stats.with_adapters[0] == 2
    assert stats.written == 2
    assert stats.written_bp[0] == 58
    assert stats.filtered["too_short"] == 1


def test_maximum_length(run):
    """-M/--maximum-length"""
    run("-M 5 -a TTAGACATATCTCCGTCG", "maxlen.fa", "lengths.fa")


def test_too_long(run, tmp_path, cores):
    """--too-long-output"""
    too_long_path = tmp_path / "toolong.fa"
    stats = run(
        [
            "--cores",
            str(cores),
            "-M",
            "5",
            "-a",
            "TTAGACATATCTCCGTCG",
            "--too-long-output",
            too_long_path,
        ],
        "maxlen.fa",
        "lengths.fa",
    )
    assert_files_equal(datapath("toolong.fa"), too_long_path)
    assert stats.filtered["too_long"] == 5


def test_length_tag(run):
    """454 data; -n and --length-tag"""
    run(
        "-n 3 -e 0.1 --length-tag length= "
        "-b TGAGACACGCAACAGGGGAAAGGCAAGGCACACAGGGGATAGG "
        "-b TCCATCTCATCCCTGCGTGTCCCATCTGTTCCCTCCCTGTCTCA",
        "454.fa",
        "454.fa",
    )


@pytest.mark.parametrize("length", list(range(3, 11)))
def test_overlap_a(tmp_path, length):
    """-O/--overlap with -a"""
    adapter = "catatctccg"
    record = ">read\nGAGACCATTCCAATG" + adapter[:length] + "\n"
    input = tmp_path / "overlap.fasta"
    input.write_text(record)
    if length < 7:
        expected = record
    else:
        expected = ">read\nGAGACCATTCCAATG\n"
    output = tmp_path / "overlap-trimmed.fasta"
    main(["-O", "7", "-e", "0", "-a", adapter, "-o", str(output), str(input)])
    assert expected == output.read_text()


def test_overlap_b(run):
    """-O/--overlap with -b"""
    run("-O 10 -b TTAGACATATCTCCGTCG", "overlapb.fa", "overlapb.fa")


def test_trim_n(run):
    run("--trim-n", "trim-n.fasta", "trim-n.fasta")


def test_qualtrim(run):
    """-q with low qualities"""
    run("-q 10 -a XXXXXX", "lowqual.fastq", "lowqual.fastq")


def test_qualbase(run):
    """-q with low qualities, using ascii(quality+64) encoding"""
    run("-q 10 --quality-base 64 -a XXXXXX", "illumina64.fastq", "illumina64.fastq")


def test_quality_trim_only(run):
    """only trim qualities, do not remove adapters"""
    run("-q 10 --quality-base 64", "illumina64.fastq", "illumina64.fastq")


def test_twoadapters(run):
    """two adapters"""
    run("-a AATTTCAGGAATT -a GTTCTCTAGTTCT", "twoadapters.fasta", "twoadapters.fasta")


def test_polya(run):
    """poly-A tails"""
    run(
        "-m 24 -O 10 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
        "polya.fasta",
        "polya.fasta",
    )


def test_polya_brace_notation(run):
    """poly-A tails"""
    run("-m 24 -O 10 -a A{35}", "polya.fasta", "polya.fasta")


# the same as --action=none
def test_no_trim(run):
    run("--no-trim --discard-untrimmed -a CCCTAGTTAAAC", "no-trim.fastq", "small.fastq")


def test_action_none(run):
    run(
        "--action=none --discard-untrimmed -a CCCTAGTTAAAC",
        "no-trim.fastq",
        "small.fastq",
    )


# the same as --action=mask
def test_mask_adapter(run):
    """mask adapter with N (reads maintain the same length)"""
    run("-b CAAG -n 3 --mask-adapter", "anywhere_repeat.fastq", "anywhere_repeat.fastq")


def test_action_mask(run):
    """mask adapter with N (reads maintain the same length)"""
    run("-b CAAG -n 3 --action=mask", "anywhere_repeat.fastq", "anywhere_repeat.fastq")


def test_action_lowercase(run):
    run(
        "-b CAAG -n 3 --action=lowercase",
        "action_lowercase.fasta",
        "action_lowercase.fasta",
    )


def test_action_retain(run):
    run(
        "-g GGTTAACC -a CAAG --action=retain",
        "action_retain.fasta",
        "action_retain.fasta",
    )


def test_action_retain_times():
    with pytest.raises(SystemExit):
        main(["-a", "ACGT", "--times=2", "--action=retain", datapath("small.fastq")])


def test_gz_multiblock(run):
    """compressed gz file with multiple blocks (created by concatenating two .gz files)"""
    run("-b TTAGACATATCTCCGTCG", "small.fastq", "multiblock.fastq.gz")


def test_read_wildcard(run):
    """test wildcards in reads"""
    run("--match-read-wildcards -b ACGTACGT", "wildcard.fa", "wildcard.fa")


@pytest.mark.parametrize(
    "adapter_type,expected",
    [
        ("-a", "wildcard_adapter.fa"),
        ("-b", "wildcard_adapter_anywhere.fa"),
    ],
)
def test_adapter_wildcard(adapter_type, expected, run, tmp_path, cores):
    """wildcards in adapter"""
    wildcard_path = tmp_path / "wildcards.txt"
    run(
        [
            "--cores",
            str(cores),
            "--wildcard-file",
            wildcard_path,
            adapter_type,
            "ACGTNNNACGT",
        ],
        expected,
        "wildcard_adapter.fa",
    )
    with open(wildcard_path) as wct:
        lines = wct.readlines()
    lines = [line.strip() for line in lines]
    assert lines == ["AAA 1", "GGG 2", "CCC 3b", "TTT 4b"]


def test_wildcard_N(run):
    """test 'N' wildcard matching with no allowed errors"""
    run("-e 0 -a GGGGGGG --match-read-wildcards", "wildcardN.fa", "wildcardN.fa")


def test_illumina_adapter_wildcard(run):
    run("-a VCCGAMCYUCKHRKDCUBBCNUWNSGHCGU", "illumina.fastq", "illumina.fastq.gz")


def test_adapter_front(run):
    """test adapter in front"""
    run("--front ADAPTER -N", "examplefront.fa", "example.fa")


def test_literal_N(run):
    """test matching literal 'N's"""
    run("-N -e 0.2 -a NNNNNNNNNNNNNN", "trimN3.fasta", "trimN3.fasta")


def test_literal_N2(run):
    run("-N -O 1 -g NNNNNNNNNNNNNN", "trimN5.fasta", "trimN5.fasta")


def test_literal_N_brace_notation(run):
    """test matching literal 'N's"""
    run("-N -e 0.2 -a N{14}", "trimN3.fasta", "trimN3.fasta")


def test_literal_N2_brace_notation(run):
    run("-N -O 1 -g N{14}", "trimN5.fasta", "trimN5.fasta")


def test_anchored_front(run):
    run("-g ^FRONTADAPT -N", "anchored.fasta", "anchored.fasta")


def test_anchored_front_ellipsis_notation(run):
    run("-a ^FRONTADAPT... -N", "anchored.fasta", "anchored.fasta")


def test_anchored_back(run):
    run("-a BACKADAPTER$ -N", "anchored-back.fasta", "anchored-back.fasta")


def test_anchored_back_ellipsis_notation(run):
    run("-a ...BACKADAPTER$ -N", "anchored-back.fasta", "anchored-back.fasta")


def test_anchored_back_no_indels(run):
    run("-a BACKADAPTER$ -N --no-indels", "anchored-back.fasta", "anchored-back.fasta")


def test_no_indels(run):
    run("-a TTAGACATAT -g GAGATTGCCA --no-indels", "no_indels.fasta", "no_indels.fasta")


def test_ellipsis_notation(run):
    run(
        "-a ...TTAGACATAT -g GAGATTGCCA --no-indels",
        "no_indels.fasta",
        "no_indels.fasta",
    )


def test_issue_46(run, tmp_path):
    """issue 46 - IndexError with --wildcard-file"""
    run(
        "--anywhere=AACGTN --wildcard-file={}".format(tmp_path / "wildcards.txt"),
        "issue46.fasta",
        "issue46.fasta",
    )


def test_strip_suffix(run):
    run("--strip-suffix _sequence -a XXXXXXX", "stripped.fasta", "simple.fasta")


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


def test_named_adapter(run):
    run(
        "-a MY_ADAPTER=GCCGAACTTCTTAGACTGCCTTAAGGACGT",
        "illumina.fastq",
        "illumina.fastq.gz",
    )


def test_adapter_with_u(run):
    run("-a GCCGAACUUCUUAGACUGCCUUAAGGACGU", "illumina.fastq", "illumina.fastq.gz")


def test_bzip2_input(run, cores):
    run(
        ["--cores", str(cores), "-a", "TTAGACATATCTCCGTCG"],
        "small.fastq",
        "small.fastq.bz2",
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


def test_no_args():
    with pytest.raises(SystemExit):
        main([])


def test_two_fastqs():
    with pytest.raises(SystemExit):
        main([datapath("paired.1.fastq"), datapath("paired.2.fastq")])


def test_anchored_no_indels(run):
    """anchored 5' adapter, mismatches only (no indels)"""
    run(
        "-g ^TTAGACATAT --no-indels -e 0.1",
        "anchored_no_indels.fasta",
        "anchored_no_indels.fasta",
    )


def test_anchored_no_indels_wildcard_read(run):
    """anchored 5' adapter, mismatches only (no indels), but wildcards in the read count as matches"""
    run(
        "-g ^TTAGACATAT --match-read-wildcards --no-indels -e 0.1",
        "anchored_no_indels_wildcard.fasta",
        "anchored_no_indels.fasta",
    )


def test_anchored_no_indels_wildcard_adapt(run):
    """anchored 5' adapter, mismatches only (no indels), but wildcards in the adapter count as matches"""
    run(
        "-g ^TTAGACANAT --no-indels -e 0.12",
        "anchored_no_indels.fasta",
        "anchored_no_indels.fasta",
    )


def test_non_iupac_characters(run):
    with pytest.raises(SystemExit):
        main(["-a", "ZACGT", datapath("small.fastq")])


def test_unconditional_cut_front(run):
    run("-u 5", "unconditional-front.fastq", "small.fastq")


def test_unconditional_cut_back(run):
    run("-u -5", "unconditional-back.fastq", "small.fastq")


def test_unconditional_cut_both(run):
    run("-u -5 -u 5", "unconditional-both.fastq", "small.fastq")


def test_unconditional_cut_too_many_commas():
    with pytest.raises(SystemExit):
        main(["-u", "5,7,8", datapath("small.fastq")])


def test_unconditional_cut_invalid_number():
    with pytest.raises(SystemExit):
        main(["-u", "a,b", datapath("small.fastq")])


def test_untrimmed_output(run, cores, tmp_path):
    path = tmp_path / "untrimmed.fastq"
    stats = run(
        ["--cores", str(cores), "-a", "TTAGACATATCTCCGTCG", "--untrimmed-output", path],
        "small.trimmed.fastq",
        "small.fastq",
    )
    assert_files_equal(cutpath("small.untrimmed.fastq"), path)
    assert stats.with_adapters[0] == 2
    assert stats.written == 2
    assert stats.written_bp[0] == 46


def test_adapter_file(run):
    run("-a file:" + datapath("adapter.fasta"), "illumina.fastq", "illumina.fastq.gz")


def test_adapter_file_5p_anchored(run):
    run(
        "-N -g file:" + datapath("prefix-adapter.fasta"),
        "anchored.fasta",
        "anchored.fasta",
    )


def test_adapter_file_3p_anchored(run):
    run(
        "-N -a file:" + datapath("suffix-adapter.fasta"),
        "anchored-back.fasta",
        "anchored-back.fasta",
    )


def test_adapter_file_5p_anchored_no_indels(run):
    run(
        "-N --no-indels -g file:" + datapath("prefix-adapter.fasta"),
        "anchored.fasta",
        "anchored.fasta",
    )


def test_adapter_file_3p_anchored_no_indels(run):
    run(
        "-N --no-indels -a file:" + datapath("suffix-adapter.fasta"),
        "anchored-back.fasta",
        "anchored-back.fasta",
    )


def test_adapter_file_empty_name(run):
    run(
        "-N -a file:" + datapath("adapter-empty-name.fasta"),
        "illumina.fastq",
        "illumina.fastq.gz",
    )


@pytest.mark.parametrize("ext", ["", ".gz"])
def test_demultiplex(cores, tmp_path, ext):
    multiout = str(tmp_path / "tmp-demulti.{name}.fasta") + ext
    params = [
        "--cores",
        str(cores),
        "-a",
        "first=AATTTCAGGAATT",
        "-a",
        "second=GTTCTCTAGTTCT",
        "-o",
        multiout,
        datapath("twoadapters.fasta"),
    ]
    main(params)
    for name in ("first", "second", "unknown"):
        actual = multiout.format(name=name)
        if ext == ".gz":
            subprocess.run(["gzip", "-d", actual], check=True)
            actual = actual[:-3]
        expected = cutpath("twoadapters.{name}.fasta".format(name=name))
        assert_files_equal(expected, actual)


def test_multiple_fake_anchored_adapters(run):
    run(
        "-g ^CGTCCGAAGTAGC -g ^ATTGCCCTAG "
        "-a TTCCATGCAGCATT$ -a CCAGTCCCCCC$ "
        "-a GCCGAACTTCTTAGACTGCCTTAAGGACGT",
        "illumina.fastq",
        "illumina.fastq.gz",
    )


def test_multiple_prefix_adapters(run):
    run("-g ^GTACGGATTGTTCAGTA -g ^TATTAAGCTCATTC", "multiprefix.fasta", "multi.fasta")


def test_multiple_prefix_adapters_noindels(run):
    run(
        "--no-indels -g ^GTACGGATTGTTCAGTA -g ^TATTAAGCTCATTC",
        "multiprefix.fasta",
        "multi.fasta",
    )


def test_multiple_suffix_adapters_noindels(run):
    run(
        "--no-indels -a CGTGATTATCTTGC$ -a CCTATTAGTGGTTGAAC$",
        "multisuffix.fasta",
        "multi.fasta",
    )


def test_max_n(run):
    assert run("--max-n 0", "maxn0.fasta", "maxn.fasta").filtered["too_many_n"] == 4
    assert run("--max-n 1", "maxn1.fasta", "maxn.fasta").filtered["too_many_n"] == 2
    assert run("--max-n 2", "maxn2.fasta", "maxn.fasta").filtered["too_many_n"] == 1
    assert run("--max-n 0.2", "maxn0.2.fasta", "maxn.fasta").filtered["too_many_n"] == 3
    assert run("--max-n 0.4", "maxn0.4.fasta", "maxn.fasta").filtered["too_many_n"] == 2


def test_quiet_is_quiet():
    captured_standard_output = StringIO()
    captured_standard_error = StringIO()
    setattr(captured_standard_output, "buffer", BytesIO())
    setattr(captured_standard_error, "buffer", BytesIO())
    old_stdout = sys.stdout
    old_stderr = sys.stderr
    try:
        sys.stdout = captured_standard_output
        sys.stderr = captured_standard_error
        main(["-o", os.devnull, "--quiet", datapath("small.fastq")])
    finally:
        sys.stdout = old_stdout
        sys.stderr = old_stderr
    assert captured_standard_output.getvalue() == ""
    assert captured_standard_error.getvalue() == ""
    assert getattr(captured_standard_output, "buffer").getvalue() == b""
    assert getattr(captured_standard_output, "buffer").getvalue() == b""


def test_x_brace_notation():
    main(["-o", os.devnull, "--quiet", "-a", "X{5}", datapath("small.fastq")])


def test_nextseq(run):
    run("--nextseq-trim 22", "nextseq.fastq", "nextseq.fastq")


def test_linked_explicitly_anchored(run):
    run("-a ^AAAAAAAAAA...TTTTTTTTTT", "linked.fasta", "linked.fasta")


def test_linked_multiple(run):
    run(
        "-a ^AAAAAAAAAA...TTTTTTTTTT -a ^AAAAAAAAAA...GCGCGCGCGC",
        "linked.fasta",
        "linked.fasta",
    )


def test_linked_both_anchored(run):
    run("-a ^AAAAAAAAAA...TTTTT$", "linked-anchored.fasta", "linked.fasta")


def test_linked_5p_not_anchored(run):
    run("-g AAAAAAAAAA...TTTTTTTTTT", "linked-not-anchored.fasta", "linked.fasta")


def test_linked_discard_untrimmed(run):
    run(
        "-a ^AAAAAAAAAA...TTTTTTTTTT --discard-untrimmed",
        "linked-discard.fasta",
        "linked.fasta",
    )


def test_linked_discard_untrimmed_g(run):
    run(
        "-g AAAAAAAAAA...TTTTTTTTTT --discard-untrimmed",
        "linked-discard-g.fasta",
        "linked.fasta",
    )


def test_linked_lowercase(run):
    run(
        "-a ^AACCGGTTTT...GGGGGGG$ -a ^AAAA...TTTT$ --times=2 --action=lowercase",
        "linked-lowercase.fasta",
        "linked.fasta",
    )


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


def test_linked_anywhere():
    with pytest.raises(SystemExit):
        main(["-b", "AAA...TTT", datapath("linked.fasta")])


def test_anywhere_anchored_5p():
    with pytest.raises(SystemExit):
        main(["-b", "^AAA", datapath("small.fastq")])


def test_anywhere_anchored_3p():
    with pytest.raises(SystemExit):
        main(["-b", "TTT$", datapath("small.fastq")])


def test_fasta(run):
    run("-a TTAGACATATCTCCGTCG", "small.fasta", "small.fastq")


def test_fasta_no_trim(run):
    run([], "small-no-trim.fasta", "small.fastq")


def test_length(run):
    run("--length 5", "shortened.fastq", "small.fastq")


def test_negative_length(run):
    run("--length -5", "shortened-negative.fastq", "small.fastq")


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


def test_xadapter(run):
    run("-g XTCCGAATAGA", "xadapter.fasta", "xadapterx.fasta")


def test_adapterx(run):
    run("-a TCCGAATAGAX", "adapterx.fasta", "xadapterx.fasta")


def test_not_rightmost(tmp_path):
    path = tmp_path / "reads.fasta"
    path.write_text(">r\nGGCTGAATTGGACTGAATTGGGT\n")
    trimmed = tmp_path / "trimmed.fasta"
    main(["-g", "CTGAATT", "-o", str(trimmed), str(path)])
    assert trimmed.read_text() == ">r\nGGACTGAATTGGGT\n"


def test_rightmost(tmp_path):
    path = tmp_path / "reads.fasta"
    path.write_text(">r\nGGCTGAATTGGACTGAATTGGGT\n")
    trimmed = tmp_path / "trimmed.fasta"
    main(["-g", "CTGAATT;rightmost", "-o", str(trimmed), str(path)])
    assert trimmed.read_text() == ">r\nGGGT\n"


def test_discard_casava(run):
    stats = run("--discard-casava", "casava.fastq", "casava.fastq")
    assert stats.filtered["casava_filtered"] == 1


def test_underscore(run):
    """File name ending in _fastq.gz (issue #275)"""
    run("-b TTAGACATATCTCCGTCG", "small.fastq", "underscore_fastq.gz")


def test_cores_autodetect(run):
    # Just make sure that it runs; functionality is not tested
    run("--cores 0 -b TTAGACATATCTCCGTCG", "small.fastq", "underscore_fastq.gz")


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


def test_minimal_report(run):
    run("-b TTAGACATATCTCCGTCG --report=minimal", "small.fastq", "small.fastq")


def test_paired_separate(run):
    """test separate trimming of paired-end reads"""
    run("-a TTAGACATAT", "paired-separate.1.fastq", "paired.1.fastq")
    run("-a CAGTGGAGTA", "paired-separate.2.fastq", "paired.2.fastq")


def test_empty_read_with_wildcard_in_adapter(run):
    run("-g CWC", "empty.fastq", "empty.fastq")


def test_print_progress_to_tty(tmp_path, mocker):
    mocker.patch("cutadapt.utils.sys.stderr").isatty.return_value = True
    main(["-o", str(tmp_path / "out.fastq"), datapath("small.fastq")])


def test_adapter_order(run):
    run("-g ^AAACC -a CCGGG", "adapterorder-ga.fasta", "adapterorder.fasta")
    run("-a CCGGG -g ^AAACC", "adapterorder-ag.fasta", "adapterorder.fasta")


def test_reverse_complement_no_rc_suffix(run, tmp_path):
    out_path = tmp_path / "out.fastq"
    main(
        [
            "-o",
            str(out_path),
            "--revcomp",
            "--no-index",
            "--rename",
            "{header}",
            "-g",
            "^TTATTTGTCT",
            "-g",
            "^TCCGCACTGG",
            datapath("revcomp.1.fastq"),
        ]
    )
    with dnaio.open(out_path) as f:
        reads = list(f)
    assert len(reads) == 6
    assert reads[1].name == "read2/1"
    assert reads[1].sequence == "ACCATCCGATATGTCTAATGTGGCCTGTTG"


def test_reverse_complement_normalized(run):
    stats = run(
        "--revcomp --no-index -g ^TTATTTGTCT -g ^TCCGCACTGG",
        "revcomp-single-normalize.fastq",
        "revcomp.1.fastq",
    )
    assert stats.n == 6
    assert stats.reverse_complemented == 2


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


def test_max_expected_errors(run, cores):
    stats = run("--max-ee=0.9", "maxee.fastq", "maxee.fastq")
    assert stats.filtered["too_many_expected_errors"] == 2


def test_max_expected_errors_fasta(tmp_path):
    path = tmp_path / "input.fasta"
    path.write_text(">read\nACGTACGT\n")
    main(["--max-ee=0.001", "-o", os.devnull, str(path)])


def test_warn_if_en_dashes_used():
    with pytest.raises(SystemExit):
        main(["–q", "25", "-o", os.devnull, "in.fastq"])


@pytest.mark.parametrize("opt", ["-y", "--suffix"])
def test_suffix(opt, run):
    """-y/--suffix parameter"""
    run(
        [opt, " {name}", "-e", "0", "-a", "OnlyT=TTTTTTTT", "-a", "OnlyG=GGGGGGGG"],
        "suffix.fastq",
        "suffix.fastq",
    )


@pytest.mark.parametrize("opt", ["--prefix", "--suffix"])
def test_rename_cannot_be_combined_with_other_renaming_options(opt):
    with pytest.raises(SystemExit):
        main(
            [
                opt,
                "something",
                "--rename='{id} {comment} extrainfo'",
                "-o",
                os.devnull,
                datapath("empty.fastq"),
            ]
        )


@pytest.mark.skipif(sys.platform == "win32", reason="Disabled on Windows")
def test_duplicate_output_paths(tmp_path):
    path = str(tmp_path / "discard.fastq")
    with pytest.raises(SystemExit):
        main(
            [
                "--untrimmed-output",
                path,
                "--too-long-output",
                path,
                "-o",
                os.devnull,
                datapath("empty.fastq"),
            ]
        )
    # "specified more than once as an output file"


def test_rename(run):
    run(
        [
            "--rename={id}_{cut_suffix} {header} {adapter_name}",
            "--cut=-4",
            "-a",
            "OnlyT=TTTTTT",
            "-a",
            "OnlyG=GGGGGG",
        ],
        "rename.fastq",
        "suffix.fastq",
    )


def test_terminates_correctly_on_error_in_subprocess(tmp_path):
    params = [
        "-j",
        "2",
        "-o",
        str(tmp_path / "out.fastq.gz"),
        datapath("format-error.fastq"),
    ]
    with pytest.raises(SystemExit):
        main(params)


def test_json_report_with_demultiplexing_and_discard_untrimmed(tmp_path):
    stats = main(
        [
            "--json",
            str(tmp_path / "demux.cutadapt.json"),
            "--discard-untrimmed",
            "-a",
            "name=ACGT",
            "-o",
            str(tmp_path / "{name}.fastq"),
            datapath("illumina.fastq.gz"),
        ]
    )
    assert stats.n == 100
    assert stats.written == 64
