# TODO
# test reading from standard input
import os
import shutil
import subprocess
import sys
import tempfile
from io import StringIO
import pytest

from cutadapt.__main__ import main
from utils import assert_files_equal, datapath, cutpath, redirect_stderr


# pytest.mark.timeout will not fail even if pytest-timeout is not installed
try:
    import pytest_timeout as _unused
except ImportError:
    raise ImportError("pytest_timeout needs to be installed")
del _unused


def test_example(run):
    run('-N -b ADAPTER', 'example.fa', 'example.fa')


def test_small(run):
    run('-a TTAGACATATCTCCGTCG', 'small.fastq', 'small.fastq')


def test_empty(run):
    """empty input"""
    run('-a TTAGACATATCTCCGTCG', 'empty.fastq', 'empty.fastq')


def test_newlines(run):
    """DOS/Windows newlines"""
    run('-e 0.12 -a TTAGACATATCTCCGTCG', 'dos.fastq', 'dos.fastq')


def test_lowercase(run):
    """lowercase adapter"""
    run('-a ttagacatatctccgtcg', 'lowercase.fastq', 'small.fastq')


def test_rest(run, tmpdir):
    """-r/--rest-file"""
    rest = str(tmpdir.join("rest.tmp"))
    run(['-b', 'ADAPTER', '-N', '-r', rest], "rest.fa", "rest.fa")
    assert_files_equal(datapath('rest.txt'), rest)


def test_restfront(run, tmpdir):
    path = str(tmpdir.join("rest.txt"))
    run(['-g', 'ADAPTER', '-N', '-r', path], "restfront.fa", "rest.fa")
    assert_files_equal(datapath('restfront.txt'), path)


def test_discard(run):
    """--discard"""
    run("-b TTAGACATATCTCCGTCG --discard", "discard.fastq", "small.fastq")


def test_discard_untrimmed(run):
    """--discard-untrimmed"""
    run('-b CAAGAT --discard-untrimmed', 'discard-untrimmed.fastq', 'small.fastq')


@pytest.mark.skip(reason='Regression since switching to dnaio')
def test_second_header_retained(run, cores):
    """test if sequence name after the "+" is retained"""
    run("--cores {} -e 0.12 -b TTAGACATATCTCCGTCG".format(cores), "plus.fastq", "plus.fastq")


@pytest.mark.skip(reason='Regression since switching to dnaio')
def test_length_tag_second_header(run, cores):
    """Ensure --length-tag= also modifies the second header line"""
    run("--cores {} -a GGCTTC --length-tag=length=".format(cores),
        'SRR2040271_1.fastq', 'SRR2040271_1.fastq')


def test_extensiontxtgz(run):
    """automatic recognition of "_sequence.txt.gz" extension"""
    run("-b TTAGACATATCTCCGTCG", "s_1_sequence.txt", "s_1_sequence.txt.gz")


def test_format(run):
    """the -f/--format parameter"""
    run("-f fastq -b TTAGACATATCTCCGTCG", "small.fastq", "small.myownextension")


def test_minimum_length(run):
    """-m/--minimum-length"""
    run("-m 5 -a TTAGACATATCTCCGTCG", "minlen.fa", "lengths.fa")


def test_too_short(run, tmpdir):
    """--too-short-output"""
    too_short_path = str(tmpdir.join('tooshort.fa'))
    run("-m 5 -a TTAGACATATCTCCGTCG --too-short-output " + too_short_path, "minlen.fa", "lengths.fa")
    assert_files_equal(datapath('tooshort.fa'), too_short_path)


def test_maximum_length(run):
    """-M/--maximum-length"""
    run("-M 5 -a TTAGACATATCTCCGTCG", "maxlen.fa", "lengths.fa")


def test_too_long(run, tmpdir):
    """--too-long-output"""
    too_long_path = str(tmpdir.join('toolong.fa'))
    run("-M 5 -a TTAGACATATCTCCGTCG --too-long-output " + too_long_path, "maxlen.fa", "lengths.fa")
    assert_files_equal(datapath('toolong.fa'), too_long_path)


def test_length_tag(run):
    """454 data; -n and --length-tag"""
    run("-n 3 -e 0.1 --length-tag length= "
        "-b TGAGACACGCAACAGGGGAAAGGCAAGGCACACAGGGGATAGG "
        "-b TCCATCTCATCCCTGCGTGTCCCATCTGTTCCCTCCCTGTCTCA", '454.fa', '454.fa')


@pytest.mark.parametrize("length", list(range(3, 11)))
def test_overlap_a(tmpdir, length):
    """-O/--overlap with -a"""
    adapter = "catatctccg"
    record = ">read\nGAGACCATTCCAATG" + adapter[:length] + '\n'
    input = tmpdir.join("overlap.fasta")
    input.write(record)
    if length < 7:
        expected = record
    else:
        expected = '>read\nGAGACCATTCCAATG\n'
    output = tmpdir.join("overlap-trimmed.fasta")
    main(["-O", "7", "-e", "0", "-a", adapter, "-o", str(output), str(input)])
    assert expected == output.read()


def test_overlap_b(run):
    """-O/--overlap with -b"""
    run("-O 10 -b TTAGACATATCTCCGTCG", "overlapb.fa", "overlapb.fa")


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
    run("-m 24 -O 10 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "polya.fasta", "polya.fasta")


def test_polya_brace_notation(run):
    """poly-A tails"""
    run("-m 24 -O 10 -a A{35}", "polya.fasta", "polya.fasta")


# the same as --action=none
def test_no_trim(run):
    run("--no-trim --discard-untrimmed -a CCCTAGTTAAAC", 'no-trim.fastq', 'small.fastq')


def test_action_none(run):
    run("--action=none --discard-untrimmed -a CCCTAGTTAAAC", 'no-trim.fastq', 'small.fastq')


# the same as --action=mask
def test_mask_adapter(run):
    """mask adapter with N (reads maintain the same length)"""
    run("-b CAAG -n 3 --mask-adapter", "anywhere_repeat.fastq", "anywhere_repeat.fastq")


def test_action_mask(run):
    """mask adapter with N (reads maintain the same length)"""
    run("-b CAAG -n 3 --action=mask", "anywhere_repeat.fastq", "anywhere_repeat.fastq")


def test_action_lowercase(run):
    run("-b CAAG -n 3 --action=lowercase", "action_lowercase.fasta", "action_lowercase.fasta")


def test_gz_multiblock(run):
    """compressed gz file with multiple blocks (created by concatenating two .gz files)"""
    run("-b TTAGACATATCTCCGTCG", "small.fastq", "multiblock.fastq.gz")


@pytest.mark.parametrize("opt", ["-y", "--suffix"])
def test_suffix(opt, run):
    """-y/--suffix parameter"""
    run([opt, ' {name}', '-e', '0', '-a', 'OnlyT=TTTTTTTT', '-a', 'OnlyG=GGGGGGGG'], "suffix.fastq", "suffix.fastq")


def test_read_wildcard(run):
    """test wildcards in reads"""
    run("--match-read-wildcards -b ACGTACGT", "wildcard.fa", "wildcard.fa")


@pytest.mark.parametrize("adapter_type,expected", [
    ("-a", "wildcard_adapter.fa"),
    ("-b", "wildcard_adapter_anywhere.fa"),
])
def test_adapter_wildcard(adapter_type, expected, run, tmpdir):
    """wildcards in adapter"""
    wildcard_path = str(tmpdir.join("wildcards.txt"))
    run("--wildcard-file {} {} ACGTNNNACGT".format(wildcard_path, adapter_type),
        expected, "wildcard_adapter.fa")
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
    run('-a TTAGACATAT -g GAGATTGCCA --no-indels', 'no_indels.fasta', 'no_indels.fasta')


def test_ellipsis_notation(run):
    run('-a ...TTAGACATAT -g GAGATTGCCA --no-indels', 'no_indels.fasta', 'no_indels.fasta')


def test_issue_46(run, tmpdir):
    """issue 46 - IndexError with --wildcard-file"""
    run("--anywhere=AACGTN --wildcard-file={}".format(
        tmpdir.join("wildcards.txt")), "issue46.fasta", "issue46.fasta")


def test_strip_suffix(run):
    run("--strip-suffix _sequence -a XXXXXXX", "stripped.fasta", "simple.fasta")


def test_info_file(run, tmpdir):
    # The true adapter sequence in the illumina.fastq.gz data set is
    # GCCTAACTTCTTAGACTGCCTTAAGGACGT (fourth base is different from the sequence shown here)
    info_path = str(tmpdir.join("info.txt"))
    run(["--info-file", info_path, "-a", "adapt=GCCGAACTTCTTAGACTGCCTTAAGGACGT"],
        "illumina.fastq", "illumina.fastq.gz")
    assert_files_equal(cutpath("illumina.info.txt"), info_path)


def test_info_file_times(run, tmpdir):
    info_path = str(tmpdir.join("info.txt"))
    run(["--info-file", info_path, "--times", "2", "-a", "adapt=GCCGAACTTCTTA",
        "-a", "adapt2=GACTGCCTTAAGGACGT"], "illumina5.fastq", "illumina5.fastq")
    assert_files_equal(cutpath('illumina5.info.txt'), info_path)


def test_info_file_fasta(run, tmpdir):
    info_path = str(tmpdir.join("info.txt"))
    # Just make sure that it runs
    run(["--info-file", info_path, "-a", "TTAGACATAT", "-g", "GAGATTGCCA", "--no-indels"],
        "no_indels.fasta", "no_indels.fasta")


def test_named_adapter(run):
    run("-a MY_ADAPTER=GCCGAACTTCTTAGACTGCCTTAAGGACGT", "illumina.fastq", "illumina.fastq.gz")


def test_adapter_with_u(run):
    run("-a GCCGAACUUCUUAGACUGCCUUAAGGACGU", "illumina.fastq", "illumina.fastq.gz")


def test_bzip2(run):
    run('-b TTAGACATATCTCCGTCG', 'small.fastq', 'small.fastq.bz2')


if sys.version_info[:2] >= (3, 3):
    def test_bzip2_multiblock(run):
        run('-b TTAGACATATCTCCGTCG', 'small.fastq', 'multiblock.fastq.bz2')


try:
    import lzma

    def test_xz(run):
        run('-b TTAGACATATCTCCGTCG', 'small.fastq', 'small.fastq.xz')
except ImportError:
    pass


def test_no_args():
    with pytest.raises(SystemExit):
        with redirect_stderr():
            main([])


def test_two_fastqs():
    with pytest.raises(SystemExit):
        with redirect_stderr():
            main([datapath('paired.1.fastq'), datapath('paired.2.fastq')])


def test_anchored_no_indels(run):
    """anchored 5' adapter, mismatches only (no indels)"""
    run('-g ^TTAGACATAT --no-indels -e 0.1', 'anchored_no_indels.fasta', 'anchored_no_indels.fasta')


def test_anchored_no_indels_wildcard_read(run):
    """anchored 5' adapter, mismatches only (no indels), but wildcards in the read count as matches"""
    run('-g ^TTAGACATAT --match-read-wildcards --no-indels -e 0.1', 'anchored_no_indels_wildcard.fasta', 'anchored_no_indels.fasta')


def test_anchored_no_indels_wildcard_adapt(run):
    """anchored 5' adapter, mismatches only (no indels), but wildcards in the adapter count as matches"""
    run('-g ^TTAGACANAT --no-indels -e 0.12', 'anchored_no_indels.fasta', 'anchored_no_indels.fasta')


def test_non_iupac_characters(run):
    with pytest.raises(SystemExit):
        with redirect_stderr():
            main(['-a', 'ZACGT', datapath('small.fastq')])


def test_unconditional_cut_front(run):
    run('-u 5', 'unconditional-front.fastq', 'small.fastq')


def test_unconditional_cut_back(run):
    run('-u -5', 'unconditional-back.fastq', 'small.fastq')


def test_unconditional_cut_both(run):
    run('-u -5 -u 5', 'unconditional-both.fastq', 'small.fastq')


def test_untrimmed_output(run, tmpdir):
    path = str(tmpdir.join("untrimmed.fastq"))
    run(["-a", "TTAGACATATCTCCGTCG", "--untrimmed-output", path], "small.trimmed.fastq", "small.fastq")
    assert_files_equal(cutpath("small.untrimmed.fastq"), path)


def test_adapter_file(run):
    run('-a file:' + datapath('adapter.fasta'), 'illumina.fastq', 'illumina.fastq.gz')


def test_adapter_file_5p_anchored(run):
    run('-N -g file:' + datapath('prefix-adapter.fasta'), 'anchored.fasta', 'anchored.fasta')


def test_adapter_file_3p_anchored(run):
    run('-N -a file:' + datapath('suffix-adapter.fasta'), 'anchored-back.fasta', 'anchored-back.fasta')


def test_adapter_file_5p_anchored_no_indels(run):
    run('-N --no-indels -g file:' + datapath('prefix-adapter.fasta'), 'anchored.fasta', 'anchored.fasta')


def test_adapter_file_3p_anchored_no_indels(run):
    run('-N --no-indels -a file:' + datapath('suffix-adapter.fasta'), 'anchored-back.fasta', 'anchored-back.fasta')


def test_adapter_file_empty_name(run):
    run('-N -a file:' + datapath('adapter-empty-name.fasta'), 'illumina.fastq', 'illumina.fastq.gz')


def test_demultiplex():
    tempdir = tempfile.mkdtemp(prefix='cutadapt-tests.')
    multiout = os.path.join(tempdir, 'tmp-demulti.{name}.fasta')
    params = ['-a', 'first=AATTTCAGGAATT', '-a', 'second=GTTCTCTAGTTCT', '-o', multiout, datapath('twoadapters.fasta')]
    assert main(params) is None
    assert_files_equal(cutpath('twoadapters.first.fasta'), multiout.format(name='first'))
    assert_files_equal(cutpath('twoadapters.second.fasta'), multiout.format(name='second'))
    assert_files_equal(cutpath('twoadapters.unknown.fasta'), multiout.format(name='unknown'))
    shutil.rmtree(tempdir)


def test_multiple_fake_anchored_adapters(run):
    run("-g ^CGTCCGAAGTAGC -g ^ATTGCCCTAG "
        "-a TTCCATGCAGCATT$ -a CCAGTCCCCCC$ "
        "-a GCCGAACTTCTTAGACTGCCTTAAGGACGT",
        "illumina.fastq",
        "illumina.fastq.gz")


def test_multiple_prefix_adapters(run):
    run("-g ^GTACGGATTGTTCAGTA -g ^TATTAAGCTCATTC", "multiprefix.fasta", "multi.fasta")


def test_multiple_prefix_adapters_noindels(run):
    run("--no-indels -g ^GTACGGATTGTTCAGTA -g ^TATTAAGCTCATTC", "multiprefix.fasta", "multi.fasta")


def test_multiple_suffix_adapters_noindels(run):
    run("--no-indels -a CGTGATTATCTTGC$ -a CCTATTAGTGGTTGAAC$", "multisuffix.fasta", "multi.fasta")


def test_max_n(run):
    run('--max-n 0', 'maxn0.fasta', 'maxn.fasta')
    run('--max-n 1', 'maxn1.fasta', 'maxn.fasta')
    run('--max-n 2', 'maxn2.fasta', 'maxn.fasta')
    run('--max-n 0.2', 'maxn0.2.fasta', 'maxn.fasta')
    run('--max-n 0.4', 'maxn0.4.fasta', 'maxn.fasta')


def test_quiet_is_quiet():
    captured_standard_output = StringIO()
    captured_standard_error = StringIO()
    old_stdout = sys.stdout
    old_stderr = sys.stderr
    try:
        sys.stdout = captured_standard_output
        sys.stderr = captured_standard_error
        main(['-o', '/dev/null', '--quiet', '-a', 'XXXX', datapath('illumina.fastq.gz')])
    finally:
        sys.stdout = old_stdout
        sys.stderr = old_stderr
    assert captured_standard_output.getvalue() == ''
    assert captured_standard_error.getvalue() == ''


def test_x_brace_notation():
    main(['-o', '/dev/null', '--quiet', '-a', 'X{5}', datapath('small.fastq')])


def test_nextseq(run):
    run('--nextseq-trim 22', 'nextseq.fastq', 'nextseq.fastq')


def test_linked_explicitly_anchored(run):
    run('-a ^AAAAAAAAAA...TTTTTTTTTT', 'linked.fasta', 'linked.fasta')


def test_linked_multiple(run):
    run('-a ^AAAAAAAAAA...TTTTTTTTTT -a ^AAAAAAAAAA...GCGCGCGCGC', 'linked.fasta', 'linked.fasta')


def test_linked_both_anchored(run):
    run('-a ^AAAAAAAAAA...TTTTT$', 'linked-anchored.fasta', 'linked.fasta')


def test_linked_5p_not_anchored(run):
    run('-g AAAAAAAAAA...TTTTTTTTTT', 'linked-not-anchored.fasta', 'linked.fasta')


def test_linked_discard_untrimmed(run):
    run('-a ^AAAAAAAAAA...TTTTTTTTTT --discard-untrimmed', 'linked-discard.fasta', 'linked.fasta')


def test_linked_discard_untrimmed_g(run):
    run('-g AAAAAAAAAA...TTTTTTTTTT --discard-untrimmed', 'linked-discard-g.fasta', 'linked.fasta')


def test_linked_anywhere():
    with pytest.raises(SystemExit):
        with redirect_stderr():
            main(['-b', 'AAA...TTT', datapath('linked.fasta')])


def test_anywhere_anchored_5p():
    with pytest.raises(SystemExit):
        with redirect_stderr():
            main(['-b', '^AAA', datapath('small.fastq')])


def test_anywhere_anchored_3p():
    with pytest.raises(SystemExit):
        with redirect_stderr():
            main(['-b', 'TTT$', datapath('small.fastq')])


def test_fasta(run):
    run('-a TTAGACATATCTCCGTCG', 'small.fasta', 'small.fastq')


def test_fasta_no_trim(run):
    run([], 'small-no-trim.fasta', 'small.fastq')


def test_length(run):
    run('--length 5', 'shortened.fastq', 'small.fastq')


def test_negative_length(run):
    run('--length -5', 'shortened-negative.fastq', 'small.fastq')


def test_run_cutadapt_process():
    subprocess.check_call(['cutadapt', '--version'])


@pytest.mark.timeout(0.5)
def test_issue_296(tmpdir):
    # Hang when using both --no-trim and --info-file together
    info_path = str(tmpdir.join('info.txt'))
    reads_path = str(tmpdir.join('reads.fasta'))
    out_path = str(tmpdir.join('out.fasta'))
    with open(reads_path, 'w') as f:
        f.write('>read\nCACAAA\n')
    main(['--info-file', info_path, '--no-trim', '-g', 'TTTCAC', '-o', out_path, reads_path])
    # Output should be unchanged because of --no-trim
    assert_files_equal(reads_path, out_path)


def test_xadapter(run):
    run('-g XTCCGAATAGA', 'xadapter.fasta', 'xadapterx.fasta')


def test_adapterx(run):
    run('-a TCCGAATAGAX', 'adapterx.fasta', 'xadapterx.fasta')


def test_discard_casava(run):
    run('--discard-casava', 'casava.fastq', 'casava.fastq')


def test_underscore(run):
    """File name ending in _fastq.gz (issue #275)"""
    run('-b TTAGACATATCTCCGTCG', 'small.fastq', 'underscore_fastq.gz')


def test_cores_autodetect(run):
    # Just make sure that it runs; functionality is not tested
    run('--cores 0 -b TTAGACATATCTCCGTCG', 'small.fastq', 'underscore_fastq.gz')


def test_write_compressed_fastq(cores, tmpdir):
    main(['--cores', str(cores), '-o', str(tmpdir.join('out.fastq.gz')), datapath('small.fastq')])


def test_minimal_report(run):
    run('-b TTAGACATATCTCCGTCG --report=minimal', 'small.fastq', 'small.fastq')


def test_paired_separate(run):
    """test separate trimming of paired-end reads"""
    run("-a TTAGACATAT", "paired-separate.1.fastq", "paired.1.fastq")
    run("-a CAGTGGAGTA", "paired-separate.2.fastq", "paired.2.fastq")


def test_run_as_module():
    """Check that "python3 -m cutadapt ..." works"""
    import subprocess
    from cutadapt import __version__
    py = subprocess.Popen([sys.executable, "-m", "cutadapt", "--version"], stdout=subprocess.PIPE)
    assert py.communicate()[0].decode().strip() == __version__


def test_standard_input_pipe(tmpdir, cores):
    """Read FASTQ from standard input"""

    import subprocess
    out_path = str(tmpdir.join("out.fastq"))
    in_path = datapath("small.fastq")
    # Use 'cat' to simulate that no file name is available for stdin
    cat = subprocess.Popen(["cat", in_path], stdout=subprocess.PIPE)
    py = subprocess.Popen([
        sys.executable, "-m", "cutadapt", "--cores", str(cores),
        "-a", "TTAGACATATCTCCGTCG", "-o", out_path, "-"],
        stdin=cat.stdout)
    _ = py.communicate()
    cat.stdout.close()
    _ = py.communicate()[0]
    assert_files_equal(cutpath("small.fastq"), out_path)


def test_standard_output(tmpdir, cores):
    """Write FASTQ to standard output (not using --output/-o option)"""

    import subprocess
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

    import subprocess
    out_path = str(tmpdir.join("out.fastq"))
    with open(out_path, "w") as out_file:
        py = subprocess.Popen([
            sys.executable, "-m", "cutadapt", "-o", "-", "--cores", str(cores),
            "-a", "TTAGACATATCTCCGTCG", datapath("small.fastq")],
            stdout=out_file)
        _ = py.communicate()
    assert_files_equal(cutpath("small.fastq"), out_path)


def test_empty_read_with_wildcard_in_adapter(run):
    run("-g CWC", "empty.fastq", "empty.fastq")
