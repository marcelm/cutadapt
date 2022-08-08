"""
Cutadapt doesn’t have a stable API, yet. This is an attempt to document how
one currently needs to use Cutadapt from Python to do certain things,
mostly in order to figure out where improvements need to be made.

The tests in this module do not check results, they are just here to
ensure that the code as shown can be executed.
"""
import copy
import json
import os

from utils import datapath


def test_command_line():
    # Call Cutadapt from Python, but pass parameters as a list of strings
    # the same way we would in the shell. The difference is that this is
    # not in a separate process, errors cause a CommandLineError instead
    # of SystemExit, and we get back a Statistics object.
    from cutadapt.__main__ import main

    stats = main(["-q", "10", "-o", os.devnull, datapath("small.fastq")])
    assert stats is not None
    json.dumps(stats.as_json())

    # TODO
    # - Should not set up logging
    # - Should not print anything
    # - still raises SystemExit if parser.error is called
    # - Should be cutadapt.run(...)
    # - Should the JSON stats be returned instead?


def test_pipeline_single(tmp_path):
    # The following is roughly equivalent to:
    # cutadapt -u 5 -a GATCGGAAGA -q 0,15 -m 10
    #   --discard-untrimmed --info-file=info.txt -o ... small.fastq
    import json
    from cutadapt.pipeline import SingleEndPipeline
    from cutadapt.utils import FileOpener
    from cutadapt.modifiers import UnconditionalCutter, QualityTrimmer, AdapterCutter
    from cutadapt.adapters import BackAdapter
    from cutadapt.pipeline import InputPaths, OutputFiles

    file_opener = FileOpener()
    adapter = BackAdapter(
        sequence="GATCGGAAGA",
        max_errors=1,
        min_overlap=3,
    )
    modifiers = [
        UnconditionalCutter(5),
        QualityTrimmer(cutoff_front=0, cutoff_back=15),
        AdapterCutter([adapter]),
    ]
    pipeline = SingleEndPipeline(modifiers)
    pipeline.minimum_length = (10,)
    pipeline.discard_untrimmed = True
    inpaths = InputPaths(datapath("small.fastq"))
    info_file = file_opener.xopen_or_none(tmp_path / "info.txt", "wb")
    out = file_opener.xopen(tmp_path / "out.fastq", "wb")
    outfiles = OutputFiles(info=info_file, out=out)
    stats = pipeline.run(inpaths, outfiles, cores=1)
    assert stats is not None
    json.dumps(stats.as_json())
    outfiles.close()
    # TODO
    # - info file isn’t written, what is missing?
    # - see next function for more TODOs that also apply here


def test_pipeline_paired(tmp_path):
    # cutadapt -u 5 -U 7 -a GATCGGAAGA -q 0,15 -m 10:0
    #   --discard-untrimmed --info-file=info.txt
    #   -o ... -p ...
    #   paired.1.fastq paired.2.fastq
    from cutadapt.pipeline import PairedEndPipeline
    from cutadapt.utils import FileOpener
    from cutadapt.modifiers import UnconditionalCutter, QualityTrimmer, AdapterCutter
    from cutadapt.adapters import BackAdapter
    from cutadapt.pipeline import InputPaths, OutputFiles

    trimmer = QualityTrimmer(cutoff_front=0, cutoff_back=15)
    adapter = BackAdapter(
        sequence="GATCGGAAGA",
        max_errors=1,
        min_overlap=3,
    )
    modifiers = [
        (UnconditionalCutter(5), UnconditionalCutter(7)),
        (trimmer, copy.copy(trimmer)),
        (AdapterCutter([adapter]), None),
    ]

    pipeline = PairedEndPipeline(modifiers, "any")

    pipeline.minimum_length = (10, None)
    pipeline.discard_untrimmed = True

    file_opener = FileOpener()
    inpaths = InputPaths(datapath("paired.1.fastq"), datapath("paired.2.fastq"))
    info_file = file_opener.xopen_or_none(tmp_path / "info.txt", "wb")
    out, out2 = file_opener.xopen_pair(
        tmp_path / "out.1.fastq", tmp_path / "out.2.fastq", "wb"
    )
    outfiles = OutputFiles(
        info=info_file,
        out=out,
        out2=out2,
    )
    stats = pipeline.run(inpaths, outfiles, cores=1, progress=True)
    assert stats is not None
    _ = stats.as_json()
    outfiles.close()

    # TODO
    # - could use += for adding modifiers
    # - allow using adapter specification strings?
    # - filters as attributes on Pipeline is awkward
    # - too many submodules (flatter namespace)
    # - more default arguments
    # - info file isn’t written, what is missing?
    # - use xopen directly instead of file_opener;
    #   possibly with myxopen = functools.partial(xopen, ...)


# How this could look in the future:
# def test_pipeline_single_new(tmp_path):
#     import json
#     from cutadapt.pipeline import SingleEndPipeline
#     from cutadapt.modifiers import UnconditionalCutter, QualityTrimmer, AdapterCutter
#     from cutadapt.adapters import BackAdapter
#     from cutadapt.steps import InfoFileWriter, SingleEndFilter, SingleEndSink
#     from cutadapt.filters import TooShort, DiscardUntrimmed
#     from contextlib import ExitStack
#     from xopen import xopen
#
#     with ExitStack() as stack:
#         # Input files
#         infile = stack.enter_context(xopen(datapath("small.fastq"), "rb"))
#
#         # Output files
#         info_file = stack.enter_context(xopen(tmp_path / "info.txt", "wb"))
#         out = stack.enter_context(xopen(tmp_path / "out.fastq", "wb"))
#         too_short = stack.enter_context(xopen(tmp_path / "tooshort.fastq", "wb"))
#         pipeline = SingleEndPipeline(
#             modifiers=[
#                 UnconditionalCutter(5),
#                 QualityTrimmer(cutoff_front=0, cutoff_back=15),
#                 AdapterCutter([BackAdapter(sequence="GATCGGAAGA", max_errors=1, min_overlap=3)]),
#             ],
#             filters=[
#                 InfoFileWriter(info_file),
#                 SingleEndFilter(TooShort(minimum_length=10), too_short),
#                 SingleEndFilter(DiscardUntrimmed()),
#             ],
#             sink=SingleEndSink(out),
#         )
#         stats = pipeline.run(infile, workers=1, progress=True)
#     json.dumps(stats.as_json())
