import typing

from dnaio import Sequence
from cutadapt.adapters import (
    BackAdapter,
    AnywhereAdapter,
    BackAdapterStatistics,
    Adapter,
)
from cutadapt.modifiers import AdapterCutter, ModificationInfo


def test_statistics() -> None:
    read = Sequence("name", "AAAACCCCAAAA")
    adapters: typing.Sequence[Adapter] = [BackAdapter("CCCC", max_errors=0.1)]
    cutter = AdapterCutter(adapters, times=3)
    cutter(read, ModificationInfo(read))
    assert isinstance(cutter.adapter_statistics[adapters[0]], BackAdapterStatistics)
    lengths = cutter.adapter_statistics[adapters[0]].end.lengths
    trimmed_bp = sum(seqlen * count for (seqlen, count) in lengths.items())
    assert trimmed_bp <= len(read), trimmed_bp


def test_end_trim_with_mismatch():
    """
    Test the not-so-obvious case where an adapter of length 13 is trimmed from
    the end of a sequence with overlap 9 and there is one deletion.
    In this case the algorithm starts with 10 bases of the adapter to get
    the hit and so the match is considered good. An insertion or substitution
    at the same spot is not a match.
    """
    adapter = BackAdapter("TCGATCGATCGAT", max_errors=0.1)

    read = Sequence("foo1", "AAAAAAAAAAATCGTCGATC")
    cutter = AdapterCutter([adapter], times=1)
    trimmed_read = cutter(read, ModificationInfo(read))

    assert trimmed_read.sequence == "AAAAAAAAAAA"
    assert cutter.adapter_statistics[adapter].end.lengths == {9: 1}
    # We see 1 error at length 9 even though the number of allowed mismatches at
    # length 9 is 0.
    assert cutter.adapter_statistics[adapter].end.errors[9][1] == 1

    read = Sequence("foo2", "AAAAAAAAAAATCGAACGA")
    cutter = AdapterCutter([adapter], times=1)
    trimmed_read = cutter(read, ModificationInfo(read))

    assert trimmed_read.sequence == read.sequence
    assert cutter.adapter_statistics[adapter].end.lengths == {}


def test_anywhere_with_errors():
    adapter = AnywhereAdapter("CCGCATTTAG", max_errors=0.1)
    for seq, expected_trimmed in (
        ("AACCGGTTccgcatttagGATC", "AACCGGTT"),
        ("AACCGGTTccgcgtttagGATC", "AACCGGTT"),  # one mismatch
        ("AACCGGTTccgcatttag", "AACCGGTT"),
        ("ccgcatttagAACCGGTT", "AACCGGTT"),
        ("ccgtatttagAACCGGTT", "AACCGGTT"),  # one mismatch
        ("ccgatttagAACCGGTT", "AACCGGTT"),  # one deletion
    ):
        read = Sequence("foo", seq)
        cutter = AdapterCutter([adapter], times=1)
        trimmed_read = cutter(read, ModificationInfo(read))
        assert trimmed_read.sequence == expected_trimmed
