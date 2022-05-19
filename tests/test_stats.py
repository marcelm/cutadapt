from cutadapt.statistics import ReadLengthStatistics


class TestReadLengthStatistics:
    def test_empty_on_init(self):
        rls = ReadLengthStatistics()
        assert rls.written_reads() == 0
        assert rls.written_bp() == (0, 0)
        lengths = rls.written_lengths()
        assert not lengths[0] and not lengths[1]

    def test_some_reads(self):
        rls = ReadLengthStatistics()
        rls.update("THEREAD")  # length: 7
        rls.update("YETANOTHER")  # length: 10
        rls.update2("FIRST", "SECOND")  # lengths: 5, 6
        rls.update("12345")

        assert rls.written_reads() == 4
        assert rls.written_bp() == (7 + 10 + 5 + 5, 6)
        lengths = rls.written_lengths()
        assert sorted(lengths[0].items()) == [(5, 2), (7, 1), (10, 1)]
        assert sorted(lengths[1].items()) == [(6, 1)]

    def test_iadd(self):
        rls = ReadLengthStatistics()
        rls.update("THEREAD")  # length: 7
        rls.update("YETANOTHER")  # length: 10
        rls.update2("FIRST", "SECOND")  # lengths: 5, 6
        rls.update("12345")

        rls2 = ReadLengthStatistics()
        rls2.update("TESTING")  # length: 7
        rls2.update2("LEFT", "RIGHT")  # lengths: 4, 5
        rls += rls2

        assert rls.written_reads() == 6
        assert rls.written_bp() == (7 + 10 + 5 + 5 + 7 + 4, 6 + 5)
        lengths = rls.written_lengths()
        assert sorted(lengths[0].items()) == [(4, 1), (5, 2), (7, 2), (10, 1)]
        assert sorted(lengths[1].items()) == [(5, 1), (6, 1)]
