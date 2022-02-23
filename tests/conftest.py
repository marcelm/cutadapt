import os

import pytest
from utils import assert_files_equal, datapath, cutpath
from cutadapt.__main__ import main
from cutadapt.report import Statistics


@pytest.fixture(params=[1, 2])
def cores(request):
    return request.param


@pytest.fixture
def run(tmp_path):
    def _run(params, expected, inpath) -> Statistics:
        if type(params) is str:
            params = params.split()
        params += ["--json", os.fspath(tmp_path / "stats.cutadapt.json")]
        tmp_fastaq = tmp_path / expected
        params += ["-o", tmp_fastaq]
        params += [datapath(inpath)]
        stats = main([str(p) for p in params])
        # TODO redirect standard output
        assert_files_equal(cutpath(expected), tmp_fastaq)
        return stats

    return _run
