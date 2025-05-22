import os

import pytest
from utils import assert_files_equal, datapath, cutpath
from cutadapt.cli import main
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


@pytest.fixture
def run_paired(tmp_path):
    def _run(params, in1, in2, expected1, expected2, cores):
        if type(params) is str:
            params = params.split()
        params += ["--cores", str(cores), "--buffer-size=512"]
        params += ["--json", os.fspath(tmp_path / "stats.cutadapt.json")]
        (tmp_path / "r1").mkdir()
        (tmp_path / "r2").mkdir()
        path1 = os.fspath(tmp_path / "r1" / expected1)
        path2 = os.fspath(tmp_path / "r2" / expected2)
        params += ["-o", path1, "-p", path2]
        params += [datapath(in1), datapath(in2)]
        stats = main(params)
        assert_files_equal(cutpath(expected1), path1)
        assert_files_equal(cutpath(expected2), path2)
        return stats

    return _run
