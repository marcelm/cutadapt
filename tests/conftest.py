import pytest
from utils import assert_files_equal, datapath, cutpath
from cutadapt.__main__ import main


@pytest.fixture(params=[1, 2])
def cores(request):
    return request.param


@pytest.fixture
def run(tmpdir):
    def _run(params, expected, inpath, inpath2=None):
        if type(params) is str:
            params = params.split()
        tmp_fastaq = str(tmpdir.join(expected))
        params += ['-o', tmp_fastaq]
        params += [datapath(inpath)]
        if inpath2:
            params += [datapath(inpath2)]
        assert main(params) is None
        # TODO redirect standard output
        assert_files_equal(cutpath(expected), tmp_fastaq)
        # TODO diff log files

    return _run
