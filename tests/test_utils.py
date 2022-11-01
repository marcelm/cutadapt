import time
from itertools import islice

from cutadapt.utils import (
    raise_open_files_limit,
    Progress,
    DummyProgress,
)


def test_raise_open_files_limit():
    try:
        raise_open_files_limit(1)
    except ValueError:
        pass


def test_progress():
    p = Progress(every=1e-6)
    p.update(100)
    time.sleep(0.001)
    p.update(0)
    p.update(900)
    p.update(10000)
    p.close()


def test_progress_scissors():
    sc = Progress.scissors(width=10)
    for i in islice(sc, 0, 30):
        next(sc)


def test_dummy_progress():
    p = DummyProgress()
    p.update(100)
    p.update(900)
    p.close()
