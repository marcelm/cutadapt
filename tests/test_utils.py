from cutadapt.utils import raise_open_files_limit, Progress, DummyProgress


def test_raise_open_files_limit():
    try:
        raise_open_files_limit(1)
    except ValueError:
        pass


def test_progress():
    p = Progress()
    p.update(100)
    p.update(1000)
    p.stop(100000)


def test_dummy_progress():
    p = DummyProgress()
    p.update(100)
    p.stop(1000)
