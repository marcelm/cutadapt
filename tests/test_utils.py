from cutadapt.utils import raise_open_files_limit, Progress

def test_raise_open_files_limit():
    raise_open_files_limit(1)


def test_progress():
    p = Progress()
    p.update(100)
    p.update(1000)
    p.stop(100000)
