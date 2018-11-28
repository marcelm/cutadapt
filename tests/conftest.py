import pytest


@pytest.fixture(params=[1, 2])
def cores(request):
    return request.param
