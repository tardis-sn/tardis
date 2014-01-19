#py.test configuration

import pytest


def pytest_addoption(parser):
    parser.addoption("--pure-kurucz", action="store", default=None,
        help="filename for pure kurucz")


@pytest.fixture(scope='session')
def pure_kurucz_filename(request):
    return request.config.getoption("--pure-kurucz")

def pytest_runtest_setup(item):
    if 'pure_kurucz' in item.keywords and item.config.getoption("--pure-kurucz") is None:
        pytest.skip("need --pure-kurucz option to run")
