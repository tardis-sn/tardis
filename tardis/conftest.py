#py.test configuration

from astropy.tests.pytest_plugins import *

def pytest_addoption(parser):
    parser.addoption("--atomic-dataset", dest='atomic-dataset', default=None, help="filename for atomic dataset")
