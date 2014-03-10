#py.test configuration

def pytest_addoption(parser):
    parser.addoption("--atomic-data-set", default=None, help="filename for atomic dataset")

