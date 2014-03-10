#py.test configuration

def pytest_addoption(parser):
    parser.addoption("--atomic-dataset", dest='atomic-dataset', default=None, help="filename for atomic dataset")
