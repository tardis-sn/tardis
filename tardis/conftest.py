import yaml
import tempfile
import tardis
from tardis.atomic import AtomData
from tardis.io.config_reader import Configuration
from tardis.io.util import yaml_load_config_file
from astropy.tests.pytest_plugins import *

# For specifying error while exception handling
from socket import gaierror

try:
    import dokuwiki
except ImportError:
    dokuwiki_available = False
else:
    dokuwiki_available = True

# -------------------------------------------------------------------------
# Initialization
# -------------------------------------------------------------------------


def pytest_addoption(parser):
    parser.addoption("--remote-data", action="store_true",
                     help="run tests with online data")
    parser.addoption("--open-files", action="store_true",
                     help="fail if any test leaves files open")

    parser.addoption("--doctest-plus", action="store_true",
                     help="enable running doctests with additional "
                     "features not found in the normal doctest "
                     "plugin")

    parser.addoption("--doctest-rst", action="store_true",
                     help="enable running doctests in .rst documentation")

    parser.addini("doctest_plus", "enable running doctests with additional "
                  "features not found in the normal doctest plugin")

    parser.addini("doctest_norecursedirs",
                  "like the norecursedirs option but applies only to doctest "
                  "collection", type="args", default=())

    parser.addini("doctest_rst",
                  "Run the doctests in the rst documentation",
                  default=False)

    parser.addoption('--repeat', action='store',
                     help='Number of times to repeat each test')

    parser.addoption("--atomic-dataset", dest='atomic-dataset', default=None,
                     help="filename for atomic dataset")

    parser.addoption("--integration-tests", dest="integration-tests",
                     help="path to configuration file for integration tests")


def pytest_configure(config):
    html_file = tempfile.NamedTemporaryFile(delete=False)
    # Html test report will be generated at this filepath by pytest-html plugin
    config.option.htmlpath = html_file.name


def pytest_unconfigure(config):
    # Html report created by pytest-html plugin is read here, uploaded to
    # dokuwiki and finally deleted.
    if dokuwiki_available:
        githash = tardis.__githash__
        report_content = open(config.option.htmlpath, 'rb').read()
        report_content = report_content.replace("<!DOCTYPE html>", "")

        report_content = (
            "Test executed on commit "
            "[[https://www.github.com/tardis-sn/tardis/commit/{0}|{0}]]\n\n"
            "{1}".format(githash, report_content)
        )

        try:
            doku_conn = dokuwiki.DokuWiki(
                url=integration_tests_config["dokuwiki"]["url"],
                user=integration_tests_config["dokuwiki"]["username"],
                password=integration_tests_config["dokuwiki"]["password"])
        except gaierror, dokuwiki.DokuWikiError:
            print "Dokuwiki connection not established, report upload failed!"
        else:
            # Upload report on dokuwiki. Temporary link due to prototyping purposes.
            doku_conn.pages.append("reports:{0}".format(githash[:7]), report_content)
            print "Uploaded report on Dokuwiki."

    # Remove the local report file. Keeping the report saved on local filesystem
    # is not desired, hence deleted.
    os.unlink(config.option.htmlpath)
    print "Deleted temporary file containing html report."

# -------------------------------------------------------------------------
# hooks for influencing reporting (invoked from _pytest_terminal)
# -------------------------------------------------------------------------


def pytest_report_header(config):

    stdoutencoding = getattr(sys.stdout, 'encoding') or 'ascii'

    s = "\n"
    if six.PY2:
        args = [x.decode('utf-8') for x in config.args]
    else:
        args = config.args
    s += "Running tests in {0}.\n\n".format(" ".join(args))

    from platform import platform
    plat = platform()
    if isinstance(plat, bytes):
        plat = plat.decode(stdoutencoding, 'replace')
    s += "Platform: {0}\n\n".format(plat)
    s += "Executable: {0}\n\n".format(sys.executable)
    s += "Full Python Version: \n{0}\n\n".format(sys.version)

    s += "encodings: sys: {0}, locale: {1}, filesystem: {2}".format(
        sys.getdefaultencoding(),
        locale.getpreferredencoding(),
        sys.getfilesystemencoding())
    if sys.version_info < (3, 3, 0):
        s += ", unicode bits: {0}".format(
            int(math.log(sys.maxunicode, 2)))
    s += '\n'

    s += "byteorder: {0}\n".format(sys.byteorder)
    s += "float info: dig: {0.dig}, mant_dig: {0.dig}\n\n".format(
        sys.float_info)

    import numpy
    s += "numpy: {0}\n".format(numpy.__version__)

    try:
        import scipy
    except ImportError:
        s += "scipy: not available\n"
    else:
        s += "scipy: {0}\n".format(scipy.__version__)

    try:
        import pandas
    except ImportError:
        s += "pandas: not available\n"
    else:
        s += "pandas: {0}\n".format(pandas.__version__)

    try:
        import astropy
    except ImportError:
        s += "astropy: not available\n"
    else:
        s += "astropy: {0}\n".format(astropy.__version__)

    try:
        import yaml
    except ImportError:
        s += "yaml: not available\n"
    else:
        s += "yaml: {0}\n".format(yaml.__version__)

    try:
        import cython
    except ImportError:
        s += "cython: not available\n"
    else:
        s += "cython: {0}\n".format(cython.__version__)

    try:
        import h5py.version
    except ImportError:
        s += "h5py: not available\n"
    else:
        s += "h5py: {0}\n".format(h5py.version.version)

    try:
        import matplotlib
    except ImportError:
        s += "matplotlib: not available\n"
    else:
        s += "matplotlib: {0}\n".format(matplotlib.__version__)

    try:
        import IPython
    except ImportError:
        s += "ipython: not available\n"
    else:
        s += "ipython: {0}\n".format(IPython.__version__)

    special_opts = ["remote_data", "pep8"]
    opts = []
    for op in special_opts:
        if getattr(config.option, op, None):
            opts.append(op)
    if opts:
        s += "Using Astropy options: {0}.\n".format(" ".join(opts))

    if six.PY3 and (config.getini('doctest_rst') or config.option.doctest_rst):
        s += "Running doctests in .rst files is not supported on Python 3.x\n"

    if not six.PY3:
        s = s.encode(stdoutencoding, 'replace')

    return s


# -------------------------------------------------------------------------
# project specific fixtures
# -------------------------------------------------------------------------


@pytest.fixture(scope="session")
def atomic_data_fname():
    atomic_data_fname = pytest.config.getvalue("atomic-dataset")
    if atomic_data_fname is None:
        pytest.skip('--atomic_database was not specified')
    else:
        return os.path.expandvars(os.path.expanduser(atomic_data_fname))


@pytest.fixture
def kurucz_atomic_data(atomic_data_fname):
    atomic_data = AtomData.from_hdf5(atomic_data_fname)

    if atomic_data.md5 != '21095dd25faa1683f4c90c911a00c3f8':
        pytest.skip('Need default Kurucz atomic dataset '
                    '(md5="21095dd25faa1683f4c90c911a00c3f8"')
    else:
        return atomic_data

@pytest.fixture
def test_data_path():
    return os.path.join(tardis.__path__[0], 'tests', 'data')

@pytest.fixture
def included_he_atomic_data(test_data_path):
    atomic_db_fname = os.path.join(test_data_path, 'chianti_he_db.h5')
    return AtomData.from_hdf5(atomic_db_fname)

@pytest.fixture
def tardis_config_verysimple():
    return yaml_load_config_file(
        'tardis/io/tests/data/tardis_configv1_verysimple.yml')


@pytest.fixture(scope="session")
def integration_tests_config(request):
    integration_tests_configpath = request.config.getvalue("integration-tests")
    if integration_tests_configpath is None:
        pytest.skip('--integration-tests was not specified')
    else:
        integration_tests_configpath = os.path.expandvars(
            os.path.expanduser(integration_tests_configpath)
        )
        return yaml.load(open(integration_tests_configpath))
