from astropy.tests.pytest_plugins import *


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

def pytest_report_header(config):

    stdoutencoding = getattr(sys.stdout, 'encoding') or 'ascii'

    s = "\n"
    if six.PY2:
        args = [x.decode('utf-8') for x in config.args]
    elif six.PY3:
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
        s += "scipy: {0}\n".format(scipy.__version__)
    except:
        s += "scipy: not available\n"

    try:
        import pandas
        s += "pandas: {0}\n".format(pandas.__version__)
    except:
        s += "pandas: not available\n"


    try:
        import astropy
    except:
        s += "astropy: not available\n"
    else:
        s += "astropy: {0}\n".format(astropy.__version__)

    try:
        import yaml
    except:
        s += "yaml: not available\n"
    else:
        s += "yaml: {0}\n".format(yaml.__version__)


    try:
        import cython
    except:
        s += "cython: not available\n"
    else:
        s += "cython: {0}\n".format(cython.__version__)



    try:
        import h5py.version
        s += "h5py: {0}\n".format(h5py.version.version)
    except:
        s += "h5py: not available\n"


    try:
        import matplotlib
        s += "matplotlib: {0}\n".format(matplotlib.__version__)
    except:
        s += "matplotlib: not available\n"

    try:
        import IPython
    except:
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