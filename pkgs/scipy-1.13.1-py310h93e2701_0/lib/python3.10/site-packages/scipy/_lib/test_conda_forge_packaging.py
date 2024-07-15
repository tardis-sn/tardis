import pytest


def test_output_separation():
    need_to_warn = False
    try:
        # check if we can import a test folder that's deleted in `scipy`
        # but present in `scipy-tests`; if this passes, skip the test
        import scipy.fft.tests
        pytest.skip("Can be ignored when `scipy-tests` is installed")
    except ModuleNotFoundError:
        # don't re-raise directly to reduce stacktrace, i.e. avoid
        # "During handling of the above exception, another exception occurred:"
        need_to_warn = True

    if need_to_warn:
        raise ImportError(
            "conda-forge builds of `scipy` do not package the tests by default "
            "to reduce the package footprint; you can ensure they're present by "
            "installing `scipy-tests` (resp. adding that to your dependencies)."
        )
