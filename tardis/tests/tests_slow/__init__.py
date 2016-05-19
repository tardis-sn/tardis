import pytest

runslow = pytest.mark.skipif(
    not pytest.config.getoption("--slow") or
    not pytest.config.getvalue("slow-test-data") or
    not pytest.config.getvalue("atomic-dataset"),
    reason="--slow, --slow-test-data, or --atomic-dataset was not specified"
)
