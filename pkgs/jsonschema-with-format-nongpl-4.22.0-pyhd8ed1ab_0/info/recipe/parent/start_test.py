from pathlib import Path
import os
import sys

import pytest

HERE = Path(__file__).parent

os.environ.update(JSON_SCHEMA_TEST_SUITE=str(HERE / "json"))

PYTEST_ARGS = [
    "-vv",
    "--pyargs",
    "jsonschema",
    "--cov=jsonschema",
    "--cov-report=term-missing:skip-covered",
    "--no-cov-on-fail",
    *sys.argv[1:],
]

print("PYTEST_ARGS", *PYTEST_ARGS)

sys.exit(pytest.main(PYTEST_ARGS))
