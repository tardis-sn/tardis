from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def example_csvy_file_dir():
    return Path(__file__).resolve().parent / "tests/data"
