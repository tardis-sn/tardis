from pathlib import Path

import pytest


@pytest.fixture
def artis_data_dir():
    return Path("tardis/io/model/artis/tests/data")
