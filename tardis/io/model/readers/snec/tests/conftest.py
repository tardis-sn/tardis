import pytest
from pathlib import Path


@pytest.fixture
def minimal_xg_file():
    return Path(__file__).parent / "data" / "minimal.xg"

@pytest.fixture
def minimal_profile_file():
    return Path(__file__).parent / "data" / "minimal_profile.dat"
