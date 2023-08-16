import pytest
from pathlib import Path


@pytest.fixture
def example_model_file_dir():
    return (
        Path(__file__).resolve().parent / "model" / "readers" / "tests" / "data"
    )
