import numpy as np

from pathlib import Path

from pytest import fixture
from astropy import units as u
from tardis.io.model.readers.cmfgen import read_cmfgen_model

MODEL_DATA_PATH = Path(__file__).parent / "data"


@fixture
def cmfgen_model_example_file():
    return read_cmfgen_model(MODEL_DATA_PATH / "cmfgen_model.csv")


def test_read_cmfgen_model_meta(cmfgen_model_example_file):
    """
    Test reading a CMFGEN model file
    """
    metadata = cmfgen_model_example_file.metadata
    assert set(metadata.keys()).issubset(
        {
            "t0",
            "velocity_unit",
            "temperature_unit",
            "densities_unit",
            "electron_densities_unit",
        }
    )
    np.testing.assert_almost_equal(metadata["t0"].value, 0.976)


def test_read_cmfgen_model_data(cmfgen_model_example_file):
    """
    Test reading a cmfgen model file
    """
    data = cmfgen_model_example_file.data
    np.testing.assert_almost_equal(data.iloc[0, 0], 871.66905)
