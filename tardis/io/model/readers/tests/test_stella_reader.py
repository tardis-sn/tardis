import numpy as np

from pathlib import Path

from pytest import fixture
from astropy import units as u
from tardis.io.model import read_stella_model

MODEL_DATA_PATH = Path(__file__).parent / "data"


@fixture
def stella_model_example_file1():
    return read_stella_model(MODEL_DATA_PATH / "mesa.stella.dat")


def test_read_stella_model_meta(stella_model_example_file1):
    """
    Test reading a STELLA model file
    """
    assert stella_model_example_file1.metadata["zones"] == 400
    np.testing.assert_almost_equal(
        stella_model_example_file1.metadata["t_max"].to(u.day).value, 50.0
    )
    np.testing.assert_almost_equal(
        stella_model_example_file1.metadata["inner_boundary_mass"]
        .to(u.g)
        .value,
        5.190242521200000e33,
    )
    np.testing.assert_almost_equal(
        stella_model_example_file1.metadata["total_mass"].to(u.g).value,
        2.618867335600000e34,
    )


def test_read_stella_model_data(stella_model_example_file1):
    """
    Test reading a STELLA model file
    """
    np.testing.assert_almost_equal(
        stella_model_example_file1.data.iloc[0, 0], 6.006769337200000e29
    )
    np.testing.assert_almost_equal(
        stella_model_example_file1.data.iloc[-1, -1], 2.123224906916000e04
    )
