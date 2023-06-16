import numpy as np

from astropy import units as u

from pathlib import Path

from tardis.io.model import read_stella_model

MODEL_DATA_PATH = Path(__file__).parent / "data"


def test_read_stella_model():
    """
    Test reading a STELLA model file
    """
    fname = MODEL_DATA_PATH / "mesa.stella.dat"
    model = read_stella_model(fname)
    assert model.metadata["zones"] == 400
    np.testing.assert_almost_equal(
        model.metadata["t_max"].to(u.day).value, 50.0
    )
    np.testing.assert_almost_equal(
        model.metadata["inner_boundary_mass"].to(u.g).value,
        5.190242521200000e33,
    )
    np.testing.assert_almost_equal(
        model.metadata["total_mass"].to(u.g).value, 2.618867335600000e34
    )
