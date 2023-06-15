import pytest
import numpy as np

from astropy import units as u

from tardis.io.model import read_stella_model

# functions to use pytest to test read_stella_model from an example file


def test_read_stella_model():
    """
    Test reading a STELLA model file
    """
    fname = "tardis/io/model/tests/data/stella_model.dat"
    model = read_stella_model(fname)
    assert model.metadata["zones"] == 400
    np.testing.assert_almost_equal(model.metadata["t_max"].values, 50.0)
    np.testing.assert_almost_equal(model.metadata["inner_boundary_mass"].values, 5.190242521200000E+33 * u.g)
    np.testing.assert_almost_equal(model.metadata["total_mass"].values, 2.618867335600000E+34 * u.g)