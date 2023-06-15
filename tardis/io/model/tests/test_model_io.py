import pytest

from tardis.io.model import read_stella_model

# functions to use pytest to test read_stella_model from an example file


def test_read_stella_model():
    """
    Test reading a STELLA model file
    """
    fname = "tardis/io/model/tests/data/stella_model.dat"
    model = read_stella_model(fname)
    assert model.header["t_max"].value == 0.0
    assert model.header["zones"] == 100
    assert model.header["inner_boundary_mass"].value == 1.0e33
    assert model.header["total_mass"].value == 1.0e33
