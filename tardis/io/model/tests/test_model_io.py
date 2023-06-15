import pytest

from tardis.io.model import read_stella_model

# functions to use pytest to test read_stella_model from an example file


def test_read_stella_model():
    """
    Test reading a STELLA model file
    """
    fname = "tardis/io/model/tests/data/stella_model.dat"
    model = read_stella_model(fname)
    assert model.metadata["t_max"].value == 0.0
    assert model.metadata["zones"] == 100
    assert model.metadata["inner_boundary_mass"].value == 1.0e33
    assert model.metadata["total_mass"].value == 1.0e33
