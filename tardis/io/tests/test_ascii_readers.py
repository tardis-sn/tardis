import os
from astropy import units as u
from tardis import io

import numpy.testing as npt

import pytest

import numpy as np

test_data_directory = os.path.dirname(__file__)


def data_path(filename):
    data_dir = os.path.dirname(__file__)
    return os.path.join(data_dir, "data", filename)


def test_simple_ascii_density_reader_time():
    time_model, velocity, density = io.read_simple_ascii_density(
        data_path("tardis_simple_ascii_density_test.dat")
    )

    assert time_model.unit.physical_type == "time"
    npt.assert_almost_equal(time_model.to(u.day).value, 1.0)


def test_simple_ascii_density_reader_data():

    time_model, velocity, density = io.read_simple_ascii_density(
        data_path("tardis_simple_ascii_density_test.dat")
    )
    assert velocity.unit == u.Unit("cm/s")

    npt.assert_allclose(velocity[3].value, 1.3e4 * 1e5)


def test_simple_ascii_abundance_reader():
    index, abundances = io.read_simple_ascii_abundances(
        data_path("artis_abundances.dat")
    )
    npt.assert_almost_equal(abundances.loc[1, 0], 1.542953e-08)
    npt.assert_almost_equal(abundances.loc[14, 54], 0.21864420000000001)


def test_ascii_reader_invalid_volumes():
    with pytest.raises(io.model_reader.ConfigurationError):
        io.read_density_file(data_path("invalid_artis_model.dat"), "artis")
