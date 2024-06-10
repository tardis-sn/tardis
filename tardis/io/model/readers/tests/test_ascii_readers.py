import os
from astropy import units as u
from tardis import io
from tardis.io.model.readers.generic_readers import ConfigurationError

import numpy.testing as npt

import pytest

import numpy as np

from pathlib import Path
from tardis.io.model.readers.base import read_density_file

from tardis.io.model.readers.generic_readers import (
    read_simple_ascii_mass_fractions,
    read_simple_ascii_density,
)

test_data_directory = os.path.dirname(__file__)


def test_simple_ascii_density_reader_time(example_model_file_dir):
    (time_model, velocity, density,) = read_simple_ascii_density(
        example_model_file_dir / "tardis_simple_ascii_density_test.dat"
    )

    assert time_model.unit.physical_type == "time"
    npt.assert_almost_equal(time_model.to(u.day).value, 1.0)
    npt.assert_allclose(velocity[3].value, 1.3e4 * 1e5)
    assert velocity.unit == u.Unit("cm/s")


def test_simple_ascii_abundance_reader(example_model_file_dir):
    (index, abundances,) = read_simple_ascii_mass_fractions(
        example_model_file_dir / "artis_abundances.dat"
    )

    npt.assert_almost_equal(abundances.loc[1, 0], 1.542953e-08)
    npt.assert_almost_equal(abundances.loc[14, 54], 0.21864420000000001)


def test_ascii_reader_invalid_volumes(example_model_file_dir):
    with pytest.raises(ConfigurationError):
        read_density_file(
            example_model_file_dir / "invalid_artis_model.dat", "artis"
        )
