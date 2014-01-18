import os
from astropy import units as u
from tardis import io
from numpy.testing import assert_almost_equal
import pytest

import numpy as np

test_data_directory = os.path.dirname(__file__)


def data_path(filename):
    data_dir = os.path.dirname(__file__)
    return os.path.join(data_dir, 'data', filename)


def test_simple_ascii_density_reader_time():
    time_model, index, v_inner, v_outer, density = io.read_simple_ascii_density(data_path('tardis_simple_ascii_density_test.dat'))

    assert time_model.unit.physical_type == 'time'
    assert_almost_equal(time_model.to(u.day).value, 1.0)

def test_simple_ascii_density_reader_data():

    time_model, index, v_inner, v_outer, density = io.read_simple_ascii_density(data_path('tardis_simple_ascii_density_test.dat'))
    assert index[0] == 0
    assert v_inner.unit == u.Unit('cm/s')

    assert_almost_equal(v_inner[3].value, 1.3e4*1e5)


def test_simple_ascii_abundance_reader():
    index, abundances = io.read_simple_ascii_abundances(data_path('artis_abundances.dat'))
    assert_almost_equal(abundances.ix[1, 0], 1.542953e-08)
    assert_almost_equal(abundances.ix[14, 54], 0.21864420000000001)

@pytest.mark.parametrize("v_inner_boundary, v_outer_boundary, actual_v_inner, actual_v_outer, inner_index, outer_index", [
    (0.0 * u.km/u.s, np.inf * u.km/u.s, np.nan, np.nan, None, None),
    (500 * u.km/u.s, 6000 * u.km/u.s, np.nan, 6000 * u.km/u.s, None, 16),
    (1300 * u.km/u.s, 6000 * u.km/u.s, 1300 * u.km/u.s, 6000 * u.km/u.s, 0, 16),
    (1600 * u.km/u.s, 6000 * u.km/u.s, 1600 * u.km/u.s, 6000 * u.km/u.s, 1, 16)
])
def test_ascii_reader_density_boundaries(v_inner_boundary, v_outer_boundary, actual_v_inner, actual_v_outer,
                                         inner_index, outer_index):
    v_inner, v_outer, mean_densities, inner_boundary_index, outer_boundary_index = \
        io.read_density_file(data_path('artis_model.dat'), 'artis', 19 * u.day, v_inner_boundary, v_outer_boundary)

    assert inner_boundary_index == inner_index
    assert outer_boundary_index == outer_index

    if not np.isnan(actual_v_inner):
        assert_almost_equal(v_inner[0], actual_v_inner)

    if not np.isnan(actual_v_outer):
        assert_almost_equal(v_outer[-1], actual_v_outer)

