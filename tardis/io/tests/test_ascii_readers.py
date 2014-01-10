import os
from astropy import units as u
from tardis import io
from numpy.testing import assert_almost_equal

test_data_directory = os.path.dirname(__file__)


def data_path(filename):
    data_dir = os.path.dirname(__file__)
    return os.path.join(data_dir, 'data', filename)


def test_simple_ascii_density_reader():
    time_model, data = io.read_simple_ascii_density(data_path('tardis_simple_ascii_density_test.dat'))

    assert time_model.unit.physical_type == 'time'
    assert_almost_equal(time_model.to(u.day).value, 1.0)

    assert data['index'][0] == 0

    assert_almost_equal(data['velocity'][3], 1.3e4)



