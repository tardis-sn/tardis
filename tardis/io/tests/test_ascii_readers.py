import os
from astropy import units as u
from tardis import io
from numpy.testing import assert_almost_equal

test_data_directory = os.path.join(os.path.dirname(__file__), 'data')

def test_simple_ascii_density_reader():
    time_model, data = io.read_simple_ascii_density(os.path.join(test_data_directory, 'tardis_simple_ascii_density_test.dat'))

    assert time_model.unit.physical_type == 'time'
    assert_almost_equal(time_model.to(u.day).value, 1.0)

    assert data['index'][0] == 0

    assert_almost_equal(data['velocity'][3], 1.3e4)



