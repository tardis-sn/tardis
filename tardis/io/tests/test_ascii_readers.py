import os
from astropy import units as u
from tardis import io
from numpy.testing import assert_almost_equal

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




