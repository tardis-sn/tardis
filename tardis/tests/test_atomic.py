__author__ = 'maryam'

#this is the test environment for atomic

from tardis import atomic

from numpy import testing

def test_atomic_h5_readin():
    data = atomic.read_atomic_data()
    assert data['symbol'][25] == "Fe"



def test_ionization_h5_readin():
    data = atomic.read_ionization_data()
    testing.assert_almost_equal(data['ionization_energy'][0],13.)