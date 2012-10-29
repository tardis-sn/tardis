__author__ = 'maryam'

#this is the test environment for atomic

from tardis import atomic

def test_atomic_h5_readin():
    data = atomic.read_atomic_data()
    assert data['symbol'][25] == "Fe"



