from .. import photon
from numpy import testing

def test_blackbody():
    testing.assert_almost_equal(photon.blackbody_nu(3e5/(1e-13*5000), 10000), 0.00018923672429264153)