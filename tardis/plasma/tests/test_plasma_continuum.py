import numpy as np

import numpy.testing as npt
from tardis.plasma.properties import YgData


def test_exp1_times_exp(regression_data):
    x = np.array([499.0, 501.0, 710.0])
    actual = YgData.exp1_times_exp(x)
    expected = regression_data.sync_ndarray(actual)
    npt.assert_allclose(actual, expected)
