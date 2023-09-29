import pytest
import numpy as np
from numpy.testing import assert_allclose
from tardis.plasma.properties import YgData


def test_exp1_times_exp(snapshot_np):
    x = np.array([499.0, 501.0, 710.0])
    # desired = np.array([0.00200000797, 0.0019920397, 0.0014064725])
    actual = YgData.exp1_times_exp(x)
    # assert_allclose(actual, desired)
    assert snapshot_np == actual
