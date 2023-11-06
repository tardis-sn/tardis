import numpy as np

from tardis.plasma.properties import YgData


def test_exp1_times_exp(snapshot_np):
    x = np.array([499.0, 501.0, 710.0])
    actual = YgData.exp1_times_exp(x)
    assert snapshot_np == actual
