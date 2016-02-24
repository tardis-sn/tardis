import os

import numpy as np
import pytest

import tardis
from tardis.montecarlo.packet_source import BlackBodySimpleSource

@pytest.fixture
def data_path():
    return os.path.join(tardis.__path__[0], 'montecarlo', 'tests', 'data')




def test_bb_packet_source(data_path):
    bb = BlackBodySimpleSource(2508)
    packets = bb.create_packet_nus(10000, 100)
    reference_packets = np.load(os.path.join(data_path, 'mc_packets_100_t10000.npy'))
    assert np.all(np.isclose(packets, reference_packets))


