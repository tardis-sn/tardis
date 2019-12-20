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
    packets = bb.create_blackbody_packet_nus(10000, 100)
    reference_packets = np.load(os.path.join(data_path, 'mc_packets_100_t10000.npy'))
    assert np.all(np.isclose(packets, reference_packets))

def test_uniform_packet_energies(data_path):
    bb = BlackBodySimpleSource(2508)
    unif_energies = bb.create_uniform_packet_energies(100)
    ref_energies = 0.01 * np.ones(100)
    assert np.all(np.isclose(unif_energies, ref_energies))

