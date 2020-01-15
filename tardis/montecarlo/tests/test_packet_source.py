import os

import numpy as np
import pandas as pd
import pytest

import tardis
from tardis.montecarlo.packet_source import BlackBodySimpleSource


@pytest.fixture
def data_path():
    return os.path.join(tardis.__path__[0], 'montecarlo', 'tests', 'data')


def test_bb_packet_sampling(tardis_ref_data):
    bb = BlackBodySimpleSource(2508)
    #ref_df = pd.read_hdf('test_bb_sampling.h5')
    ref_df = tardis_ref_data['/packet_unittest/blackbody']
    nus = bb.create_blackbody_packet_nus(10000, 100)
    mus = bb.create_zero_limb_darkening_packet_mus(100)
    unif_energies = bb.create_uniform_packet_energies(100)
    assert np.all(np.isclose(nus, ref_df['nus']))
    assert np.all(np.isclose(mus, ref_df['mus']))
    assert np.all(np.isclose(unif_energies, ref_df['energies']))