import os

import numpy as np
import pandas as pd
import pytest

import tardis
from tardis.montecarlo.packet_source import BlackBodySimpleSource


@pytest.fixture
def data_path():
    return os.path.join(tardis.__path__[0], "montecarlo", "tests", "data")


@pytest.fixture
def packet_unit_test_fpath(tardis_ref_path):
    return os.path.abspath(os.path.join(tardis_ref_path, "packet_unittest.h5"))


@pytest.mark.array_compare(file_format='pdhdf')
def test_bb_packet_sampling():
    bb = BlackBodySimpleSource(2508)

    bb_df = pd.DataFrame(columns=['nu', 'mu', 'energy'])

    bb_df['nu'] = bb.create_blackbody_packet_nus(10000, 100)
    bb_df['mu'] =  bb.create_zero_limb_darkening_packet_mus(100)
    bb_df['energy'] = bb.create_uniform_packet_energies(100)

    return bb_df
