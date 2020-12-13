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


def test_bb_packet_sampling(request, tardis_ref_data, packet_unit_test_fpath):
    bb = BlackBodySimpleSource(2508)
    rng = np.random.default_rng(seed=1963)
    # ref_df = pd.read_hdf('test_bb_sampling.h5')
    if request.config.getoption("--generate-reference"):
        ref_bb = pd.read_hdf(packet_unit_test_fpath, key="/blackbody")
        ref_bb.to_hdf(
            tardis_ref_data, key="/packet_unittest/blackbody", mode="a"
        )
        pytest.skip("Reference data was generated during this run.")

    ref_df = tardis_ref_data["/packet_unittest/blackbody"]
    nus = bb.create_blackbody_packet_nus(10000, 100, rng)
    mus = bb.create_zero_limb_darkening_packet_mus(100, rng)
    unif_energies = bb.create_uniform_packet_energies(100, rng)
    assert np.all(np.isclose(nus, ref_df["nus"]))
    assert np.all(np.isclose(mus, ref_df["mus"]))
    assert np.all(np.isclose(unif_energies, ref_df["energies"]))
