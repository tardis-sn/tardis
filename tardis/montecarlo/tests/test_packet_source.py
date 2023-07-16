import os

import numpy as np
import pandas as pd
import pytest

from tardis.montecarlo.packet_source import BlackBodySimpleSource
from tardis.montecarlo import (
    montecarlo_configuration as montecarlo_configuration,
)


@pytest.fixture
def packet_unit_test_fpath(tardis_ref_path):
    return os.path.abspath(os.path.join(tardis_ref_path, "packet_unittest.h5"))


def test_bb_packet_sampling(request, tardis_ref_data, packet_unit_test_fpath):
    montecarlo_configuration.LEGACY_MODE_ENABLED = True
    bb = BlackBodySimpleSource(base_seed=1963, legacy_second_seed=2508)
    # ref_df = pd.read_hdf('test_bb_sampling.h5')
    if request.config.getoption("--generate-reference"):
        ref_bb = pd.read_hdf(packet_unit_test_fpath, key="/blackbody")
        ref_bb.to_hdf(
            tardis_ref_data, key="/packet_unittest/blackbody", mode="a"
        )
        pytest.skip("Reference data was generated during this run.")

    ref_df = tardis_ref_data["/packet_unittest/blackbody"]
    bb.temperature = 10000
    nus = bb.create_packet_nus(100)
    mus = bb.create_packet_mus(100)
    unif_energies = bb.create_packet_energies(100)
    assert np.all(np.isclose(nus, ref_df["nus"]))
    assert np.all(np.isclose(mus, ref_df["mus"]))
    assert np.all(np.isclose(unif_energies, ref_df["energies"]))
