import pytest
import pandas as pd
import os
import numpy.testing as npt
import numpy as np
import copy

import tardis.montecarlo.montecarlo_numba.base as base
from tardis.montecarlo.montecarlo_numba.numba_interface import PacketCollection, VPacketCollection
from tardis.montecarlo import montecarlo_configuration as montecarlo_configuration
import tardis.montecarlo.montecarlo_numba.numba_interface as numba_interface
import tardis.montecarlo.montecarlo_numba.numba_config as numba_config
import tardis.montecarlo.montecarlo_numba.r_packet as r_packet
import tardis.montecarlo.montecarlo_numba.single_packet_loop as spl


@pytest.mark.xfail(reason='To be implemented')
def test_montecarlo_radial1d():
    assert False


def test_montecarlo_main_loop(nb_simulation_verysimple, tardis_ref_path, tmpdir):
    fname = str(tmpdir.mkdir("data").join("test.hdf"))
    C_fname = os.path.join(tardis_ref_path, "montecarlo_one_packet_compare_data.h5")

    sim = nb_simulation_verysimple
    sim.to_hdf(fname)

    actual_nu = pd.read_hdf(fname, key="/simulation/runner/output_nu").values
    expected_nu = pd.read_hdf(C_fname, key="/simulation/runner/output_nu").values

    actual_energy = pd.read_hdf(fname, key="/simulation/runner/output_energy").values
    expected_energy = pd.read_hdf(C_fname, key="/simulation/runner/output_energy").values

    npt.assert_array_almost_equal(actual_nu, expected_nu)
    npt.assert_array_almost_equal(actual_energy, expected_energy)