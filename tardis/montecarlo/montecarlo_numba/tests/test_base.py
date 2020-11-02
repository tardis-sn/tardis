import pytest
import pandas as pd
import os
import numpy.testing as npt
import numpy as np
from copy import deepcopy

from tardis.simulation import Simulation


@pytest.mark.xfail(reason='To be implemented')
def test_montecarlo_radial1d():
    assert False


def test_montecarlo_main_loop(config_verysimple, atomic_dataset, tardis_ref_path, tmpdir):
    fname = str(tmpdir.mkdir("data").join("test.hdf"))
    C_fname = os.path.join(tardis_ref_path, "montecarlo_one_packet_compare_data.h5")

    atomic_data = deepcopy(atomic_dataset)
    sim = Simulation.from_config(config_verysimple, atom_data=atomic_data)
    sim.iterate(100000)
    sim.to_hdf(fname)

    actual_nu = pd.read_hdf(fname, key="/simulation/runner/output_nu").values
    expected_nu = pd.read_hdf(C_fname, key="/simulation/runner/output_nu").values

    actual_energy = pd.read_hdf(fname, key="/simulation/runner/output_energy").values
    expected_energy = pd.read_hdf(C_fname, key="/simulation/runner/output_energy").values

    npt.assert_array_almost_equal(actual_nu, expected_nu)
    npt.assert_array_almost_equal(actual_energy, expected_energy)