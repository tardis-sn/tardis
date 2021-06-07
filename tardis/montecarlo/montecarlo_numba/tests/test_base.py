import pytest
import pandas as pd
import os
import numpy.testing as npt
import numpy as np
from copy import deepcopy
from astropy import units as u

from tardis.montecarlo import (
    montecarlo_configuration as montecarlo_configuration,
)
import tardis.montecarlo.montecarlo_numba.base as base
from tardis.simulation import Simulation
import tardis.montecarlo.montecarlo_numba.r_packet as r_packet
import tardis.montecarlo.montecarlo_numba.single_packet_loop as spl
from tardis.montecarlo.montecarlo_numba.numba_interface import (
    PacketCollection,
    VPacketCollection,
    NumbaModel,
    numba_plasma_initialize,
    Estimators,
    configuration_initialize,
)


@pytest.mark.xfail(reason="To be implemented")
def test_montecarlo_radial1d():
    assert False

# Currently, this essentially tests the full simulation, would probably be better to
# test individual elements
def test_montecarlo_main_loop(
    config_verysimple,
    atomic_dataset,
    tardis_ref_path,
    tmpdir,
    set_seed_fixture,
    random_call_fixture,
    request
):

    montecarlo_configuration.LEGACY_MODE_ENABLED = True
    # Setup model config from verysimple
    atomic_data = deepcopy(atomic_dataset)
    config_verysimple.montecarlo.last_no_of_packets = 1e5
    config_verysimple.montecarlo.no_of_virtual_packets = 0
    config_verysimple.montecarlo.iterations = 1
    config_verysimple.montecarlo.single_packet_seed = 0
    config_verysimple.plasma.line_interaction_type = 'macroatom'
    del config_verysimple["config_dirname"]

    sim = Simulation.from_config(config_verysimple, atom_data=atomic_data)
    sim.run()


    compare_fname = os.path.join(tardis_ref_path, "montecarlo_1e5_compare_data.h5")
    if request.config.getoption("--generate-reference"):
        sim.to_hdf(compare_fname, overwrite=True)

    # Load compare data from refdata
    else:
        expected_nu = pd.read_hdf(
            compare_fname, key="/simulation/runner/output_nu"
        ).values
        expected_energy = pd.read_hdf(
            compare_fname, key="/simulation/runner/output_energy"
        ).values
        expected_nu_bar_estimator = pd.read_hdf(
            compare_fname, key="/simulation/runner/nu_bar_estimator"
        ).values
        expected_j_estimator = pd.read_hdf(
            compare_fname, key="/simulation/runner/j_estimator"
        ).values


        actual_energy = sim.runner.output_energy
        actual_nu = sim.runner.output_nu
        actual_nu_bar_estimator = sim.runner.nu_bar_estimator
        actual_j_estimator = sim.runner.j_estimator

        # Compare
        npt.assert_allclose(
            actual_nu_bar_estimator, expected_nu_bar_estimator, rtol=1e-13
        )
        npt.assert_allclose(actual_j_estimator, expected_j_estimator, rtol=1e-13)
        npt.assert_allclose(actual_energy.value, expected_energy, rtol=1e-13)
        npt.assert_allclose(actual_nu.value, expected_nu, rtol=1e-13)
