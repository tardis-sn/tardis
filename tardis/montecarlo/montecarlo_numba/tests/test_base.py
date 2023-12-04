import pytest
import pandas as pd
import os
import numpy.testing as npt
from copy import deepcopy
from tardis.base import run_tardis
from pandas.testing import assert_frame_equal

from tardis.montecarlo import (
    montecarlo_configuration as montecarlo_configuration,
)
from tardis.simulation import Simulation


@pytest.mark.xfail(reason="To be implemented")
def test_montecarlo_radial1d():
    assert False


@pytest.fixture(scope="function")
def montecarlo_main_loop_config(
    config_montecarlo_1e5_verysimple,
    atomic_dataset,
):
    montecarlo_configuration.LEGACY_MODE_ENABLED = True
    # Setup model config from verysimple
    atomic_data = deepcopy(atomic_dataset)
    config_montecarlo_1e5_verysimple.montecarlo.last_no_of_packets = 1e5
    config_montecarlo_1e5_verysimple.montecarlo.no_of_virtual_packets = 0
    config_montecarlo_1e5_verysimple.montecarlo.iterations = 1
    config_montecarlo_1e5_verysimple.plasma.line_interaction_type = "macroatom"

    del config_montecarlo_1e5_verysimple["config_dirname"]
    return config_montecarlo_1e5_verysimple


def test_montecarlo_main_loop(
    montecarlo_main_loop_config,
    tardis_ref_path,
    request,
    atomic_dataset,
):
    montecarlo_main_loop_simulation = Simulation.from_config(
        montecarlo_main_loop_config,
        atom_data=atomic_dataset,
        virtual_packet_logging=True,
    )
    montecarlo_main_loop_simulation.run_convergence()
    montecarlo_main_loop_simulation.run_final()

    compare_fname = os.path.join(
        tardis_ref_path, "montecarlo_1e5_compare_data.h5"
    )
    if request.config.getoption("--generate-reference"):
        montecarlo_main_loop_simulation.to_hdf(compare_fname, overwrite=True)

    # Load compare data from refdata
    expected_nu = pd.read_hdf(
        compare_fname, key="/simulation/transport/output_nu"
    ).values
    expected_energy = pd.read_hdf(
        compare_fname, key="/simulation/transport/output_energy"
    ).values
    expected_nu_bar_estimator = pd.read_hdf(
        compare_fname, key="/simulation/transport/nu_bar_estimator"
    ).values
    expected_j_estimator = pd.read_hdf(
        compare_fname, key="/simulation/transport/j_estimator"
    ).values

    actual_energy = montecarlo_main_loop_simulation.transport.output_energy
    actual_nu = montecarlo_main_loop_simulation.transport.output_nu
    actual_nu_bar_estimator = (
        montecarlo_main_loop_simulation.transport.nu_bar_estimator
    )
    actual_j_estimator = montecarlo_main_loop_simulation.transport.j_estimator

    # Compare
    npt.assert_allclose(
        actual_nu_bar_estimator, expected_nu_bar_estimator, rtol=1e-13
    )
    npt.assert_allclose(actual_j_estimator, expected_j_estimator, rtol=1e-13)
    npt.assert_allclose(actual_energy.value, expected_energy, rtol=1e-13)
    npt.assert_allclose(actual_nu.value, expected_nu, rtol=1e-13)


def donot_test_montecarlo_main_loop_vpacket_log(
    montecarlo_main_loop_simulation,
    tardis_ref_path,
    request,
):
    montecarlo_main_loop_simulation = Simulation.from_config(
        montecarlo_main_loop_config,
        atom_data=atomic_data,
        virtual_packet_logging=False,
    )

    montecarlo_main_loop_simulation.run_convergence()
    montecarlo_main_loop_simulation.run_final()

    compare_fname = os.path.join(
        tardis_ref_path, "montecarlo_1e5_compare_data.h5"
    )
    if request.config.getoption("--generate-reference"):
        montecarlo_main_loop_simulation.to_hdf(compare_fname, overwrite=True)

    # Load compare data from refdata
    expected_nu = pd.read_hdf(
        compare_fname, key="/simulation/transport/output_nu"
    ).values
    expected_energy = pd.read_hdf(
        compare_fname, key="/simulation/transport/output_energy"
    ).values
    expected_nu_bar_estimator = pd.read_hdf(
        compare_fname, key="/simulation/transport/nu_bar_estimator"
    ).values
    expected_j_estimator = pd.read_hdf(
        compare_fname, key="/simulation/transport/j_estimator"
    ).values

    actual_energy = montecarlo_main_loop_simulation.transport.output_energy
    actual_nu = montecarlo_main_loop_simulation.transport.output_nu
    actual_nu_bar_estimator = (
        montecarlo_main_loop_simulation.transport.nu_bar_estimator
    )
    actual_j_estimator = montecarlo_main_loop_simulation.transport.j_estimator

    # Compare
    npt.assert_allclose(
        actual_nu_bar_estimator, expected_nu_bar_estimator, rtol=1e-13
    )
    npt.assert_allclose(actual_j_estimator, expected_j_estimator, rtol=1e-13)
    npt.assert_allclose(actual_energy.value, expected_energy, rtol=1e-13)
    npt.assert_allclose(actual_nu.value, expected_nu, rtol=1e-13)
