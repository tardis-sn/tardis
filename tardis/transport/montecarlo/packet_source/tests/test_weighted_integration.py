from copy import deepcopy

import numpy.testing as npt
import pandas as pd
import pytest

from tardis.simulation import Simulation


@pytest.mark.xfail(reason="Relies on later test")
def test_montecarlo_main_loop_weighted(
    montecarlo_main_loop_config,
    regression_data,
    atomic_dataset,
    simple_weighted_packet_source,
):
    atomic_dataset = deepcopy(atomic_dataset)
    montecarlo_main_loop_simulation_weighted = Simulation.from_config(
        montecarlo_main_loop_config,
        atom_data=atomic_dataset,
        virtual_packet_logging=False,
        legacy_mode_enabled=True,
    )
    montecarlo_main_loop_simulation_weighted.packet_source = (
        simple_weighted_packet_source
    )
    montecarlo_main_loop_simulation_weighted.run_convergence()
    montecarlo_main_loop_simulation_weighted.run_final()

    # Get the montecarlo simple regression data
    # Directly use the regression data path to construct the path
    regression_data_dir = (
        regression_data.regression_data_path
        / "tardis/transport/montecarlo/tests/test_montecarlo_main_loop/test_montecarlo_main_loop.h5"
    )
    expected_hdf_store = pd.HDFStore(regression_data_dir, mode="r")

    # Load compare data from refdata

    expected_nu = expected_hdf_store[
        "/simulation/transport/transport_state/output_nu"
    ]
    expected_energy = expected_hdf_store[
        "/simulation/transport/transport_state/output_energy"
    ]
    expected_nu_bar_estimator = expected_hdf_store[
        "/simulation/transport/transport_state/nu_bar_estimator"
    ]
    expected_j_estimator = expected_hdf_store[
        "/simulation/transport/transport_state/j_estimator"
    ]
    expected_hdf_store.close()
    transport_state = (
        montecarlo_main_loop_simulation_weighted.transport.transport_state
    )
    actual_energy = transport_state.packet_collection.output_energies
    actual_nu = transport_state.packet_collection.output_nus
    actual_nu_bar_estimator = transport_state.estimators_bulk.mean_frequency
    actual_j_estimator = transport_state.estimators_bulk.mean_intensity_total

    # Compare
    npt.assert_allclose(
        actual_nu_bar_estimator, expected_nu_bar_estimator, rtol=1e-2
    )
    npt.assert_allclose(actual_j_estimator, expected_j_estimator, rtol=1e-2)
    npt.assert_allclose(actual_energy, expected_energy, rtol=1e-2)
    npt.assert_allclose(actual_nu, expected_nu, rtol=1e-2)
