from copy import deepcopy

import numpy.testing as npt
import pytest

from tardis.simulation import Simulation


@pytest.fixture(scope="function")
def weighted_transport_state(
    montecarlo_main_loop_config,
    atomic_dataset,
    simple_weighted_packet_source,
):
    atomic_dataset = deepcopy(atomic_dataset)

    simulation = Simulation.from_config(
        montecarlo_main_loop_config,
        atom_data=atomic_dataset,
        virtual_packet_logging=False,
        legacy_mode_enabled=True,
    )

    simulation.packet_source = simple_weighted_packet_source
    simulation.run_convergence()
    simulation.run_final()

    return simulation.transport.transport_state


@pytest.mark.parametrize(
    "attr_name, getter",
    [
        ("output_energies", lambda ts: ts.packet_collection.output_energies),
        ("output_nus", lambda ts: ts.packet_collection.output_nus),
        ("nu_bar_estimator", lambda ts: ts.radfield_mc_estimators.nu_bar_estimator),
        ("j_estimator", lambda ts: ts.radfield_mc_estimators.j_estimator),
    ],
    ids=[
        "output_energies",
        "output_nus",
        "nu_bar_estimator",
        "j_estimator",
    ],
)
def test_montecarlo_main_loop_weighted(
    weighted_transport_state,
    regression_data,
    attr_name,
    getter,
):
    actual_value = getter(weighted_transport_state)

    regression_data.fname = f"test_montecarlo_main_loop_weighted__{attr_name}"
    expected_value = regression_data.sync_ndarray(actual_value)

    npt.assert_allclose(
        actual_value,
        expected_value,
        rtol=1e-2,
        err_msg=f"Mismatch in {attr_name}",
    )
