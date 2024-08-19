import astropy.units as u
import numpy as np
import pandas as pd
import pytest

import tardis
from tardis.io.configuration.config_reader import Configuration
from tardis.simulation import Simulation


@pytest.fixture(scope="module")
def config(example_configuration_dir):
    return Configuration.from_yaml(
        example_configuration_dir / "tardis_configv1_verysimple.yml"
    )


@pytest.fixture(scope="module")
def simulation_one_loop(config, atomic_data_fname):
    config.atom_data = atomic_data_fname
    config.montecarlo.iterations = 2
    config.montecarlo.no_of_packets = int(4e4)
    config.montecarlo.last_no_of_packets = int(4e4)

    sim = Simulation.from_config(config)
    sim.run_convergence()
    sim.run_final()
    return sim


@pytest.mark.parametrize(
    "attr",
    [
        "iterations_w",
        "iterations_t_rad",
        "iterations_electron_densities",
        "iterations_t_inner",
    ],
)
def test_plasma_state_iterations(simulation_one_loop, attr, regression_data):
    actual = getattr(simulation_one_loop, attr)
    if hasattr(actual, "value"):
        actual = actual.value
    actual = pd.DataFrame(actual)
    expected = regression_data.sync_dataframe(actual)
    pd.testing.assert_frame_equal(actual, expected, rtol=1e-5, atol=1e-8)


@pytest.mark.parametrize(
    "attr",
    [
        "nu_bar_estimator",
        "j_estimator",
        "t_radiative",
        "dilution_factor",
        "output_nus",
        "output_energies",
    ],
)
def test_plasma_estimates(simulation_one_loop, attr, regression_data):
    if attr in ["nu_bar_estimator", "j_estimator"]:
        actual = getattr(
            simulation_one_loop.transport.transport_state.radfield_mc_estimators,
            attr,
        )
    elif attr in ["t_radiative", "dilution_factor"]:
        actual = getattr(simulation_one_loop.simulation_state, attr)
    elif attr in ["output_nus", "output_energies"]:
        actual = getattr(
            simulation_one_loop.transport.transport_state.packet_collection,
            attr,
        )
    else:
        actual = getattr(simulation_one_loop.transport, attr)

    if hasattr(actual, "value"):
        actual = actual.value
    actual = pd.Series(actual)
    expected = regression_data.sync_dataframe(actual)
    pd.testing.assert_series_equal(actual, expected, rtol=1e-5, atol=1e-8)


@pytest.fixture(scope="module")
def simulation_without_loop(atomic_data_fname, config):
    config.atom_data = atomic_data_fname
    config.montecarlo.iterations = 2
    return Simulation.from_config(config)


def test_plasma_state_storer_store(
    atomic_data_fname, config, simulation_without_loop
):
    simulation = simulation_without_loop

    w_test = np.linspace(0, 1, 20)
    t_rad_test = np.linspace(12000, 9000, 20) * u.K
    electron_densities_test = pd.Series(np.linspace(1e7, 1e6, 20))
    t_inner_test = 12500 * u.K

    simulation.store_plasma_state(
        1, w_test, t_rad_test, electron_densities_test, t_inner_test
    )

    np.testing.assert_allclose(simulation.iterations_w[1, :], w_test)
    np.testing.assert_allclose(simulation.iterations_t_rad[1, :], t_rad_test)
    np.testing.assert_allclose(
        simulation.iterations_electron_densities[1, :], electron_densities_test
    )
    np.testing.assert_allclose(simulation.iterations_t_inner[1], t_inner_test)


def test_plasma_state_storer_reshape(
    atomic_data_fname, config, simulation_without_loop
):
    simulation = simulation_without_loop
    simulation.reshape_plasma_state_store(0)

    assert simulation.iterations_t_rad.shape == (1, 20)
    assert simulation.iterations_w.shape == (1, 20)
    assert simulation.iterations_electron_densities.shape == (1, 20)
    assert simulation.iterations_t_inner.shape == (1,)


def test_version_tag(simulation_without_loop):
    assert simulation_without_loop.version == tardis.__version__
