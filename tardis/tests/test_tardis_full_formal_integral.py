from pathlib import Path

import numpy.testing as npt
import pytest
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose

from tardis.io.configuration.config_reader import Configuration
from tardis.io.util import HDFWriterMixin
from tardis.simulation.base import Simulation
from tardis.tests.fixtures.regression_data import RegressionData

config_line_modes = ["downbranch", "macroatom"]
interpolate_shells = [-1, 30]


@pytest.fixture(scope="module", params=config_line_modes)
def base_config(request, example_configuration_dir: Path):
    config = Configuration.from_yaml(
        example_configuration_dir / "tardis_configv1_verysimple.yml"
    )
    config["plasma"]["line_interaction_type"] = request.param
    config["montecarlo"]["no_of_packets"] = 4.0e4
    config["montecarlo"]["last_no_of_packets"] = 1.0e5
    config["montecarlo"]["no_of_virtual_packets"] = 0
    config["spectrum"]["method"] = "integrated"
    config["spectrum"]["integrated"]["points"] = 200
    print("config", config)

    return config


@pytest.fixture(scope="module", params=interpolate_shells)
def config(base_config, request):
    base_config["spectrum"]["integrated"]["interpolate_shells"] = request.param
    return base_config


class SimulationContainer(HDFWriterMixin):
    hdf_properties = ["spectrum_solver", "transport"]

    def __init__(self, simulation):
        self.spectrum_solver = simulation.spectrum_solver
        self.transport = simulation.transport


class TestTransportSimpleFormalIntegral:
    """
    Very simple run with the formal integral spectral synthesis method
    """

    _name = "test_transport_simple_integral"

    @pytest.fixture(scope="class")
    def simulation(self, config, atomic_data_fname):
        config.atom_data = atomic_data_fname

        self.name = self._name + f"_{config.plasma.line_interaction_type:s}"
        if config.spectrum.integrated.interpolate_shells > 0:
            self.name += "_interp"

        simulation = Simulation.from_config(config)
        simulation.run_convergence()
        simulation.run_final()
        simulation.spectrum_solver.hdf_properties = [
            "spectrum_real_packets",
            "spectrum_integrated",
        ]
        simulation.transport.hdf_properties = ["transport_state"]

        return simulation

    def test_simulation(self, simulation, request):
        regression_data = RegressionData(request)
        container = SimulationContainer(simulation)
        regression_data.sync_hdf_store(container)

    def test_j_blue_estimators(self, simulation, request):
        regression_data = RegressionData(request)
        j_blue_estimator = (
            simulation.transport.transport_state.radfield_mc_estimators.j_blue_estimator
        )
        expected = regression_data.sync_ndarray(j_blue_estimator)
        npt.assert_allclose(j_blue_estimator, expected)

    def test_spectrum(self, simulation, request):
        regression_data = RegressionData(request)
        luminosity = simulation.spectrum_solver.spectrum_real_packets.luminosity
        expected = regression_data.sync_ndarray(luminosity.cgs.value)
        expected = u.Quantity(expected, "erg /s")
        assert_quantity_allclose(luminosity, expected)

    def test_spectrum_integrated(self, simulation, request):
        regression_data = RegressionData(request)
        luminosity = simulation.spectrum_solver.spectrum_integrated.luminosity
        expected = regression_data.sync_ndarray(luminosity.cgs.value)
        expected = u.Quantity(expected, "erg /s")
        assert_quantity_allclose(luminosity, expected)
