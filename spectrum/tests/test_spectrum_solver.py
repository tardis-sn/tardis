from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
import pytest
from astropy.tests.helper import assert_quantity_allclose

from tardis.io.configuration.config_reader import Configuration
from tardis.simulation.base import Simulation
from tardis.spectrum.base import SpectrumSolver
from tardis.tests.fixtures.regression_data import RegressionData


class TestSpectrumSolver:
    regression_data: RegressionData = None

    @pytest.fixture(scope="class")
    def simulation(
        self,
        request,
        atomic_data_fname,
        generate_reference,
        example_configuration_dir: Path,
    ):
        config = Configuration.from_yaml(
            str(example_configuration_dir / "tardis_configv1_verysimple.yml")
        )
        config["atom_data"] = atomic_data_fname

        simulation = Simulation.from_config(config)
        simulation.run_final()

        request.cls.regression_data = RegressionData(request)
        data = request.cls.regression_data.sync_hdf_store(simulation)

        yield simulation
        data.close()

    def get_expected_data(self, key: str):
        return pd.read_hdf(self.regression_data.fpath, key)

    def test_initialization(self, simulation):
        transport_state = simulation.transport.transport_state
        spectrum_frequency_grid = simulation.transport.spectrum_frequency_grid

        solver = SpectrumSolver(transport_state, spectrum_frequency_grid, None)
        assert solver.transport_state == transport_state
        assert np.array_equal(
            solver.spectrum_frequency_grid.value, spectrum_frequency_grid.value
        )
        assert np.array_equal(
            solver._montecarlo_virtual_luminosity.value,
            np.zeros_like(spectrum_frequency_grid.value),
        )
        assert solver._integrator is None
        assert solver.integrator_settings is None
        assert solver._spectrum_integrated is None

    def test_spectrum_real_packets(self, simulation):
        transport_state = simulation.transport.transport_state
        spectrum_frequency_grid = simulation.transport.spectrum_frequency_grid

        solver = SpectrumSolver(transport_state, spectrum_frequency_grid, None)
        result = solver.spectrum_real_packets.luminosity
        key = "simulation/spectrum_solver/spectrum_real_packets/luminosity"
        expected = self.get_expected_data(key)

        luminosity = u.Quantity(expected, "erg /s")

        assert_quantity_allclose(
            result,
            luminosity,
        )

    def test_spectrum_real_packets_reabsorbed(self, simulation):
        transport_state = simulation.transport.transport_state
        spectrum_frequency_grid = simulation.transport.spectrum_frequency_grid

        solver = SpectrumSolver(transport_state, spectrum_frequency_grid, None)
        result = solver.spectrum_real_packets_reabsorbed.luminosity
        key = "simulation/spectrum_solver/spectrum_real_packets_reabsorbed/luminosity"
        expected = self.get_expected_data(key)

        luminosity = u.Quantity(expected, "erg /s")

        assert_quantity_allclose(
            result,
            luminosity,
        )

    def test_solve(self, simulation):
        transport_state = simulation.transport.transport_state
        spectrum_frequency_grid = simulation.transport.spectrum_frequency_grid

        solver = SpectrumSolver(transport_state, spectrum_frequency_grid, None)
        result_real, result_virtual, result_integrated = solver.solve(
            transport_state
        )
        key_real = "simulation/spectrum_solver/spectrum_real_packets/luminosity"
        expected_real = self.get_expected_data(key_real)

        luminosity_real = u.Quantity(expected_real, "erg /s")

        assert_quantity_allclose(
            result_real.luminosity,
            luminosity_real,
        )

        assert_quantity_allclose(
            result_virtual.luminosity,
            u.Quantity(
                np.zeros_like(spectrum_frequency_grid.value)[:-1], "erg / s"
            ),
        )

        assert result_integrated is None
