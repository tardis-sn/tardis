from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
import pytest
from astropy.tests.helper import assert_quantity_allclose
from tardisbase.testing.regression_data.regression_data import RegressionData

from tardis.io.configuration.config_reader import Configuration
from tardis.simulation.base import Simulation
from tardis.spectrum.base import SpectrumSolver
from tardis.spectrum.formal_integral.base import IntegrationError
from tardis.spectrum.formal_integral.formal_integral_solver import (
    FormalIntegralSolver,
)


class TestSpectrumSolver:
    regression_data: RegressionData = None

    @pytest.fixture(scope="class")
    def simulation_regression_data(
        self,
        simulation_verysimple_default,
        request: pytest.FixtureRequest,
    ):
        request.cls.regression_data = RegressionData(request)
        data = request.cls.regression_data.sync_hdf_store(simulation_verysimple_default)

        yield simulation_verysimple_default
        data.close()

    def get_expected_data(self, key: str):
        return pd.read_hdf(self.regression_data.fpath, key)

    def test_initialization(self, simulation_regression_data):
        transport_state = simulation_regression_data.transport.transport_state
        spectrum_frequency_grid = simulation_regression_data.transport.spectrum_frequency_grid

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

    def test_spectrum_real_packets(self, simulation_regression_data):
        transport_state = simulation_regression_data.transport.transport_state
        spectrum_frequency_grid = simulation_regression_data.transport.spectrum_frequency_grid

        solver = SpectrumSolver(transport_state, spectrum_frequency_grid, None)
        result = solver.spectrum_real_packets.luminosity
        key = "simulation/spectrum_solver/spectrum_real_packets/luminosity"
        expected = self.get_expected_data(key)

        luminosity = u.Quantity(expected, "erg /s")

        assert_quantity_allclose(
            result,
            luminosity,
        )

    def test_spectrum_real_packets_reabsorbed(self, simulation_regression_data):
        transport_state = simulation_regression_data.transport.transport_state
        spectrum_frequency_grid = simulation_regression_data.transport.spectrum_frequency_grid

        solver = SpectrumSolver(transport_state, spectrum_frequency_grid, None)
        result = solver.spectrum_real_packets_reabsorbed.luminosity
        key = "simulation/spectrum_solver/spectrum_real_packets_reabsorbed/luminosity"
        expected = self.get_expected_data(key)

        luminosity = u.Quantity(expected, "erg /s")

        assert_quantity_allclose(
            result,
            luminosity,
        )

    def test_solve(self, simulation_regression_data):
        transport_state = simulation_regression_data.transport.transport_state
        spectrum_frequency_grid = simulation_regression_data.transport.spectrum_frequency_grid

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


def test_nonuniform_frequency_grid_raises_integration_error(
    simulation_verysimple: Simulation,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The real formal solver reports unsupported frequency spacing."""
    simulation = simulation_verysimple
    monkeypatch.setattr(simulation.transport, "line_interaction_type", "macroatom")
    frequencies = np.array([4.0, 5.0, 7.0]) * 1e14 * u.Hz
    integrator = FormalIntegralSolver(points=10, interpolate_shells=-1, method="numba")

    with pytest.raises(IntegrationError, match="requires a uniform frequency grid"):
        integrator.solve(
            frequencies,
            simulation.simulation_state,
            simulation.transport,
            simulation.opacity_state,
            simulation.plasma.atomic_data,
            simulation.plasma.electron_densities,
            simulation.macro_atom_state,
        )


def test_integrated_spectrum_warns_with_integration_error(
    simulation_verysimple: Simulation,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A real nonuniform integration produces a diagnostic and NaN fallback."""
    simulation = simulation_verysimple
    monkeypatch.setattr(simulation.transport, "line_interaction_type", "macroatom")
    transport_state = simulation.transport.transport_state
    frequencies = np.array([4.0, 5.0, 7.0, 8.0]) * 1e14 * u.Hz
    solver = SpectrumSolver(transport_state, frequencies, None)
    integrator = FormalIntegralSolver(points=10, interpolate_shells=-1, method="numba")
    solver.setup_optional_spectra(
        transport_state,
        integrator=integrator,
        simulation_state=simulation.simulation_state,
        transport=simulation.transport,
        plasma=simulation.plasma,
        opacity_state=simulation.opacity_state,
        macro_atom_state=simulation.macro_atom_state,
    )

    with pytest.warns(UserWarning, match="requires a uniform frequency grid"):
        spectrum = solver.spectrum_integrated

    assert np.isnan(spectrum.frequency.value).all()
    assert np.isnan(spectrum.luminosity.value).all()
