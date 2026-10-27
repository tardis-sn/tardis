from collections.abc import Callable
from types import SimpleNamespace
from unittest.mock import Mock

import astropy.units as u
import numpy as np
import pandas as pd
import pytest
from astropy.tests.helper import assert_quantity_allclose
from tardisbase.testing.regression_data.regression_data import RegressionData

from tardis.spectrum.base import SpectrumSolver
from tardis.transport.montecarlo.configuration.base import (
    MonteCarloConfiguration,
)
from tardis.transport.montecarlo.modes.nonhomologous.montecarlo_transport_state import (
    MonteCarloTransportStateNonhomologous,
)
from tardis.transport.montecarlo.modes.nonhomologous.virtual_packet import (
    trace_vpacket_volley as trace_vpacket_volley_nonhomologous,
)
from tardis.transport.montecarlo.montecarlo_transport_state import (
    MonteCarloTransportState,
)
from tardis.transport.montecarlo.packets.virtual_packet import (
    trace_vpacket_volley as trace_vpacket_volley_classic,
)


@pytest.mark.parametrize(
    (
        "transport_state_type",
        "expected_trace_function",
        "expected_time_explosion",
    ),
    [
        (MonteCarloTransportState, trace_vpacket_volley_classic, 123.0),
        (
            MonteCarloTransportStateNonhomologous,
            trace_vpacket_volley_nonhomologous,
            0.0,
        ),
    ],
)
def test_virtual_packet_postprocessing_owned_by_spectrum_solver(
    monkeypatch: pytest.MonkeyPatch,
    transport_state_type: (
        type[MonteCarloTransportState]
        | type[MonteCarloTransportStateNonhomologous]
    ),
    expected_trace_function: Callable[..., None],
    expected_time_explosion: float,
) -> None:
    spectrum_frequency_grid = np.arange(4) * u.Hz
    virtual_packet_luminosity = np.arange(4, dtype=float)
    rpacket_tracker = SimpleNamespace()
    vpacket_tracker = SimpleNamespace()
    transport_state = object.__new__(transport_state_type)
    transport_state.geometry_state_numba = SimpleNamespace()
    transport_state.opacity_state_numba = SimpleNamespace()
    transport_state.time_explosion = 123.0 * u.s
    montecarlo_configuration = MonteCarloConfiguration()
    montecarlo_configuration.NUMBER_OF_VPACKETS = 3
    montecarlo_configuration.ENABLE_VPACKET_TRACKING = True
    calculate_virtual_packet_spectrum = Mock(
        return_value=(virtual_packet_luminosity, vpacket_tracker)
    )

    monkeypatch.setattr(
        "tardis.spectrum.base.calculate_virtual_packet_spectrum",
        calculate_virtual_packet_spectrum,
    )

    solver = SpectrumSolver(
        transport_state, spectrum_frequency_grid, integrator_settings=None
    )
    solver.solve_virtual_packet_spectrum(
        transport_state,
        montecarlo_configuration,
        rpacket_tracker,
    )

    calculate_arguments = calculate_virtual_packet_spectrum.call_args.args
    assert calculate_arguments[1] == expected_time_explosion
    assert calculate_arguments[5] is rpacket_tracker
    assert calculate_arguments[7] is expected_trace_function
    np.testing.assert_array_equal(
        solver._montecarlo_virtual_luminosity.value,
        virtual_packet_luminosity,
    )
    assert transport_state.vpacket_tracker is vpacket_tracker
    assert transport_state.virt_logging is True


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
        spectrum_frequency_grid = (
            simulation_regression_data.spectrum_solver.spectrum_frequency_grid
        )

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
        spectrum_frequency_grid = (
            simulation_regression_data.spectrum_solver.spectrum_frequency_grid
        )

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
        spectrum_frequency_grid = (
            simulation_regression_data.spectrum_solver.spectrum_frequency_grid
        )

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
        spectrum_frequency_grid = (
            simulation_regression_data.spectrum_solver.spectrum_frequency_grid
        )

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
