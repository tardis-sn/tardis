from pathlib import Path

import numpy.testing as npt
import pandas as pd
import pytest
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose
from tardisbase.testing.regression_data.regression_data import RegressionData

from tardis import run_tardis
from tardis.io.configuration.config_reader import Configuration
from tardis.simulation.base import Simulation


def test_run_tardis_from_config_obj(
    atomic_data_fname, example_configuration_dir: Path
):
    """
    Tests whether the run_tardis function can take in the Configuration object
    as arguments
    """
    config = Configuration.from_yaml(
        example_configuration_dir / "tardis_configv1_verysimple.yml"
    )
    config["atom_data"] = atomic_data_fname

    try:
        sim = run_tardis(config)
    except Exception as e:
        pytest.fail(str(e.args[0]))


def test_run_tardis_simulation_callbacks_none(
    atomic_data_fname, example_configuration_dir: Path
):
    """
    Test that run_tardis handles simulation_callbacks=None correctly
    """
    config = Configuration.from_yaml(
        example_configuration_dir / "tardis_configv1_verysimple.yml"
    )
    config["atom_data"] = atomic_data_fname

    try:
        sim = run_tardis(config, simulation_callbacks=None)
        assert isinstance(sim, Simulation)
    except Exception as e:
        pytest.fail(f"run_tardis failed with simulation_callbacks=None: {e}")


class TestTransportSimple:
    """
    Very simple run
    """

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
        simulation.run_convergence()
        simulation.run_final()

        request.cls.regression_data = RegressionData(request)
        data = request.cls.regression_data.sync_hdf_store(simulation)

        yield simulation
        data.close()

    def get_expected_data(self, key: str):
        return pd.read_hdf(self.regression_data.fpath, key)

    def test_j_blue_estimators(self, simulation):
        key = "simulation/transport/transport_state/j_blue_estimator"
        expected = self.get_expected_data(key)

        npt.assert_allclose(
            simulation.transport.transport_state.estimators_line.mean_intensity_blueward,
            expected.values,
        )

    def test_spectrum(self, simulation):
        key = "simulation/spectrum_solver/spectrum_real_packets/luminosity"
        expected = self.get_expected_data(key)

        luminosity = u.Quantity(expected, "erg /s")

        assert_quantity_allclose(
            simulation.spectrum_solver.spectrum_real_packets.luminosity,
            luminosity,
        )

    def test_virtual_spectrum(self, simulation):
        key = "simulation/spectrum_solver/spectrum_virtual_packets/luminosity"
        expected = self.get_expected_data(key)

        luminosity = u.Quantity(expected, "erg /s")

        assert_quantity_allclose(
            simulation.spectrum_solver.spectrum_virtual_packets.luminosity,
            luminosity,
        )
