from pathlib import Path

import numpy.testing as npt
import pandas as pd
import pytest
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose

from tardis import run_tardis
from tardis.io.configuration.config_reader import Configuration
from tardis.simulation.base import Simulation
from tardis.tests.fixtures.regression_data import RegressionData


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


class TestTransportSimple:
    """
    Very simple run
    """

    regression_data: RegressionData = None

    @pytest.fixture(scope="class")
    def transport_state(
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

        transport_state = simulation.transport.transport_state
        request.cls.regression_data = RegressionData(request)
        request.cls.regression_data.sync_hdf_store(transport_state)

        return transport_state

    def get_expected_data(self, key: str):
        return pd.read_hdf(self.regression_data.fpath, key)

    def test_j_blue_estimators(self, transport_state):
        key = "transport_state/j_blue_estimator"
        expected = self.get_expected_data(key)

        npt.assert_allclose(
            transport_state.radfield_mc_estimators.j_blue_estimator,
            expected.values,
        )

    def test_spectrum(self, transport_state):
        key = "transport_state/spectrum/luminosity"
        expected = self.get_expected_data(key)

        luminosity = u.Quantity(expected, "erg /s")

        assert_quantity_allclose(
            transport_state.spectrum.luminosity, luminosity
        )

    def test_virtual_spectrum(self, transport_state):
        key = "transport_state/spectrum_virtual/luminosity"
        expected = self.get_expected_data(key)

        luminosity = u.Quantity(expected, "erg /s")

        assert_quantity_allclose(
            transport_state.spectrum_virtual.luminosity, luminosity
        )
