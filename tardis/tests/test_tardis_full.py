from pathlib import Path

import pytest
import numpy as np
import numpy.testing as npt
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose

from tardis.simulation.base import Simulation
from tardis.io.configuration.config_reader import Configuration

from tardis import run_tardis


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

    name = "test_transport_simple"

    @pytest.fixture(scope="class")
    def transport(
        self,
        atomic_data_fname,
        tardis_ref_data,
        generate_reference,
        example_configuration_dir: Path,
    ):
        config = Configuration.from_yaml(
            example_configuration_dir / "tardis_configv1_verysimple.yml"
        )
        config["atom_data"] = atomic_data_fname

        simulation = Simulation.from_config(config)
        simulation.run_convergence()
        simulation.run_final()
        if not generate_reference:
            return simulation.transport
        else:
            simulation.transport.hdf_properties = [
                "transport_state",
            ]
            simulation.transport.to_hdf(
                tardis_ref_data, "", self.name, overwrite=True
            )
            pytest.skip("Reference data was generated during this run.")

    @pytest.fixture(scope="class")
    def refdata(self, tardis_ref_data):
        def get_ref_data(key):
            return tardis_ref_data[f"{self.name}/{key}"]

        return get_ref_data

    def test_j_blue_estimators(self, transport, refdata):
        j_blue_estimator = refdata("j_blue_estimator").values

        npt.assert_allclose(
            transport.transport_state.radfield_mc_estimators.j_blue_estimator,
            j_blue_estimator,
        )

    def test_spectrum(self, transport, refdata):
        luminosity = u.Quantity(refdata("transport_state/spectrum/luminosity"), "erg /s")

        assert_quantity_allclose(
            transport.transport_state.spectrum.luminosity, luminosity
        )

    def test_virtual_spectrum(self, transport, refdata):
        luminosity = u.Quantity(
            refdata("transport_state/spectrum_virtual/luminosity"), "erg /s"
        )

        assert_quantity_allclose(
            transport.transport_state.spectrum_virtual.luminosity, luminosity
        )
