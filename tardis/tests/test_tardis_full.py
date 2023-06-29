import os
import pytest
import numpy as np
import numpy.testing as npt
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose

from tardis.simulation.base import Simulation
from tardis.io.config_reader import Configuration

from tardis import run_tardis


def test_run_tardis_from_config_obj(atomic_data_fname):
    """
    Tests whether the run_tardis function can take in the Configuration object
    as arguments
    """
    config = Configuration.from_yaml(
        "tardis/io/tests/data/tardis_configv1_verysimple.yml"
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
    def transport(self, atomic_data_fname, tardis_ref_data, generate_reference):
        config = Configuration.from_yaml(
            "tardis/io/tests/data/tardis_configv1_verysimple.yml"
        )
        config["atom_data"] = atomic_data_fname

        simulation = Simulation.from_config(config)
        simulation.run_convergence()
        simulation.run_final()

        if not generate_reference:
            return simulation.transport
        else:
            simulation.transport.hdf_properties = [
                "j_blue_estimator",
                "spectrum",
                "spectrum_virtual",
            ]
            simulation.transport.to_hdf(
                tardis_ref_data, "", self.name, overwrite=True
            )
            pytest.skip("Reference data was generated during this run.")

    @pytest.fixture(scope="class")
    def refdata(self, tardis_ref_data):
        def get_ref_data(key):
            return tardis_ref_data[os.path.join(self.name, key)]

        return get_ref_data

    def test_j_blue_estimators(self, transport, refdata):
        j_blue_estimator = refdata("j_blue_estimator").values

        npt.assert_allclose(transport.j_blue_estimator, j_blue_estimator)

    def test_spectrum(self, transport, refdata):
        luminosity = u.Quantity(refdata("spectrum/luminosity"), "erg /s")

        assert_quantity_allclose(transport.spectrum.luminosity, luminosity)

    def test_virtual_spectrum(self, transport, refdata):
        luminosity = u.Quantity(
            refdata("spectrum_virtual/luminosity"), "erg /s"
        )

        assert_quantity_allclose(
            transport.spectrum_virtual.luminosity, luminosity
        )

    def test_transport_properties(self, transport):
        """
        Tests whether a number of transport attributes exist and also verifies
        their types

        Currently, transport attributes needed to call the model routine to_hdf5
        are checked.
        """

        virt_type = np.ndarray

        props_required_by_modeltohdf5 = dict(
            [
                ("virt_packet_last_interaction_type", virt_type),
                ("virt_packet_last_line_interaction_in_id", virt_type),
                ("virt_packet_last_line_interaction_in_shell_id", virt_type),
                ("virt_packet_last_line_interaction_out_id", virt_type),
                ("virt_packet_last_line_interaction_out_shell_id", virt_type),
                ("virt_packet_last_interaction_in_nu", virt_type),
                ("virt_packet_nus", virt_type),
                ("virt_packet_energies", virt_type),
            ]
        )

        required_props = props_required_by_modeltohdf5.copy()

        for prop, prop_type in required_props.items():
            actual = getattr(transport, prop)
            assert type(actual) == prop_type, (
                f"wrong type of attribute '{prop}':"
                f"expected {prop_type}, found {type(actual)}"
            )
