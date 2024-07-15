from pathlib import Path

import numpy.testing as npt
import pytest
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose

from tardis.io.configuration.config_reader import Configuration
from tardis.simulation.base import Simulation

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


class TestTransportSimpleFormalIntegral:
    """
    Very simple run with the formal integral spectral synthesis method
    """

    _name = "test_transport_simple_integral"

    @pytest.fixture(scope="class")
    def transport(
        self, config, atomic_data_fname, tardis_ref_data, generate_reference
    ):
        config.atom_data = atomic_data_fname

        self.name = self._name + f"_{config.plasma.line_interaction_type:s}"
        if config.spectrum.integrated.interpolate_shells > 0:
            self.name += "_interp"

        simulation = Simulation.from_config(config)
        simulation.run_convergence()
        simulation.run_final()

        if not generate_reference:
            return simulation.transport
        else:
            simulation.transport.hdf_properties = ["transport_state"]
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
        j_blue_estimator = refdata("transport_state/j_blue_estimator").values

        npt.assert_allclose(
            transport.transport_state.radfield_mc_estimators.j_blue_estimator,
            j_blue_estimator,
        )

    def test_spectrum(self, transport, refdata):
        luminosity = u.Quantity(
            refdata("transport_state/spectrum/luminosity"), "erg /s"
        )

        assert_quantity_allclose(
            transport.transport_state.spectrum.luminosity, luminosity
        )

    def test_spectrum_integrated(self, transport, refdata):
        luminosity = u.Quantity(
            refdata("transport_state/spectrum_integrated/luminosity"), "erg /s"
        )

        assert_quantity_allclose(
            transport.transport_state.spectrum_integrated.luminosity, luminosity
        )
