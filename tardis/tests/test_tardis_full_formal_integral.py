import os
import pytest
import numpy.testing as npt
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose

from tardis.simulation.base import Simulation
from tardis.io.config_reader import Configuration

pytestmark = pytest.mark.skip(reason="memory problem")

config_line_modes = ["downbranch", "macroatom"]
interpolate_shells = [-1, 30]


@pytest.fixture(scope="module", params=config_line_modes)
def base_config(request):
    config = Configuration.from_yaml(
        "tardis/io/tests/data/tardis_configv1_verysimple.yml"
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


class TestRunnerSimpleFormalInegral:
    """
    Very simple run with the formal integral spectral synthesis method
    """

    _name = "test_runner_simple_integral"

    @pytest.fixture(scope="class")
    def runner(
        self, config, atomic_data_fname, tardis_ref_data, generate_reference
    ):
        config.atom_data = atomic_data_fname

        self.name = self._name + "_{:s}".format(
            config.plasma.line_interaction_type
        )
        if config.spectrum.integrated.interpolate_shells > 0:
            self.name += "_interp"

        simulation = Simulation.from_config(config)
        simulation.run()

        if not generate_reference:
            return simulation.runner
        else:
            simulation.runner.hdf_properties = [
                "j_blue_estimator",
                "spectrum",
                "spectrum_integrated",
            ]
            simulation.runner.to_hdf(tardis_ref_data, "", self.name)
            pytest.skip("Reference data was generated during this run.")

    @pytest.fixture(scope="class")
    def refdata(self, tardis_ref_data):
        def get_ref_data(key):
            return tardis_ref_data[os.path.join(self.name, key)]

        return get_ref_data

    def test_j_blue_estimators(self, runner, refdata):
        j_blue_estimator = refdata("j_blue_estimator").values

        npt.assert_allclose(runner.j_blue_estimator, j_blue_estimator)

    def test_spectrum(self, runner, refdata):
        luminosity = u.Quantity(refdata("spectrum/luminosity"), "erg /s")

        assert_quantity_allclose(runner.spectrum.luminosity, luminosity)

    def test_spectrum_integrated(self, runner, refdata):
        luminosity = u.Quantity(
            refdata("spectrum_integrated/luminosity"), "erg /s"
        )

        print(
            "actual, desired: ",
            luminosity,
            runner.spectrum_integrated.luminosity,
        )
        assert_quantity_allclose(
            runner.spectrum_integrated.luminosity, luminosity
        )
