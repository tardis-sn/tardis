import os
import pytest
import numpy as np
import numpy.testing as npt
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose

from tardis.simulation.base import Simulation
from tardis.io.config_reader import Configuration


class TestRunnerSimpleFormalInegral():
    """
    Very simple run with the formal integral spectral synthesis method
    """
    name = 'test_runner_simple_integral'

    @pytest.fixture(scope="class")
    def runner(
            self, atomic_data_fname,
            tardis_ref_data, generate_reference):
        config = Configuration.from_yaml(
                'tardis/io/tests/data/tardis_configv1_verysimple.yml')
        config['atom_data'] = atomic_data_fname
        config["plasma"]["line_interaction_type"] = "downbranch"
        config["montecarlo"]["no_of_virtual_packets"] = 0
        config["spectrum"]["method"] = "integrated"

        simulation = Simulation.from_config(config)
        simulation.run()

        if not generate_reference:
            return simulation.runner
        else:
            simulation.runner.hdf_properties = [
                    'j_blue_estimator',
                    'spectrum',
                    'spectrum_integrated'
                    ]
            simulation.runner.to_hdf(
                    tardis_ref_data,
                    '',
                    self.name)
            pytest.skip(
                    'Reference data was generated during this run.')

    @pytest.fixture(scope='class')
    def refdata(self, tardis_ref_data):
        def get_ref_data(key):
            return tardis_ref_data[os.path.join(
                    self.name, key)]
        return get_ref_data

    def test_j_blue_estimators(self, runner, refdata):
        j_blue_estimator = refdata('j_blue_estimator').values

        npt.assert_allclose(
                runner.j_blue_estimator,
                j_blue_estimator)

    def test_spectrum(self, runner, refdata):
        luminosity = u.Quantity(refdata('spectrum/luminosity'), 'erg /s')

        assert_quantity_allclose(
            runner.spectrum.luminosity,
            luminosity)

    def test_virtual_integrated(self, runner, refdata):
        luminosity = u.Quantity(
                refdata('spectrum_virtual/luminosity'), 'erg /s')

        assert_quantity_allclose(
            runner.spectrum_integrated.luminosity,
            luminosity)
