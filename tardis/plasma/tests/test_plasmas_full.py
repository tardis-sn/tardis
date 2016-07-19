import os

import h5py
import numpy as np
import numpy.testing as nptesting
import pytest

from astropy import units as u
from astropy.tests.helper import remote_data

import tardis
from tardis.model import Radial1DModel
from tardis.io.config_reader import Configuration
from tardis.simulation.base import run_radial1d


def data_path(fname):
    return os.path.join(tardis.__path__[0], 'plasma', 'tests', 'data', fname)


@pytest.fixture
def plasma_compare_data_fname():
    return data_path('plasma_test_data.h5')


@pytest.fixture
def plasma_compare_data(plasma_compare_data_fname):
    return h5py.File(plasma_compare_data_fname, 'r')


@remote_data
class TestPlasmas(object):
    @classmethod
    @pytest.fixture(scope="class", autouse=True)
    def setup(self, atom_data):
        self.lte_config = Configuration.from_yaml(
            'tardis/plasma/tests/data/plasma_test_config_lte.yml',
            atom_data=atom_data
        )
        self.lte_model = Radial1DModel(self.lte_config)
        run_radial1d(self.lte_model)

        self.nlte_config = Configuration.from_yaml(
            'tardis/plasma/tests/data/plasma_test_config_nlte.yml',
            atom_data=atom_data
        )
        self.nlte_model = Radial1DModel(self.nlte_config)
        run_radial1d(self.nlte_model)

    def test_lte_plasma(self, plasma_compare_data):
        old_plasma_t_rads = plasma_compare_data['test_lte1/t_rad']
        old_plasma_levels = plasma_compare_data['test_lte1/levels']

        new_plasma_t_rads = self.lte_model.t_rads / u.Unit('K')
        new_plasma_levels = \
            self.lte_model.plasma.get_value(
            'level_number_density').ix[8].ix[1][10].values
        np.testing.assert_allclose(
            new_plasma_t_rads, old_plasma_t_rads, atol=100)
        np.testing.assert_allclose(
            new_plasma_levels, old_plasma_levels, rtol=0.1)

    def test_nlte_plasma(self, plasma_compare_data):
        old_plasma_t_rads = plasma_compare_data['test_nlte1/t_rad']
        old_plasma_levels = plasma_compare_data['test_nlte1/levels']
        new_plasma_t_rads = self.nlte_model.t_rads / u.Unit('K')
        new_plasma_levels = \
            self.nlte_model.plasma.get_value(
            'level_number_density').ix[2].ix[1][10].values
        np.testing.assert_allclose(
            new_plasma_t_rads, old_plasma_t_rads, atol=150)
        np.testing.assert_allclose(
            new_plasma_levels, old_plasma_levels, rtol=0.1)
