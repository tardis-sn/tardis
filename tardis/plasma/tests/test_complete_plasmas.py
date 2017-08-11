import os

import pandas as pd
import pytest
from pandas.util import testing as pdt
from numpy.testing import assert_almost_equal
from tardis.io.config_reader import Configuration
from tardis.simulation import Simulation


class BasePlasmaTest(object):
    #Class defining all common tests for different setups of Plasma
    #This can then be inherited for different Plasma setup
    @classmethod
    @pytest.fixture(scope="class", autouse=True)
    def setup(cls, atomic_data_fname, tardis_ref_path, reference_file_path, config_path):
        cls.config = Configuration.from_yaml(config_path)
        cls.config['atom_data'] = atomic_data_fname
        cls.sim = Simulation.from_config(cls.config)

        if pytest.config.getvalue("--generate-reference"):
            if os.path.exists(reference_file_path):
                pytest.skip(
                    'Reference data {0} does exist and tests will not '
                    'proceed generating new data'.format(reference_file_path))
            cls.sim.plasma.to_hdf(reference_file_path)
            pytest.skip("Reference data saved at {0}".format(
                reference_file_path))
        cls.plasma = cls.sim.plasma

    @pytest.fixture(scope="class")
    def reference_file_path(self):
        pass

    @pytest.fixture(scope="class")
    def config_path(self):
        pass

    general_properties = ['beta_rad', 'g_electron', 'selected_atoms',
                          'number_density', 't_electrons', 'w', 't_rad']
    partiton_properties = ['level_boltzmann_factor', 'partition_function']
    atomic_properties = ['excitation_energy', 'lines', 'lines_lower_level_index',
                         'lines_upper_level_index', 'atomic_mass', 'ionization_data',
                         'nu', 'wavelength_cm', 'f_lu', 'metastability']
    ion_population_properties = [
        'phi', 'ion_number_density', 'electron_densities']
    level_population_properties = ['level_number_density']
    radiative_properties = ['stimulated_emission_factor',
                            'tau_sobolevs', 'beta_sobolev', 'transition_probabilities']
    j_blues_properties = ['j_blues']

    combined_properties = general_properties + partiton_properties + atomic_properties + ion_population_properties + \
        level_population_properties + radiative_properties + \
        j_blues_properties

    @pytest.mark.parametrize("attr", combined_properties)
    def test_plasma_properties(self, reference_file_path, attr):
        actual = getattr(self.plasma, attr)
        expected = pd.read_hdf(reference_file_path,
                               os.path.join('plasma', attr))
        assert_almost_equal(actual, expected.values)

    def test_levels(self, reference_file_path):
        actual = pd.DataFrame(self.plasma.levels)
        expected = pd.read_hdf(reference_file_path,
                               os.path.join('plasma', 'levels'))
        pdt.assert_almost_equal(actual, expected)

    scalars_properties = ['time_explosion', 'link_t_rad_t_electron']

    @pytest.mark.parametrize("attr", scalars_properties)
    def test_scalars_properties(self, reference_file_path, attr):
        actual = getattr(self.plasma, attr)
        if hasattr(actual, 'cgs'):
            actual = actual.cgs.value
        expected = pd.read_hdf(reference_file_path, os.path.join(
            'plasma', 'scalars'))[attr]
        assert_almost_equal(actual, expected)

    def test_helium_treatment(self, reference_file_path):
        actual = self.plasma.helium_treatment
        expected = pd.read_hdf(reference_file_path, os.path.join(
            'plasma', 'scalars'))['helium_treatment']
        assert actual == expected


class TestLTEPlasma(BasePlasmaTest):

    @pytest.fixture(scope="class")
    def reference_file_path(self, tardis_ref_path):
        return os.path.join(tardis_ref_path, 'plasma_reference', 'plasma_lte_reference.h5')

    @pytest.fixture(scope="class")
    def config_path(self):
        return os.path.join('tardis', 'plasma', 'tests', 'data', 'plasma_test_config_lte.yml')


class TestNLTEPlasma(BasePlasmaTest):

    @pytest.fixture(scope="class")
    def reference_file_path(self, tardis_ref_path):
        return os.path.join(tardis_ref_path, 'plasma_reference', 'plasma_nlte_reference.h5')

    @pytest.fixture(scope="class")
    def config_path(self):
        return os.path.join('tardis', 'plasma', 'tests', 'data', 'plasma_test_config_nlte.yml')

    #Additional Tests for NLTE Plasma

    nlte_general_properties = ['beta_electron']
    nlte_atomic_properties = ['zeta_data']
    nlte_ion_population_properties = ['delta', 'previous_electron_densities']
    nlte_radiative_properties = ['previous_beta_sobolev']
    nlte_properties = nlte_atomic_properties + nlte_general_properties + \
        nlte_ion_population_properties + nlte_radiative_properties

    @pytest.mark.parametrize("attr", nlte_properties)
    def test_nlte_properties(self, reference_file_path, attr):
        actual = getattr(self.plasma, attr)
        expected = pd.read_hdf(reference_file_path,
                               os.path.join('plasma', attr))
        assert_almost_equal(actual, expected.values)
